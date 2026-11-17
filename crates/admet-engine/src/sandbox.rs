use std::fs::File;
use std::io::Read;
use std::path::{Component, Path, PathBuf};
use std::sync::Arc;

#[cfg(unix)]
use std::os::fd::OwnedFd;

#[cfg(unix)]
use rustix::fs::{AtFlags, FileType, Mode, OFlags, ResolveFlags, mkdirat, openat, openat2, statat};
use thiserror::Error;

/// Default maximum bytes read from a single sandboxed input file (10 MiB).
pub const MAX_READ_BYTES: usize = 10 * 1024 * 1024;

/// Errors raised when validating or opening sandboxed filesystem paths.
#[derive(Debug, Error)]
pub enum SandboxError {
    #[error("path '{0}' is outside authorized root '{1}'")]
    OutsideRoot(String, String),
    #[error("symbolic links are forbidden inside sandbox: '{0}'")]
    SymlinkNotAllowed(String),
    #[error("path '{0}' is not a regular file (FIFOs, sockets, and directories are rejected)")]
    NotARegularFile(String),
    #[error("file '{0}' exceeds maximum allowed size of {1} bytes")]
    FileTooLarge(String, usize),
    #[error("no root directory configured for this operation")]
    NoRootConfigured,
    #[error("I/O error resolving path: {0}")]
    Io(#[from] std::io::Error),
}

/// Directory sandbox enforcing strict `--read-root` and `--write-root` isolation using
/// pinned root directory file descriptors and kernel `openat2(RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS)`
/// so intermediate-directory symlink swaps cannot escape the sandbox.
#[derive(Debug, Clone)]
pub struct PathSandbox {
    read_root: Option<PathBuf>,
    write_root: Option<PathBuf>,
    #[cfg(unix)]
    read_fd: Option<Arc<OwnedFd>>,
    #[cfg(unix)]
    write_fd: Option<Arc<OwnedFd>>,
}

impl PathSandbox {
    /// Create a sandbox with optional canonicalized read and write roots.
    pub fn new(read_root: Option<&Path>, write_root: Option<&Path>) -> Result<Self, SandboxError> {
        let r = match read_root {
            Some(p) => Some(std::fs::canonicalize(p)?),
            None => None,
        };
        let w = match write_root {
            Some(p) => {
                if !p.exists() {
                    std::fs::create_dir_all(p)?;
                }
                Some(std::fs::canonicalize(p)?)
            }
            None => None,
        };

        #[cfg(unix)]
        let read_fd = match &r {
            Some(p) => Some(Arc::new(open_root_dir_fd(p)?)),
            None => None,
        };
        #[cfg(unix)]
        let write_fd = match &w {
            Some(p) => Some(Arc::new(open_root_dir_fd(p)?)),
            None => None,
        };

        Ok(Self {
            read_root: r,
            write_root: w,
            #[cfg(unix)]
            read_fd,
            #[cfg(unix)]
            write_fd,
        })
    }

    /// Resolve and verify read access (note: prefer `open_read` / `read_to_string_limited`
    /// to operate directly on the kernel-opened file descriptor).
    pub fn resolve_read(&self, candidate: &Path) -> Result<PathBuf, SandboxError> {
        let (path, _file) = self.open_read(candidate)?;
        Ok(path)
    }

    /// Securely open a regular file for reading inside `read_root` atomically via `openat2`
    /// (`RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS`), rejecting symlinks at any component, FIFOs,
    /// sockets, directories, and path escapes.
    pub fn open_read(&self, candidate: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self
            .read_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let rel = normalize_relative_under_root(root, candidate)?;
        let target = root.join(&rel);

        #[cfg(unix)]
        let file = {
            let root_fd = self
                .read_fd
                .as_ref()
                .ok_or(SandboxError::NoRootConfigured)?;
            let fd = open_beneath_no_symlinks(
                root_fd.as_ref(),
                &rel,
                OFlags::RDONLY | OFlags::NOFOLLOW | OFlags::NONBLOCK | OFlags::CLOEXEC,
                Mode::empty(),
            )
            .map_err(|e| map_unix_open_err(e, &target))?;
            File::from(fd)
        };
        #[cfg(not(unix))]
        let file = File::open(&target)?;

        let meta = file.metadata()?;
        if !meta.is_file() {
            return Err(SandboxError::NotARegularFile(target.display().to_string()));
        }
        Ok((target, file))
    }

    /// Securely read a bounded UTF-8 file inside `read_root` (protects against FIFOs and OOM).
    pub fn read_to_string_limited(
        &self,
        candidate: &Path,
        max_bytes: usize,
    ) -> Result<(PathBuf, String), SandboxError> {
        let (path, file) = self.open_read(candidate)?;
        let mut buf = Vec::new();
        file.take((max_bytes + 1) as u64).read_to_end(&mut buf)?;
        if buf.len() > max_bytes {
            return Err(SandboxError::FileTooLarge(
                path.display().to_string(),
                max_bytes,
            ));
        }
        let s = String::from_utf8(buf)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        Ok((path, s))
    }

    /// Validate a relative write target path lexically under `write_root`.
    pub fn resolve_write(&self, candidate: &Path) -> Result<PathBuf, SandboxError> {
        let root = self
            .write_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let rel = normalize_relative_under_root(root, candidate)?;
        Ok(root.join(rel))
    }

    /// Atomically open/create a writable regular file inside `write_root` using `dirfd`-pinned
    /// parent traversal (`openat2` + `mkdirat`) and `RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS`,
    /// truncating only after confirming the opened fd is a regular file.
    pub fn open_write(&self, candidate: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self
            .write_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let rel = normalize_relative_under_root(root, candidate)?;
        let target = root.join(&rel);

        #[cfg(unix)]
        let file = {
            let root_fd = self
                .write_fd
                .as_ref()
                .ok_or(SandboxError::NoRootConfigured)?;
            let comps: Vec<&std::ffi::OsStr> = rel.iter().collect();
            let Some((file_name, parent_comps)) = comps.split_last() else {
                return Err(SandboxError::OutsideRoot(
                    candidate.display().to_string(),
                    root.display().to_string(),
                ));
            };

            let mut parent_fd: Option<OwnedFd> = None;
            for comp in parent_comps {
                let base_fd = parent_fd.as_ref().map_or_else(|| root_fd.as_ref(), |f| f);
                let next_fd = open_or_create_subdir(base_fd, Path::new(comp), &target)?;
                parent_fd = Some(next_fd);
            }
            let dir_ref = parent_fd.as_ref().map_or_else(|| root_fd.as_ref(), |f| f);

            let fd = open_beneath_no_symlinks(
                dir_ref,
                Path::new(file_name),
                OFlags::WRONLY
                    | OFlags::CREATE
                    | OFlags::NOFOLLOW
                    | OFlags::NONBLOCK
                    | OFlags::CLOEXEC,
                Mode::from_bits_truncate(0o644),
            )
            .map_err(|e| map_unix_open_err(e, &target))?;
            File::from(fd)
        };
        #[cfg(not(unix))]
        let file = {
            if let Some(parent) = target.parent()
                && !parent.exists()
            {
                std::fs::create_dir_all(parent)?;
            }
            std::fs::OpenOptions::new()
                .write(true)
                .create(true)
                .open(&target)?
        };

        let meta = file.metadata()?;
        if !meta.is_file() {
            return Err(SandboxError::NotARegularFile(target.display().to_string()));
        }
        file.set_len(0)?;
        Ok((target, file))
    }
}

#[cfg(unix)]
fn open_root_dir_fd(root: &Path) -> Result<OwnedFd, SandboxError> {
    let f = std::os::unix::fs::OpenOptionsExt::custom_flags(
        std::fs::OpenOptions::new().read(true),
        libc::O_DIRECTORY | libc::O_NOFOLLOW | libc::O_CLOEXEC,
    )
    .open(root)?;
    Ok(OwnedFd::from(f))
}

#[cfg(unix)]
fn open_beneath_no_symlinks(
    dirfd: &OwnedFd,
    rel: &Path,
    flags: OFlags,
    mode: Mode,
) -> Result<OwnedFd, rustix::io::Errno> {
    match openat2(
        dirfd,
        rel,
        flags,
        mode,
        ResolveFlags::BENEATH | ResolveFlags::NO_SYMLINKS,
    ) {
        Ok(fd) => Ok(fd),
        Err(rustix::io::Errno::NOSYS) => {
            // Component-by-component dirfd walk fallback on pre-5.6 kernels
            let comps: Vec<&std::ffi::OsStr> = rel.iter().collect();
            let Some((last, parents)) = comps.split_last() else {
                return Err(rustix::io::Errno::INVAL);
            };
            let mut cur: Option<OwnedFd> = None;
            for comp in parents {
                let base = cur.as_ref().unwrap_or(dirfd);
                let next = openat(
                    base,
                    Path::new(comp),
                    OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
                    Mode::empty(),
                )?;
                cur = Some(next);
            }
            let base = cur.as_ref().unwrap_or(dirfd);
            openat(base, Path::new(last), flags | OFlags::NOFOLLOW, mode)
        }
        Err(e) => Err(e),
    }
}

#[cfg(unix)]
fn open_or_create_subdir(
    dirfd: &OwnedFd,
    comp: &Path,
    full_target: &Path,
) -> Result<OwnedFd, SandboxError> {
    let dir_flags = OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC;
    match open_beneath_no_symlinks(dirfd, comp, dir_flags, Mode::empty()) {
        Ok(fd) => Ok(fd),
        Err(rustix::io::Errno::NOENT) => {
            match mkdirat(dirfd, comp, Mode::from_bits_truncate(0o755)) {
                Ok(()) | Err(rustix::io::Errno::EXIST) => {}
                Err(e) => return Err(map_unix_open_err(e, full_target)),
            }
            open_beneath_no_symlinks(dirfd, comp, dir_flags, Mode::empty())
                .map_err(|e| map_unix_open_err(e, full_target))
        }
        Err(e) => {
            if e == rustix::io::Errno::NOTDIR
                && let Ok(st) = statat(dirfd, comp, AtFlags::SYMLINK_NOFOLLOW)
                && FileType::from_raw_mode(st.st_mode).is_symlink()
            {
                return Err(SandboxError::SymlinkNotAllowed(
                    full_target.display().to_string(),
                ));
            }
            Err(map_unix_open_err(e, full_target))
        }
    }
}

#[cfg(unix)]
fn map_unix_open_err(err: rustix::io::Errno, target: &Path) -> SandboxError {
    if err == rustix::io::Errno::LOOP {
        SandboxError::SymlinkNotAllowed(target.display().to_string())
    } else if err == rustix::io::Errno::XDEV {
        SandboxError::OutsideRoot(target.display().to_string(), "sandbox root".to_string())
    } else if err == rustix::io::Errno::NOTDIR
        || err == rustix::io::Errno::ISDIR
        || err == rustix::io::Errno::NXIO
    {
        SandboxError::NotARegularFile(target.display().to_string())
    } else {
        SandboxError::Io(std::io::Error::from_raw_os_error(err.raw_os_error()))
    }
}

fn normalize_relative_under_root(root: &Path, candidate: &Path) -> Result<PathBuf, SandboxError> {
    let joined = if candidate.is_absolute() {
        normalize_lexical(candidate)
    } else {
        normalize_lexical(&root.join(candidate))
    };
    let rel = joined
        .strip_prefix(root)
        .map_err(|_| {
            SandboxError::OutsideRoot(candidate.display().to_string(), root.display().to_string())
        })?
        .to_path_buf();
    if rel.as_os_str().is_empty() {
        return Err(SandboxError::OutsideRoot(
            candidate.display().to_string(),
            root.display().to_string(),
        ));
    }
    Ok(rel)
}

fn normalize_lexical(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for comp in path.components() {
        match comp {
            Component::CurDir => {}
            Component::ParentDir => {
                out.pop();
            }
            other => out.push(other.as_os_str()),
        }
    }
    out
}
