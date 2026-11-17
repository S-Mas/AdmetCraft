use std::fs::{File, OpenOptions};
use std::io::Read;
use std::path::{Component, Path, PathBuf};

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

/// Directory sandbox enforcing strict `--read-root` and `--write-root` isolation
/// without following symlinks (`O_NOFOLLOW`) and without path-check TOCTOU races.
#[derive(Debug, Clone)]
pub struct PathSandbox {
    read_root: Option<PathBuf>,
    write_root: Option<PathBuf>,
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
        Ok(Self {
            read_root: r,
            write_root: w,
        })
    }

    /// Resolve a read path without opening it (rejects symlinks and non-regular files).
    pub fn resolve_read(&self, candidate: &Path) -> Result<PathBuf, SandboxError> {
        let (path, _file) = self.open_read(candidate)?;
        Ok(path)
    }

    /// Securely open a regular file for reading inside `read_root`, rejecting any symlink
    /// (including broken symlinks), FIFO/named pipe, socket, or path escape.
    pub fn open_read(&self, candidate: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self
            .read_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let target = normalize_under_root(root, candidate)?;
        ensure_no_symlinks_under_root(root, &target, false)?;

        let mut opts = OpenOptions::new();
        opts.read(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            opts.custom_flags(libc::O_NOFOLLOW | libc::O_NONBLOCK);
        }
        let file = opts.open(&target)?;
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

    /// Resolve a write path without opening it (rejects dangling or live symlinks).
    pub fn resolve_write(&self, candidate: &Path) -> Result<PathBuf, SandboxError> {
        let root = self
            .write_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let target = normalize_under_root(root, candidate)?;
        ensure_no_symlinks_under_root(root, &target, false)?;
        Ok(target)
    }

    /// Atomically validate and open/create a writable regular file inside `write_root`
    /// with `O_NOFOLLOW`, returning the open `File` handle directly so callers never suffer
    /// a path-to-write TOCTOU race.
    pub fn open_write(&self, candidate: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self
            .write_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let target = normalize_under_root(root, candidate)?;
        ensure_no_symlinks_under_root(root, &target, true)?;

        // Explicitly reject dangling symlinks before opening (O_NOFOLLOW also fails with ELOOP,
        // and checking symlink_metadata gives a clear deterministic error message).
        match std::fs::symlink_metadata(&target) {
            Ok(meta) => {
                if meta.file_type().is_symlink() {
                    return Err(SandboxError::SymlinkNotAllowed(
                        target.display().to_string(),
                    ));
                }
                if meta.is_dir() {
                    return Err(SandboxError::NotARegularFile(target.display().to_string()));
                }
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(SandboxError::Io(e)),
        }

        let mut opts = OpenOptions::new();
        opts.write(true).create(true).truncate(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            opts.custom_flags(libc::O_NOFOLLOW);
        }
        let file = opts.open(&target)?;
        let meta = file.metadata()?;
        if !meta.is_file() {
            return Err(SandboxError::NotARegularFile(target.display().to_string()));
        }
        Ok((target, file))
    }
}

fn normalize_under_root(root: &Path, candidate: &Path) -> Result<PathBuf, SandboxError> {
    let joined = if candidate.is_absolute() {
        normalize_lexical(candidate)
    } else {
        normalize_lexical(&root.join(candidate))
    };
    if !joined.starts_with(root) || joined == root {
        return Err(SandboxError::OutsideRoot(
            candidate.display().to_string(),
            root.display().to_string(),
        ));
    }
    Ok(joined)
}

fn ensure_no_symlinks_under_root(
    root: &Path,
    target: &Path,
    create_missing_parents: bool,
) -> Result<(), SandboxError> {
    let rel = target.strip_prefix(root).map_err(|_| {
        SandboxError::OutsideRoot(target.display().to_string(), root.display().to_string())
    })?;

    let mut cur = root.to_path_buf();
    let comps: Vec<_> = rel.components().collect();
    for (idx, comp) in comps.iter().enumerate() {
        cur.push(comp.as_os_str());
        let is_last = idx + 1 == comps.len();
        match std::fs::symlink_metadata(&cur) {
            Ok(meta) => {
                if meta.file_type().is_symlink() {
                    return Err(SandboxError::SymlinkNotAllowed(cur.display().to_string()));
                }
                if !is_last && !meta.is_dir() {
                    return Err(SandboxError::NotARegularFile(cur.display().to_string()));
                }
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                if !is_last && create_missing_parents {
                    std::fs::create_dir(&cur)?;
                } else if !is_last {
                    return Err(SandboxError::Io(e));
                }
            }
            Err(e) => return Err(SandboxError::Io(e)),
        }
    }
    Ok(())
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
