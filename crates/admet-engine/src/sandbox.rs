use std::path::{Component, Path, PathBuf};

use thiserror::Error;

/// Errors raised when validating sandboxed filesystem paths.
#[derive(Debug, Error)]
pub enum SandboxError {
    #[error("path '{0}' is outside authorized root '{1}'")]
    OutsideRoot(String, String),
    #[error("no root directory configured for this operation")]
    NoRootConfigured,
    #[error("I/O error resolving path: {0}")]
    Io(#[from] std::io::Error),
}

/// Directory sandbox enforcing strict `--read-root` and `--write-root` isolation.
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

    /// Resolve a candidate path for reading inside `read_root`.
    pub fn resolve_read(&self, candidate: &Path) -> Result<PathBuf, SandboxError> {
        let root = self
            .read_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let target = if candidate.is_absolute() {
            candidate.to_path_buf()
        } else {
            root.join(candidate)
        };
        let canon = std::fs::canonicalize(&target)?;
        if !canon.starts_with(root) {
            return Err(SandboxError::OutsideRoot(
                candidate.display().to_string(),
                root.display().to_string(),
            ));
        }
        Ok(canon)
    }

    /// Resolve a candidate path for writing inside `write_root`.
    /// Works even when the target file does not exist yet by canonicalizing its existing ancestor.
    pub fn resolve_write(&self, candidate: &Path) -> Result<PathBuf, SandboxError> {
        let root = self
            .write_root
            .as_ref()
            .ok_or(SandboxError::NoRootConfigured)?;
        let joined = if candidate.is_absolute() {
            normalize_lexical(candidate)
        } else {
            normalize_lexical(&root.join(candidate))
        };

        // Walk up until we find an existing ancestor to canonicalize (prevents symlink escape)
        let mut suffix_parts = Vec::new();
        let mut probe = joined.as_path();
        while !probe.exists() {
            let Some(name) = probe.file_name() else {
                break;
            };
            suffix_parts.push(name.to_os_string());
            let Some(parent) = probe.parent() else {
                break;
            };
            probe = parent;
        }
        let mut canon = std::fs::canonicalize(probe)?;
        for part in suffix_parts.into_iter().rev() {
            canon.push(part);
        }

        if !canon.starts_with(root) {
            return Err(SandboxError::OutsideRoot(
                candidate.display().to_string(),
                root.display().to_string(),
            ));
        }
        Ok(canon)
    }
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
