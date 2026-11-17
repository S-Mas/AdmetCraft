use thiserror::Error;

/// Errors produced by `chem-core` parsing and molecular graph operations.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ChemError {
    /// Input string is empty or contains no molecular data.
    #[error("empty molecular input")]
    EmptyInput,

    /// Invalid SMILES syntax at a specific byte offset.
    #[error("invalid SMILES at position {pos}: {msg}")]
    InvalidSmiles {
        /// Byte offset in the input SMILES string.
        pos: usize,
        /// Human-readable reason.
        msg: String,
    },

    /// Unclosed ring closure tag in SMILES.
    #[error("unclosed ring closure '{ring}' in SMILES")]
    UnclosedRing {
        /// Ring closure number.
        ring: u16,
    },

    /// Unmatched parenthesis in SMILES branch.
    #[error("unmatched parenthesis in SMILES at position {pos}")]
    UnmatchedParen {
        /// Byte offset in the input SMILES string.
        pos: usize,
    },

    /// Unknown or unsupported element symbol.
    #[error("unknown element symbol: '{symbol}'")]
    UnknownElement {
        /// Element symbol encountered.
        symbol: String,
    },

    /// Invalid MOL/SDF V2000 record.
    #[error("invalid MOL/SDF record at line {line}: {msg}")]
    InvalidSdf {
        /// 1-based line index within the record.
        line: usize,
        /// Human-readable reason.
        msg: String,
    },

    /// Invalid SMARTS pattern.
    #[error("invalid SMARTS '{pattern}': {msg}")]
    InvalidSmarts {
        /// Original SMARTS string.
        pattern: String,
        /// Human-readable reason.
        msg: String,
    },

    /// Molecular input exceeds configured atom safety limit.
    #[error("molecule has {count} atoms, exceeding maximum allowed limit of {limit}")]
    TooManyAtoms {
        /// Atom count in the input.
        count: usize,
        /// Configured upper limit.
        limit: usize,
    },

    /// Internal graph indexing error.
    #[error("molecular graph error: {0}")]
    Graph(String),
}
