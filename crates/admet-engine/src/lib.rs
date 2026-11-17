//! Layer 3 (`admet-engine`): Session state, batch management, undo/redo history,
//! unified `CommandRegistry`, and filesystem sandboxing.

#![forbid(unsafe_code)]

pub mod sandbox;

use std::path::{Path, PathBuf};

use admet_calc::{MoleculeEvaluation, evaluate_batch};
use admet_report::{ExportFormat, ReportError, write_report_file};
use chem_core::{ChemError, Molecule, parse_csv_batch, parse_sdf, parse_smiles_batch};
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub use sandbox::{PathSandbox, SandboxError};

/// Errors produced by `admet-engine`.
#[derive(Debug, Error)]
pub enum EngineError {
    #[error("chemistry error: {0}")]
    Chem(#[from] ChemError),
    #[error("report error: {0}")]
    Report(#[from] ReportError),
    #[error("sandbox error: {0}")]
    Sandbox(#[from] SandboxError),
    #[error("nothing to undo")]
    NothingToUndo,
    #[error("nothing to redo")]
    NothingToRedo,
    #[error("invalid molecule index {0}")]
    InvalidIndex(usize),
}

/// Failed input line diagnostic.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParseFailure {
    pub line: usize,
    pub input: String,
    pub error: String,
}

/// Snapshot stored on undo/redo stacks.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct SessionSnapshot {
    smiles_input: String,
    evaluations: Vec<MoleculeEvaluation>,
    failures: Vec<ParseFailure>,
    selected_index: Option<usize>,
}

/// Active AdmetCraft session holding current molecules, evaluations, and undo/redo history.
#[derive(Debug, Clone)]
pub struct AdmetSession {
    pub smiles_input: String,
    pub evaluations: Vec<MoleculeEvaluation>,
    pub failures: Vec<ParseFailure>,
    pub selected_index: Option<usize>,
    undo_stack: Vec<SessionSnapshot>,
    redo_stack: Vec<SessionSnapshot>,
}

impl Default for AdmetSession {
    fn default() -> Self {
        Self::new()
    }
}

impl AdmetSession {
    /// Create an empty session.
    pub fn new() -> Self {
        Self {
            smiles_input: String::new(),
            evaluations: Vec::new(),
            failures: Vec::new(),
            selected_index: None,
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
        }
    }

    /// True if at least one undo step is available.
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// True if at least one redo step is available.
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    fn snapshot(&self) -> SessionSnapshot {
        SessionSnapshot {
            smiles_input: self.smiles_input.clone(),
            evaluations: self.evaluations.clone(),
            failures: self.failures.clone(),
            selected_index: self.selected_index,
        }
    }

    fn restore(&mut self, snap: SessionSnapshot) {
        self.smiles_input = snap.smiles_input;
        self.evaluations = snap.evaluations;
        self.failures = snap.failures;
        self.selected_index = snap.selected_index;
    }

    fn push_undo(&mut self) {
        self.undo_stack.push(self.snapshot());
        if self.undo_stack.len() > 64 {
            self.undo_stack.remove(0);
        }
        self.redo_stack.clear();
    }
}

/// Filter criteria for druglikeness / ADME screening.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct FilterCriteria {
    pub max_lipinski_violations: Option<usize>,
    pub require_bbb: Option<bool>,
    pub require_gi_high: Option<bool>,
    pub no_pains: bool,
    pub no_brenk: bool,
}

impl FilterCriteria {
    /// Return true if `eval` satisfies every active predicate.
    pub fn matches(&self, eval: &MoleculeEvaluation) -> bool {
        if let Some(max_v) = self.max_lipinski_violations
            && eval.druglikeness.lipinski.violations > max_v
        {
            return false;
        }
        if let Some(req_bbb) = self.require_bbb
            && eval.pharmacokinetics.bbb_permeant != req_bbb
        {
            return false;
        }
        if let Some(req_gi) = self.require_gi_high {
            let is_high = eval.pharmacokinetics.gi_absorption == admet_calc::AbsorptionLevel::High;
            if is_high != req_gi {
                return false;
            }
        }
        if self.no_pains && !eval.medicinal_chemistry.pains_alerts.is_empty() {
            return false;
        }
        if self.no_brenk && !eval.medicinal_chemistry.brenk_alerts.is_empty() {
            return false;
        }
        true
    }
}

/// Deterministic command enum shared across UI, command palette (`Ctrl+K`), CLI, and MCP.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Command {
    EvaluateBatchText { input: String },
    LoadExampleBatch,
    ImportSdf { content: String },
    ImportCsv { content: String },
    FilterBatch { criteria: FilterCriteria },
    SelectMolecule { index: usize },
    ClearSession,
    Undo,
    Redo,
}

/// Metadata describing an executable command for command palette discovery.
#[derive(Debug, Clone)]
pub struct CommandDescriptor {
    pub id: &'static str,
    pub title: &'static str,
    pub shortcut: &'static str,
    pub category: &'static str,
}

/// Outcome returned after executing a command.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommandOutcome {
    pub message: String,
    pub evaluated_count: usize,
    pub failure_count: usize,
}

/// Central deterministic command dispatcher.
pub struct CommandRegistry;

impl CommandRegistry {
    /// Return static metadata for palette-accessible commands.
    pub fn catalog() -> &'static [CommandDescriptor] {
        &[
            CommandDescriptor {
                id: "load_examples",
                title: "Load Example Drugs (15 Reference Molecules)",
                shortcut: "Ctrl+E",
                category: "Input",
            },
            CommandDescriptor {
                id: "run_batch",
                title: "Run ADMET Calculation",
                shortcut: "Ctrl+Enter",
                category: "Compute",
            },
            CommandDescriptor {
                id: "filter_lipinski_0",
                title: "Filter: Lipinski 0 Violations & No PAINS",
                shortcut: "",
                category: "Filter",
            },
            CommandDescriptor {
                id: "filter_bbb_yes",
                title: "Filter: BBB Permeant (Yolk Only)",
                shortcut: "",
                category: "Filter",
            },
            CommandDescriptor {
                id: "clear_session",
                title: "Clear All Molecules",
                shortcut: "",
                category: "Session",
            },
            CommandDescriptor {
                id: "undo",
                title: "Undo Last Action",
                shortcut: "Ctrl+Z",
                category: "History",
            },
            CommandDescriptor {
                id: "redo",
                title: "Redo Action",
                shortcut: "Ctrl+Y",
                category: "History",
            },
        ]
    }

    /// Execute `cmd` against `session`.
    pub fn execute(session: &mut AdmetSession, cmd: Command) -> Result<CommandOutcome, EngineError> {
        match cmd {
            Command::EvaluateBatchText { input } => {
                session.push_undo();
                session.smiles_input = input.clone();
                let parsed = parse_smiles_batch(&input);
                let mut mols: Vec<Molecule> = Vec::new();
                let mut failures: Vec<ParseFailure> = Vec::new();
                for (line, raw, res) in parsed {
                    match res {
                        Ok(m) => mols.push(m),
                        Err(e) => failures.push(ParseFailure {
                            line,
                            input: raw,
                            error: e.to_string(),
                        }),
                    }
                }
                let evals = evaluate_batch(&mols);
                session.evaluations = evals;
                session.failures = failures;
                session.selected_index = if session.evaluations.is_empty() {
                    None
                } else {
                    Some(0)
                };
                Ok(CommandOutcome {
                    message: format!(
                        "Evaluated {} molecule(s) ({} error(s))",
                        session.evaluations.len(),
                        session.failures.len()
                    ),
                    evaluated_count: session.evaluations.len(),
                    failure_count: session.failures.len(),
                })
            }
            Command::LoadExampleBatch => {
                let sample = EXAMPLE_BATCH.to_string();
                Self::execute(session, Command::EvaluateBatchText { input: sample })
            }
            Command::ImportSdf { content } => {
                session.push_undo();
                let parsed = parse_sdf(&content);
                let mut mols = Vec::new();
                let mut failures = Vec::new();
                for (idx, res) in parsed.into_iter().enumerate() {
                    match res {
                        Ok(m) => mols.push(m),
                        Err(e) => failures.push(ParseFailure {
                            line: idx + 1,
                            input: format!("SDF record {}", idx + 1),
                            error: e.to_string(),
                        }),
                    }
                }
                let evals = evaluate_batch(&mols);
                session.smiles_input = evals
                    .iter()
                    .map(|e| format!("{} {}", e.canonical_smiles, e.name))
                    .collect::<Vec<_>>()
                    .join("\n");
                session.evaluations = evals;
                session.failures = failures;
                session.selected_index = if session.evaluations.is_empty() {
                    None
                } else {
                    Some(0)
                };
                Ok(CommandOutcome {
                    message: format!("Imported {} SDF record(s)", session.evaluations.len()),
                    evaluated_count: session.evaluations.len(),
                    failure_count: session.failures.len(),
                })
            }
            Command::ImportCsv { content } => {
                session.push_undo();
                let parsed = parse_csv_batch(&content)?;
                let mut mols = Vec::new();
                let mut failures = Vec::new();
                for (idx, res) in parsed.into_iter().enumerate() {
                    match res {
                        Ok(m) => mols.push(m),
                        Err(e) => failures.push(ParseFailure {
                            line: idx + 2,
                            input: format!("CSV row {}", idx + 2),
                            error: e.to_string(),
                        }),
                    }
                }
                let evals = evaluate_batch(&mols);
                session.smiles_input = evals
                    .iter()
                    .map(|e| format!("{} {}", e.canonical_smiles, e.name))
                    .collect::<Vec<_>>()
                    .join("\n");
                session.evaluations = evals;
                session.failures = failures;
                session.selected_index = if session.evaluations.is_empty() {
                    None
                } else {
                    Some(0)
                };
                Ok(CommandOutcome {
                    message: format!("Imported {} CSV row(s)", session.evaluations.len()),
                    evaluated_count: session.evaluations.len(),
                    failure_count: session.failures.len(),
                })
            }
            Command::FilterBatch { criteria } => {
                session.push_undo();
                session.evaluations.retain(|e| criteria.matches(e));
                session.smiles_input = session
                    .evaluations
                    .iter()
                    .map(|e| format!("{} {}", e.canonical_smiles, e.name))
                    .collect::<Vec<_>>()
                    .join("\n");
                session.selected_index = if session.evaluations.is_empty() {
                    None
                } else {
                    Some(0)
                };
                Ok(CommandOutcome {
                    message: format!("Filtered to {} molecule(s)", session.evaluations.len()),
                    evaluated_count: session.evaluations.len(),
                    failure_count: 0,
                })
            }
            Command::SelectMolecule { index } => {
                if index >= session.evaluations.len() {
                    return Err(EngineError::InvalidIndex(index));
                }
                session.selected_index = Some(index);
                Ok(CommandOutcome {
                    message: format!("Selected '{}'", session.evaluations[index].name),
                    evaluated_count: session.evaluations.len(),
                    failure_count: session.failures.len(),
                })
            }
            Command::ClearSession => {
                session.push_undo();
                session.smiles_input.clear();
                session.evaluations.clear();
                session.failures.clear();
                session.selected_index = None;
                Ok(CommandOutcome {
                    message: "Cleared session".to_string(),
                    evaluated_count: 0,
                    failure_count: 0,
                })
            }
            Command::Undo => {
                let prev = session
                    .undo_stack
                    .pop()
                    .ok_or(EngineError::NothingToUndo)?;
                let current = session.snapshot();
                session.redo_stack.push(current);
                session.restore(prev);
                Ok(CommandOutcome {
                    message: "Undid last action".to_string(),
                    evaluated_count: session.evaluations.len(),
                    failure_count: session.failures.len(),
                })
            }
            Command::Redo => {
                let next = session
                    .redo_stack
                    .pop()
                    .ok_or(EngineError::NothingToRedo)?;
                let current = session.snapshot();
                session.undo_stack.push(current);
                session.restore(next);
                Ok(CommandOutcome {
                    message: "Redid action".to_string(),
                    evaluated_count: session.evaluations.len(),
                    failure_count: session.failures.len(),
                })
            }
        }
    }

    /// Export session evaluations safely using optional `PathSandbox`.
    pub fn export_session(
        session: &AdmetSession,
        path: &Path,
        format: ExportFormat,
        sandbox: Option<&PathSandbox>,
    ) -> Result<PathBuf, EngineError> {
        let resolved = match sandbox {
            Some(sb) => sb.resolve_write(path)?,
            None => path.to_path_buf(),
        };
        write_report_file(&session.evaluations, &resolved, format)?;
        Ok(resolved)
    }
}

const EXAMPLE_BATCH: &str = "\
CC(=O)Oc1ccccc1C(=O)O Aspirin
Cn1cnc2c1c(=O)n(C)c(=O)n2C Caffeine
CC(C)Cc1ccc(C(C)C(=O)O)cc1 Ibuprofen
CC(=O)Nc1ccc(O)cc1 Paracetamol
CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O Penicillin_G
CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O Atorvastatin
CN(C)C(=N)NC(=N)N Metformin
CN1c2ccc(Cl)cc2C(=NCC1=O)c1ccccc1 Diazepam
CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O Warfarin
COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1 Omeprazole
CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12 Sildenafil
CC/C(=C(/c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1 Tamoxifen
CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5 Morphine
O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O Ciprofloxacin
CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1 Fluoxetine";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn undo_redo_and_filter_work() -> Result<(), Box<dyn std::error::Error>> {
        let mut session = AdmetSession::new();
        let out = CommandRegistry::execute(&mut session, Command::LoadExampleBatch)?;
        assert_eq!(out.evaluated_count, 15);
        assert!(session.can_undo());

        CommandRegistry::execute(
            &mut session,
            Command::FilterBatch {
                criteria: FilterCriteria {
                    max_lipinski_violations: Some(0),
                    require_bbb: Some(true),
                    require_gi_high: Some(true),
                    no_pains: true,
                    no_brenk: false,
                },
            },
        )?;
        assert!(session.evaluations.len() < 15);
        assert!(!session.evaluations.is_empty());

        CommandRegistry::execute(&mut session, Command::Undo)?;
        assert_eq!(session.evaluations.len(), 15);

        CommandRegistry::execute(&mut session, Command::Redo)?;
        assert!(session.evaluations.len() < 15);
        Ok(())
    }

    #[test]
    fn sandbox_rejects_escape_paths() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("admetcraft_sb_test");
        let read_dir = base.join("read");
        let write_dir = base.join("write");
        std::fs::create_dir_all(&read_dir)?;
        std::fs::create_dir_all(&write_dir)?;
        let inside_file = read_dir.join("ok.smi");
        std::fs::write(&inside_file, "CCO Ethanol\n")?;

        let sb = PathSandbox::new(Some(&read_dir), Some(&write_dir))?;
        assert!(sb.resolve_read(Path::new("ok.smi")).is_ok());
        assert!(sb.resolve_read(Path::new("../write")).is_err());
        assert!(sb.resolve_write(Path::new("out.csv")).is_ok());
        assert!(sb.resolve_write(Path::new("../read/evil.csv")).is_err());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
