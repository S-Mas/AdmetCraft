//! Layer 3 (`admet-engine`): Session state, batch management, undo/redo history,
//! unified `CommandRegistry`, and filesystem sandboxing.

#![forbid(unsafe_code)]

pub mod sandbox;

use std::path::{Path, PathBuf};

use admet_calc::{MoleculeEvaluation, TppPreset};
use admet_qsar::{compute_tpp_mpo, evaluate_batch_auto};
use admet_report::{ExportFormat, ReportError, write_report_file, write_report_to_writer};
use chem_core::{ChemError, Molecule, parse_csv_batch, parse_sdf, parse_smiles_batch};
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub use sandbox::{MAX_READ_BYTES, PathSandbox, SandboxError};

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
    pub tpp_preset: TppPreset,
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
            tpp_preset: TppPreset::OralSmallMolecule,
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
        }
    }

    /// Update active Target Product Profile preset and recompute TPP MPO on all loaded molecules.
    pub fn set_tpp_preset(&mut self, preset: TppPreset) {
        self.tpp_preset = preset;
        for ev in &mut self.evaluations {
            ev.tpp_mpo = Some(compute_tpp_mpo(ev, preset));
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
    #[serde(default)]
    pub no_pains: bool,
    #[serde(default)]
    pub no_brenk: bool,
    #[serde(default)]
    pub min_qed: Option<f64>,
    #[serde(default)]
    pub min_cns_mpo: Option<f64>,
    #[serde(default)]
    pub no_high_herg: bool,
    #[serde(default)]
    pub no_ames: bool,
    #[serde(default)]
    pub no_dili: bool,
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
        if let Some(min_qed) = self.min_qed
            && eval.druglikeness.qed < min_qed
        {
            return false;
        }
        if let Some(min_mpo) = self.min_cns_mpo
            && eval.druglikeness.cns_mpo < min_mpo
        {
            return false;
        }
        if self.no_high_herg && eval.safety.herg_risk == admet_calc::HergRisk::High {
            return false;
        }
        if self.no_ames && !eval.safety.ames_alerts.is_empty() {
            return false;
        }
        if self.no_dili && !eval.safety.dili_alerts.is_empty() {
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
                id: "filter_qed_05",
                title: "Filter: Drug-like QED >= 0.50 & No PAINS",
                shortcut: "",
                category: "Filter",
            },
            CommandDescriptor {
                id: "filter_cns_mpo_4",
                title: "Filter: CNS MPO >= 4.0 & BBB+",
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
    pub fn execute(
        session: &mut AdmetSession,
        cmd: Command,
    ) -> Result<CommandOutcome, EngineError> {
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
                let evals = evaluate_batch_auto(&mols, session.tpp_preset);
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
                let evals = evaluate_batch_auto(&mols, session.tpp_preset);
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
                let evals = evaluate_batch_auto(&mols, session.tpp_preset);
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
                let prev = session.undo_stack.pop().ok_or(EngineError::NothingToUndo)?;
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
                let next = session.redo_stack.pop().ok_or(EngineError::NothingToRedo)?;
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

    /// Export session evaluations safely using optional `PathSandbox` (atomically via `open_write`
    /// when sandboxed to eliminate path TOCTOU races).
    pub fn export_session(
        session: &AdmetSession,
        path: &Path,
        format: ExportFormat,
        sandbox: Option<&PathSandbox>,
    ) -> Result<PathBuf, EngineError> {
        match sandbox {
            Some(sb) => {
                let (resolved, mut file) = sb.open_write(path)?;
                write_report_to_writer(&session.evaluations, &mut file, format)?;
                Ok(resolved)
            }
            None => {
                write_report_file(&session.evaluations, path, format)?;
                Ok(path.to_path_buf())
            }
        }
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
                    ..FilterCriteria::default()
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
    fn sandbox_rejects_escape_paths_dangling_symlinks_and_fifos()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("admetcraft_sb_test_p0");
        let _ = std::fs::remove_dir_all(&base);
        let read_dir = base.join("read");
        let write_dir = base.join("write");
        let secret_dir = base.join("secret");
        std::fs::create_dir_all(&read_dir)?;
        std::fs::create_dir_all(&write_dir)?;
        std::fs::create_dir_all(&secret_dir)?;

        let inside_file = read_dir.join("ok.smi");
        std::fs::write(&inside_file, "CCO Ethanol\n")?;

        let sb = PathSandbox::new(Some(&read_dir), Some(&write_dir))?;
        assert!(sb.open_read(Path::new("ok.smi")).is_ok());
        assert!(sb.open_read(Path::new("../write")).is_err());

        let (wp, mut wf) = sb.open_write(Path::new("out.csv"))?;
        use std::io::Write;
        wf.write_all(b"ok")?;
        assert!(wp.starts_with(std::fs::canonicalize(&write_dir)?));
        assert!(sb.open_write(Path::new("../read/evil.csv")).is_err());

        #[cfg(unix)]
        {
            use std::os::unix::fs::symlink;
            // 1. Dangling symlink inside write_root pointing outside -> MUST be rejected!
            let dangling = write_dir.join("dangling_out.csv");
            symlink(secret_dir.join("pwned.csv"), &dangling)?;
            assert!(matches!(
                sb.open_write(Path::new("dangling_out.csv")),
                Err(SandboxError::SymlinkNotAllowed(_))
            ));
            assert!(!secret_dir.join("pwned.csv").exists());

            // 2. Symlinked directory inside write_root pointing outside -> MUST be rejected!
            let sym_dir = write_dir.join("linked_dir");
            symlink(&secret_dir, &sym_dir)?;
            let linked_res = sb.open_write(Path::new("linked_dir/evil.csv"));
            assert!(
                matches!(linked_res, Err(SandboxError::SymlinkNotAllowed(_))),
                "expected SymlinkNotAllowed, got: {linked_res:?}"
            );
            assert!(!secret_dir.join("evil.csv").exists());

            // 3. Symlink inside read_root pointing outside -> MUST be rejected!
            let read_sym = read_dir.join("leak.smi");
            std::fs::write(secret_dir.join("secret.smi"), "CCO Secret\n")?;
            symlink(secret_dir.join("secret.smi"), &read_sym)?;
            assert!(sb.open_read(Path::new("leak.smi")).is_err());

            // 4. FIFO inside read_root -> MUST be rejected immediately without blocking!
            let fifo_path = read_dir.join("pipe.smi");
            let c_path = std::ffi::CString::new(fifo_path.to_string_lossy().as_bytes())?;
            let mk_rc = unsafe_free_mkfifo(&c_path)?;
            assert_eq!(mk_rc, 0);
            assert!(matches!(
                sb.open_read(Path::new("pipe.smi")),
                Err(SandboxError::NotARegularFile(_))
            ));
        }

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[cfg(unix)]
    fn unsafe_free_mkfifo(path: &std::ffi::CStr) -> Result<i32, Box<dyn std::error::Error>> {
        // Call system `mkfifo` binary so we preserve #![forbid(unsafe_code)]!
        let status = std::process::Command::new("mkfifo")
            .arg(path.to_str()?)
            .status()?;
        Ok(status.code().unwrap_or(1))
    }

    #[cfg(unix)]
    #[test]
    fn sandbox_openat2_blocks_4000_intermediate_directory_symlink_swap_races()
    -> Result<(), Box<dyn std::error::Error>> {
        use std::os::unix::fs::symlink;
        use std::sync::{
            Arc,
            atomic::{AtomicBool, Ordering},
        };

        let base = std::env::temp_dir().join("admetcraft_sb_race_4k");
        let _ = std::fs::remove_dir_all(&base);
        let read_dir = base.join("read");
        let write_dir = base.join("write");
        let outside_dir = base.join("outside");
        let real_r_sub = read_dir.join("real_sub");
        let real_w_sub = write_dir.join("real_sub");
        std::fs::create_dir_all(&real_r_sub)?;
        std::fs::create_dir_all(&real_w_sub)?;
        std::fs::create_dir_all(&outside_dir)?;
        std::fs::write(real_r_sub.join("target.smi"), "INSIDE")?;
        std::fs::write(outside_dir.join("target.smi"), "SECRET_OUTSIDE")?;

        let r_link = read_dir.join("sub");
        let w_link = write_dir.join("sub");
        std::fs::create_dir_all(&r_link)?;
        std::fs::create_dir_all(&w_link)?;

        let sb = PathSandbox::new(Some(&read_dir), Some(&write_dir))?;
        let stop = Arc::new(AtomicBool::new(false));
        let stop_r = Arc::clone(&stop);
        let stop_w = Arc::clone(&stop);

        let racer = std::thread::spawn(move || {
            let mut toggle = false;
            while !stop_r.load(Ordering::Relaxed) {
                toggle = !toggle;
                let _ = std::fs::remove_file(&r_link);
                let _ = std::fs::remove_dir_all(&r_link);
                let _ = std::fs::remove_file(&w_link);
                let _ = std::fs::remove_dir_all(&w_link);
                if toggle {
                    let _ = symlink(&outside_dir, &r_link);
                    let _ = symlink(&outside_dir, &w_link);
                } else {
                    let _ = std::fs::create_dir(&r_link);
                    let _ = std::fs::write(r_link.join("target.smi"), "INSIDE");
                    let _ = std::fs::create_dir(&w_link);
                }
            }
        });

        for _ in 0..4_000 {
            if let Ok((_, content)) = sb.read_to_string_limited(Path::new("sub/target.smi"), 1024) {
                assert_ne!(
                    content, "SECRET_OUTSIDE",
                    "intermediate-dir symlink race leaked outside file!"
                );
            }
            if let Ok((_, mut f)) = sb.open_write(Path::new("sub/pwned.csv")) {
                use std::io::Write;
                let _ = f.write_all(b"x");
            }
        }
        stop_w.store(true, Ordering::Relaxed);
        let _ = racer.join();

        assert!(
            !base.join("outside/pwned.csv").exists(),
            "intermediate-dir symlink race wrote outside sandbox!"
        );
        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
