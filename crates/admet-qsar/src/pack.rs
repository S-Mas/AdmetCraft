#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::OnceLock;
use std::time::Instant;

use admet_calc::{
    ApplicabilityDomain, ApplicabilityDomainStatus, QsarAdmetReport, QsarTaskKind, Tox21Summary,
    evaluate_molecule,
};
use chem_core::parse_smiles;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::data_dir::default_model_pack_path;
use crate::fingerprint::{MolFeatures, extract_mol_features};
use crate::gbdt::{GbdtTrainConfig, TrainedEndpointModel, train_endpoint_gbdt};

/// Errors emitted during QSAR training or model-pack loading.
#[derive(Debug, Error)]
pub enum QsarError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("dataset error for '{endpoint}': {message}")]
    Dataset { endpoint: String, message: String },
}

/// Specification of one official TDC ADMET or Tox21 endpoint.
#[derive(Debug, Clone, Copy)]
pub struct TdcEndpointSpec {
    pub id: &'static str,
    pub file_name: &'static str,
    pub display_name: &'static str,
    pub category: &'static str,
    pub task: QsarTaskKind,
    pub units: &'static str,
    pub log10_transform: bool,
    pub species_filter: Option<&'static str>,
    /// Optional explicit label column name for multi-task TDC tables (e.g. `tox21.tab`);
    /// when `None`, defaults to standard TDC `"Y"` column.
    pub target_col: Option<&'static str>,
}

/// All 32 official Therapeutics Data Commons (TDC) endpoints: 20 core ADMET + 12 Tox21 pathways.
pub const TDC_ENDPOINTS: &[TdcEndpointSpec] = &[
    TdcEndpointSpec {
        id: "herg_karim",
        file_name: "herg_karim.tab",
        display_name: "hERG Cardiotoxicity (Karim)",
        category: "Toxicity",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "ames",
        file_name: "ames.tab",
        display_name: "AMES Mutagenicity",
        category: "Toxicity",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "dili",
        file_name: "dili.tab",
        display_name: "DILI (Liver Injury)",
        category: "Toxicity",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "clintox",
        file_name: "clintox.tab",
        display_name: "ClinTox (Clinical Toxicity)",
        category: "Toxicity",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "cyp1a2_veith",
        file_name: "cyp1a2_veith.tab",
        display_name: "CYP1A2 Inhibitor",
        category: "Metabolism",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "cyp2c9_veith",
        file_name: "cyp2c9_veith.tab",
        display_name: "CYP2C9 Inhibitor",
        category: "Metabolism",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "cyp2c19_veith",
        file_name: "cyp2c19_veith.tab",
        display_name: "CYP2C19 Inhibitor",
        category: "Metabolism",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "cyp2d6_veith",
        file_name: "cyp2d6_veith.tab",
        display_name: "CYP2D6 Inhibitor",
        category: "Metabolism",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "cyp3a4_veith",
        file_name: "cyp3a4_veith.tab",
        display_name: "CYP3A4 Inhibitor",
        category: "Metabolism",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "pgp_broccatelli",
        file_name: "pgp_broccatelli.tab",
        display_name: "P-gp Inhibitor / Efflux",
        category: "Absorption",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "hia_hou",
        file_name: "hia_hou.tab",
        display_name: "HIA (Human Intestinal Abs.)",
        category: "Absorption",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "bioavailability_ma",
        file_name: "bioavailability_ma.tab",
        display_name: "Oral Bioavailability (F20%)",
        category: "Absorption",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "caco2_wang",
        file_name: "caco2_wang.tab",
        display_name: "Caco-2 Permeability",
        category: "Absorption",
        task: QsarTaskKind::Regression,
        units: "log(cm/s)",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "bbb_martins",
        file_name: "bbb_martins.tab",
        display_name: "BBB Penetration",
        category: "Distribution",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "ppbr_az",
        file_name: "ppbr_az.tab",
        display_name: "PPB Human (% Bound)",
        category: "Distribution",
        task: QsarTaskKind::Regression,
        units: "%",
        log10_transform: false,
        species_filter: Some("Homo sapiens"),
        target_col: None,
    },
    TdcEndpointSpec {
        id: "vdss_lombardo",
        file_name: "vdss_lombardo.tab",
        display_name: "VDss (Volume of Distr.)",
        category: "Distribution",
        task: QsarTaskKind::Regression,
        units: "L/kg",
        log10_transform: true,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "clearance_microsome_az",
        file_name: "clearance_microsome_az.tab",
        display_name: "Microsomal CLint (AZ)",
        category: "Excretion",
        task: QsarTaskKind::Regression,
        units: "uL/min/mg",
        log10_transform: true,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "half_life_obach",
        file_name: "half_life_obach.tab",
        display_name: "Human Half-Life (t1/2)",
        category: "Excretion",
        task: QsarTaskKind::Regression,
        units: "h",
        log10_transform: true,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "solubility_aqsoldb",
        file_name: "solubility_aqsoldb.tab",
        display_name: "Aqueous Solubility (AqSolDB)",
        category: "Physicochem",
        task: QsarTaskKind::Regression,
        units: "log(mol/L)",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    TdcEndpointSpec {
        id: "lipophilicity_astrazeneca",
        file_name: "lipophilicity_astrazeneca.tab",
        display_name: "Lipophilicity LogD7.4 (AZ)",
        category: "Physicochem",
        task: QsarTaskKind::Regression,
        units: "LogD7.4",
        log10_transform: false,
        species_filter: None,
        target_col: None,
    },
    // --- Tox21 Nuclear Receptors (Endocrine / Metabolic Disruption, 7 assays) ---
    TdcEndpointSpec {
        id: "tox21_nr_ar",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-AR (Androgen Receptor)",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-AR"),
    },
    TdcEndpointSpec {
        id: "tox21_nr_ar_lbd",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-AR-LBD (Androgen LBD)",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-AR-LBD"),
    },
    TdcEndpointSpec {
        id: "tox21_nr_ahr",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-AhR (Aryl Hydrocarbon R.)",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-AhR"),
    },
    TdcEndpointSpec {
        id: "tox21_nr_aromatase",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-Aromatase Inhibition",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-Aromatase"),
    },
    TdcEndpointSpec {
        id: "tox21_nr_er",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-ER (Estrogen Receptor)",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-ER"),
    },
    TdcEndpointSpec {
        id: "tox21_nr_er_lbd",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-ER-LBD (Estrogen LBD)",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-ER-LBD"),
    },
    TdcEndpointSpec {
        id: "tox21_nr_ppar_gamma",
        file_name: "tox21.tab",
        display_name: "Tox21 NR-PPAR-gamma",
        category: "Tox21-NR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("NR-PPAR-gamma"),
    },
    // --- Tox21 Cellular Stress Response Pathways (5 assays) ---
    TdcEndpointSpec {
        id: "tox21_sr_are",
        file_name: "tox21.tab",
        display_name: "Tox21 SR-ARE (Oxidative Stress / Nrf2)",
        category: "Tox21-SR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("SR-ARE"),
    },
    TdcEndpointSpec {
        id: "tox21_sr_atad5",
        file_name: "tox21.tab",
        display_name: "Tox21 SR-ATAD5 (Genotoxicity / DNA)",
        category: "Tox21-SR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("SR-ATAD5"),
    },
    TdcEndpointSpec {
        id: "tox21_sr_hse",
        file_name: "tox21.tab",
        display_name: "Tox21 SR-HSE (Heat Shock Response)",
        category: "Tox21-SR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("SR-HSE"),
    },
    TdcEndpointSpec {
        id: "tox21_sr_mmp",
        file_name: "tox21.tab",
        display_name: "Tox21 SR-MMP (Mitochondrial Toxicity)",
        category: "Tox21-SR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("SR-MMP"),
    },
    TdcEndpointSpec {
        id: "tox21_sr_p53",
        file_name: "tox21.tab",
        display_name: "Tox21 SR-p53 (p53 DNA Damage)",
        category: "Tox21-SR",
        task: QsarTaskKind::BinaryClassification,
        units: "prob",
        log10_transform: false,
        species_filter: None,
        target_col: Some("SR-p53"),
    },
];

/// Versioned container of all trained TDC QSAR models loaded locally by AdmetCraft v1.0.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdmetModelPack {
    pub format: String,
    pub version: String,
    pub trained_at: String,
    pub split_method: String,
    pub models: Vec<TrainedEndpointModel>,
}

impl AdmetModelPack {
    pub fn load_from_path(path: &Path) -> Result<Self, QsarError> {
        let raw = std::fs::read_to_string(path)?;
        let pack: Self = serde_json::from_str(&raw)?;
        Ok(pack)
    }

    pub fn save_to_path(&self, path: &Path) -> Result<(), QsarError> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let raw = serde_json::to_string(self)?;
        std::fs::write(path, raw)?;
        Ok(())
    }

    /// Evaluates all trained QSAR endpoints on `feat` and computes a global Applicability Domain summary.
    #[must_use]
    pub fn predict_all(&self, feat: &MolFeatures, source_label: &str) -> QsarAdmetReport {
        let mut endpoints = Vec::with_capacity(self.models.len());
        let mut best_max_sim = 0.0f64;
        let mut sum_top5 = 0.0f64;
        let mut sum_z = 0.0f64;
        let mut in_dom = 0usize;
        let mut border = 0usize;

        for model in &self.models {
            let pred = model.predict(feat);
            if pred.applicability_domain.max_tanimoto > best_max_sim {
                best_max_sim = pred.applicability_domain.max_tanimoto;
            }
            sum_top5 += pred.applicability_domain.top5_tanimoto;
            sum_z += pred.applicability_domain.descriptor_z_rms;
            match pred.applicability_domain.status {
                ApplicabilityDomainStatus::InDomain => in_dom += 1,
                ApplicabilityDomainStatus::Borderline => border += 1,
                ApplicabilityDomainStatus::OutOfDomain => {}
            }
            endpoints.push(pred);
        }

        let m = self.models.len().max(1) as f64;
        let mean_top5 = ((sum_top5 / m) * 1000.0).round() / 1000.0;
        let mean_z = ((sum_z / m) * 100.0).round() / 100.0;
        let status = if in_dom * 2 >= self.models.len() {
            ApplicabilityDomainStatus::InDomain
        } else if (in_dom + border) * 2 >= self.models.len() {
            ApplicabilityDomainStatus::Borderline
        } else {
            ApplicabilityDomainStatus::OutOfDomain
        };

        QsarAdmetReport {
            model_pack_version: self.version.clone(),
            model_pack_source: source_label.to_string(),
            overall_ad: ApplicabilityDomain {
                status,
                max_tanimoto: best_max_sim,
                top5_tanimoto: mean_top5,
                descriptor_z_rms: mean_z,
            },
            endpoints,
        }
    }

    /// Renders a clean Markdown scientific validation summary table of all endpoints.
    #[must_use]
    pub fn render_validation_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str("# AdmetCraft v1.1 — TDC QSAR & Tox21 Validation Report (Bemis-Murcko Scaffold Split 70/10/20)\n\n");
        out.push_str(&format!(
            "- **Model Pack Version:** `{}` | **Split:** `{}` | **Endpoints:** `{}`\n\n",
            self.version,
            self.split_method,
            self.models.len()
        ));
        out.push_str("| Endpoint | Category | Task | Train / Cal / Test | Trees | Held-Out Scaffold Metric | Secondary Metric | 90% Conformal Cov. |\n");
        out.push_str("|---|---|---|---:|---:|---:|---:|---:|\n");

        let mut cls_aucs = Vec::new();
        let mut reg_rhos = Vec::new();
        let mut covs = Vec::new();

        for m in &self.models {
            let mm = &m.metrics;
            covs.push(mm.conformal_90_coverage);
            let (primary, secondary, task_s) = match m.task {
                QsarTaskKind::BinaryClassification => {
                    let auc = mm.roc_auc.unwrap_or(0.5);
                    cls_aucs.push(auc);
                    (
                        format!("ROC AUC **{:.3}**", auc),
                        format!(
                            "PR {:.3} / Brier {:.3}",
                            mm.pr_auc.unwrap_or(0.0),
                            mm.brier_score.unwrap_or(0.25)
                        ),
                        "Binary",
                    )
                }
                QsarTaskKind::Regression => {
                    let rho = mm.spearman_rho.unwrap_or(0.0);
                    reg_rhos.push(rho);
                    (
                        format!("Spearman **{:.3}**", rho),
                        format!(
                            "RMSE {:.3} / MAE {:.3}",
                            mm.rmse.unwrap_or(0.0),
                            mm.mae.unwrap_or(0.0)
                        ),
                        if m.log10_transform {
                            "Reg (log10)"
                        } else {
                            "Regression"
                        },
                    )
                }
            };
            out.push_str(&format!(
                "| `{}` ({}) | {} | {} | {} / {} / {} | {} | {} | {} | **{:.1}%** |\n",
                m.endpoint_id,
                m.display_name,
                m.category,
                task_s,
                mm.train_n,
                mm.cal_n,
                mm.test_n,
                mm.num_trees_used,
                primary,
                secondary,
                mm.conformal_90_coverage * 100.0
            ));
        }

        let mean_auc = if cls_aucs.is_empty() {
            0.0
        } else {
            cls_aucs.iter().sum::<f64>() / (cls_aucs.len() as f64)
        };
        let mean_rho = if reg_rhos.is_empty() {
            0.0
        } else {
            reg_rhos.iter().sum::<f64>() / (reg_rhos.len() as f64)
        };
        let mean_cov = if covs.is_empty() {
            0.0
        } else {
            covs.iter().sum::<f64>() / (covs.len() as f64)
        };
        out.push_str(&format!(
            "\n**Summary:** Mean Classification Scaffold ROC AUC = **{mean_auc:.3}** ({} endpoints) | Mean Regression Scaffold Spearman rho = **{mean_rho:.3}** ({} endpoints) | Mean 90% Conformal Coverage = **{:.1}%**.\n",
            cls_aucs.len(),
            reg_rhos.len(),
            mean_cov * 100.0
        ));
        out
    }
}

static DEFAULT_PACK_CACHE: OnceLock<Option<AdmetModelPack>> = OnceLock::new();

/// Loads and caches the default local TDC model pack from `CRAFT_DATA_DIR/models/admetcraft_v1/admetcraft_tdc_v1.json`.
/// Returns `None` gracefully if no model pack has been trained/installed yet.
#[must_use]
pub fn load_default_model_pack() -> Option<&'static AdmetModelPack> {
    DEFAULT_PACK_CACHE
        .get_or_init(|| {
            let p = default_model_pack_path();
            if p.is_file() {
                AdmetModelPack::load_from_path(&p).ok()
            } else {
                None
            }
        })
        .as_ref()
}

fn unquote(s: &str) -> &str {
    s.trim().trim_matches('"').trim()
}

/// Summarizes the 12 Tox21 Nuclear Receptor (`Tox21-NR`) and Cellular Stress (`Tox21-SR`)
/// assays from a QSAR report, returning `None` if no Tox21 endpoints are loaded.
#[must_use]
pub fn summarize_tox21(report: &QsarAdmetReport) -> Option<Tox21Summary> {
    let mut evaluated = 0usize;
    let mut endocrine_alerts = Vec::new();
    let mut stress_alerts = Vec::new();
    let mut max_prob = 0.0_f64;

    for ep in &report.endpoints {
        if !ep.endpoint_id.starts_with("tox21_") {
            continue;
        }
        evaluated += 1;
        if ep.value > max_prob {
            max_prob = ep.value;
        }
        // Flag assay when predicted probability >= 0.50 and conformal set is not confidently Negative
        if ep.value >= 0.50 && ep.conformal_call != "Negative" {
            let tag = format!("{} (P={:.2})", ep.display_name, ep.value);
            if ep.category == "Tox21-NR" {
                endocrine_alerts.push(tag);
            } else {
                stress_alerts.push(tag);
            }
        }
    }

    if evaluated == 0 {
        None
    } else {
        Some(Tox21Summary {
            assays_evaluated: evaluated,
            endocrine_alerts,
            stress_alerts,
            max_probability: (max_prob * 1000.0).round() / 1000.0,
        })
    }
}

/// Parses a TDC `.tab` file into `(smiles, target_y)` rows with optional species filtering,
/// multi-task column selection (`spec.target_col`), and `log10` transformation.
pub fn parse_tdc_tab(path: &Path, spec: &TdcEndpointSpec) -> Result<Vec<(String, f64)>, QsarError> {
    let raw = std::fs::read_to_string(path)?;
    let mut lines = raw.lines();
    let Some(header) = lines.next() else {
        return Err(QsarError::Dataset {
            endpoint: spec.id.to_string(),
            message: "empty .tab file".to_string(),
        });
    };
    let cols: Vec<&str> = header.split('\t').map(unquote).collect();
    let smi_idx = cols
        .iter()
        .position(|&c| {
            c.eq_ignore_ascii_case("Drug")
                || c.eq_ignore_ascii_case("X")
                || c.eq_ignore_ascii_case("SMILES")
        })
        .unwrap_or(1);
    let y_idx = if let Some(req_col) = spec.target_col {
        cols.iter()
            .position(|&c| c.eq_ignore_ascii_case(req_col))
            .ok_or_else(|| QsarError::Dataset {
                endpoint: spec.id.to_string(),
                message: format!("required multi-task column '{req_col}' not found in header"),
            })?
    } else {
        cols.iter()
            .position(|&c| c.eq_ignore_ascii_case("Y"))
            .unwrap_or(2)
    };
    let sp_idx = cols.iter().position(|&c| c.eq_ignore_ascii_case("Species"));

    let mut out = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= smi_idx.max(y_idx) {
            continue;
        }
        if let (Some(req_sp), Some(si)) = (spec.species_filter, sp_idx) {
            let got = fields.get(si).map(|s| unquote(s)).unwrap_or("");
            if !got.eq_ignore_ascii_case(req_sp) {
                continue;
            }
        }
        let smi = unquote(fields[smi_idx]);
        if smi.is_empty() {
            continue;
        }
        let raw_y = unquote(fields[y_idx]);
        if raw_y.is_empty() || raw_y.eq_ignore_ascii_case("nan") || raw_y.eq_ignore_ascii_case("na")
        {
            continue;
        }
        let Ok(mut y) = raw_y.parse::<f64>() else {
            continue;
        };
        if !y.is_finite() {
            continue;
        }
        if spec.log10_transform {
            if y <= 0.0 {
                continue;
            }
            y = y.log10();
        }
        out.push((smi.to_string(), y));
    }
    Ok(out)
}

/// Deterministic Bemis-Murcko 70% Train / 10% Cal / 20% Test scaffold split with ZERO scaffold leakage.
pub fn scaffold_split_indices(
    feats: &[MolFeatures],
    ys: &[f64],
    is_binary: bool,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let n = feats.len();
    let mut groups: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, f) in feats.iter().enumerate() {
        groups.entry(f.scaffold_hash).or_default().push(i);
    }

    // Sort scaffold groups descending by size (standard DeepChem / TDC scaffold split),
    // breaking ties deterministically by scaffold hash.
    let mut buckets: Vec<(u64, Vec<usize>)> = groups.into_iter().collect();
    buckets.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));

    let target_train = ((n as f64) * 0.70).round() as usize;
    let target_cal = ((n as f64) * 0.10).round().max(10.0) as usize;

    let mut train = Vec::with_capacity(target_train + 32);
    let mut cal = Vec::with_capacity(target_cal + 32);
    let mut test = Vec::with_capacity(n / 4 + 32);

    for (_, idxs) in buckets {
        // Assign entire scaffold bucket atomically (never split a scaffold across folds!)
        if train.len() + idxs.len() <= target_train || train.is_empty() {
            train.extend(idxs);
        } else if cal.len() + idxs.len() <= target_cal || cal.is_empty() {
            cal.extend(idxs);
        } else {
            test.extend(idxs);
        }
    }

    // Safety check for tiny imbalanced binary datasets: if cal or test ended up single-class,
    // swap one small opposite-class scaffold bucket from train while preserving zero scaffold overlap.
    if is_binary {
        for fold in [&mut cal, &mut test] {
            let has_pos = fold.iter().any(|&i| ys[i] > 0.5);
            let has_neg = fold.iter().any(|&i| ys[i] <= 0.5);
            if !has_pos || !has_neg {
                let need_pos = !has_pos;
                if let Some(&cand_idx) = train.iter().rev().find(|&&i| (ys[i] > 0.5) == need_pos) {
                    let scaf = feats[cand_idx].scaffold_hash;
                    let moved: Vec<usize> = train
                        .iter()
                        .copied()
                        .filter(|&i| feats[i].scaffold_hash == scaf)
                        .collect();
                    if moved.len() * 3 < train.len() {
                        train.retain(|&i| feats[i].scaffold_hash != scaf);
                        fold.extend(moved);
                    }
                }
            }
        }
    }

    (train, cal, test)
}

/// Trains all available TDC ADMET endpoints found in `tdc_dir`, writes the compact model pack
/// to `out_pack_path`, and optionally writes a Markdown validation report to `out_report_md`.
pub fn train_tdc_model_pack(
    tdc_dir: &Path,
    out_pack_path: &Path,
    out_report_md: Option<&Path>,
    cfg: &GbdtTrainConfig,
) -> Result<AdmetModelPack, QsarError> {
    let t_all = Instant::now();
    let mut trained_models = Vec::with_capacity(TDC_ENDPOINTS.len());

    for (idx, spec) in TDC_ENDPOINTS.iter().enumerate() {
        let tab_path = tdc_dir.join(spec.file_name);
        if !tab_path.is_file() {
            continue;
        }
        let t0 = Instant::now();
        let raw_rows = parse_tdc_tab(&tab_path, spec)?;

        // Parallel SMILES parsing & feature extraction across all CPU threads
        let parsed: Vec<(MolFeatures, f64)> = raw_rows
            .par_iter()
            .filter_map(|(smi, y)| {
                let mol = parse_smiles(smi).ok()?;
                if mol.atoms.is_empty() {
                    return None;
                }
                let eval = evaluate_molecule(&mol);
                let feat = extract_mol_features(&mol, &eval);
                Some((feat, *y))
            })
            .collect();

        if parsed.len() < 40 {
            continue;
        }

        let feats: Vec<MolFeatures> = parsed.iter().map(|(f, _)| f.clone()).collect();
        let ys: Vec<f64> = parsed.iter().map(|(_, y)| *y).collect();
        let is_binary = spec.task == QsarTaskKind::BinaryClassification;
        let (tr_idx, cal_idx, te_idx) = scaffold_split_indices(&feats, &ys, is_binary);

        let tr_f: Vec<&MolFeatures> = tr_idx.iter().map(|&i| &feats[i]).collect();
        let tr_y: Vec<f64> = tr_idx.iter().map(|&i| ys[i]).collect();
        let cal_f: Vec<&MolFeatures> = cal_idx.iter().map(|&i| &feats[i]).collect();
        let cal_y: Vec<f64> = cal_idx.iter().map(|&i| ys[i]).collect();
        let te_f: Vec<&MolFeatures> = te_idx.iter().map(|&i| &feats[i]).collect();
        let te_y: Vec<f64> = te_idx.iter().map(|&i| ys[i]).collect();

        let mut ep_cfg = cfg.clone();
        ep_cfg.seed = cfg.seed.wrapping_add((idx as u64 + 1).wrapping_mul(1009));
        // Smaller endpoints benefit from gentler depth & min_child_hess to avoid overfitting
        if tr_f.len() < 800 {
            ep_cfg.max_depth = ep_cfg.max_depth.min(4);
            ep_cfg.min_child_hess = 2.0;
            ep_cfg.l2_reg = 4.0;
        }

        let model = train_endpoint_gbdt(
            spec.id,
            spec.display_name,
            spec.category,
            spec.task,
            spec.units,
            spec.log10_transform,
            &tr_f,
            &tr_y,
            &cal_f,
            &cal_y,
            &te_f,
            &te_y,
            &ep_cfg,
        );

        let metric_summary = match spec.task {
            QsarTaskKind::BinaryClassification => format!(
                "AUC={:.3} PR={:.3} Cov90={:.1}%",
                model.metrics.roc_auc.unwrap_or(0.5),
                model.metrics.pr_auc.unwrap_or(0.0),
                model.metrics.conformal_90_coverage * 100.0
            ),
            QsarTaskKind::Regression => format!(
                "rho={:.3} RMSE={:.3} MAE={:.3} Cov90={:.1}%",
                model.metrics.spearman_rho.unwrap_or(0.0),
                model.metrics.rmse.unwrap_or(0.0),
                model.metrics.mae.unwrap_or(0.0),
                model.metrics.conformal_90_coverage * 100.0
            ),
        };
        eprintln!(
            "[{:2}/{}] {:25} n={:5} ({}/{}/{}) trees={:3} | {} ({} ms)",
            idx + 1,
            TDC_ENDPOINTS.len(),
            spec.id,
            feats.len(),
            tr_f.len(),
            cal_f.len(),
            te_f.len(),
            model.metrics.num_trees_used,
            metric_summary,
            t0.elapsed().as_millis()
        );
        trained_models.push(model);
    }

    let pack = AdmetModelPack {
        format: "admetcraft-qsar-pack".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        trained_at: format!("local-rust-gbdt-{}s", t_all.elapsed().as_secs()),
        split_method: "Bemis-Murcko Scaffold Split 70/10/20 + Split Conformal (90%)".to_string(),
        models: trained_models,
    };

    pack.save_to_path(out_pack_path)?;
    if let Some(md_path) = out_report_md {
        if let Some(parent) = md_path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(md_path, pack.render_validation_markdown())?;
    }
    Ok(pack)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_tdc_tab_reads_multitask_tox21_columns_and_skips_missing_cells()
    -> Result<(), Box<dyn std::error::Error>> {
        let p = std::env::temp_dir().join("admetcraft_tox21_parse_test.tab");
        let content = "Drug_ID\tDrug\tNR-AR\tSR-MMP\n\
                       c1\tCCO\t0\t1\n\
                       c2\tc1ccccc1\t\t0\n\
                       c3\tCC(=O)O\tNaN\t1\n\
                       c4\tNc1ccccc1\t1\tNA\n";
        std::fs::write(&p, content)?;

        let nr_ar = TDC_ENDPOINTS
            .iter()
            .find(|e| e.id == "tox21_nr_ar")
            .ok_or("missing tox21_nr_ar spec")?;
        let sr_mmp = TDC_ENDPOINTS
            .iter()
            .find(|e| e.id == "tox21_sr_mmp")
            .ok_or("missing tox21_sr_mmp spec")?;

        let rows_ar = parse_tdc_tab(&p, nr_ar)?;
        let rows_mmp = parse_tdc_tab(&p, sr_mmp)?;
        let _ = std::fs::remove_file(&p);

        assert_eq!(rows_ar.len(), 2);
        assert_eq!(rows_ar[0].0, "CCO");
        assert!((rows_ar[0].1 - 0.0).abs() < 1e-9);
        assert_eq!(rows_ar[1].0, "Nc1ccccc1");
        assert!((rows_ar[1].1 - 1.0).abs() < 1e-9);

        assert_eq!(rows_mmp.len(), 3);
        Ok(())
    }
}
