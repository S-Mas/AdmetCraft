use serde::{Deserialize, Serialize};

/// Gastrointestinal passive absorption class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AbsorptionLevel {
    High,
    Low,
}

impl AbsorptionLevel {
    /// SwissADME display label.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::High => "High",
            Self::Low => "Low",
        }
    }
}

/// Qualitative aqueous solubility class (SwissADME scale).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolubilityClass {
    Insoluble,
    PoorlySoluble,
    ModeratelySoluble,
    Soluble,
    VerySoluble,
    HighlySoluble,
}

impl SolubilityClass {
    /// Classify a `log S` (mol/L) value.
    pub fn from_log_s(log_s: f64) -> Self {
        if log_s < -10.0 {
            Self::Insoluble
        } else if log_s < -6.0 {
            Self::PoorlySoluble
        } else if log_s < -4.0 {
            Self::ModeratelySoluble
        } else if log_s < -2.0 {
            Self::Soluble
        } else if log_s < 0.0 {
            Self::VerySoluble
        } else {
            Self::HighlySoluble
        }
    }

    /// SwissADME display label.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Insoluble => "Insoluble",
            Self::PoorlySoluble => "Poorly soluble",
            Self::ModeratelySoluble => "Moderately soluble",
            Self::Soluble => "Soluble",
            Self::VerySoluble => "Very soluble",
            Self::HighlySoluble => "Highly soluble",
        }
    }
}

/// Physicochemical descriptors section (SwissADME parity).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhysicochemicalProperties {
    pub formula: String,
    pub molecular_weight: f64,
    pub heavy_atoms: usize,
    pub aromatic_heavy_atoms: usize,
    pub fraction_csp3: f64,
    pub rotatable_bonds: usize,
    pub h_bond_acceptors: usize,
    pub h_bond_donors: usize,
    pub molar_refractivity: f64,
    pub tpsa: f64,
}

/// Lipophilicity predictions across 5 models plus consensus average.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lipophilicity {
    pub ilogp: f64,
    pub xlogp3: f64,
    pub wlogp: f64,
    pub mlogp: f64,
    pub silicos_it: f64,
    pub consensus_log_p: f64,
}

/// Single water solubility model output.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterSolubilityModel {
    pub log_s: f64,
    pub solubility_mg_ml: f64,
    pub solubility_mol_l: f64,
    pub class: SolubilityClass,
}

/// Water solubility section (ESOL, Ali, Silicos-IT).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterSolubility {
    pub esol: WaterSolubilityModel,
    pub ali: WaterSolubilityModel,
    pub silicos_it: WaterSolubilityModel,
}

/// Pharmacokinetics (ADME) predictions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pharmacokinetics {
    pub gi_absorption: AbsorptionLevel,
    pub bbb_permeant: bool,
    pub pgp_substrate: bool,
    pub cyp1a2_inhibitor: bool,
    pub cyp2c19_inhibitor: bool,
    pub cyp2c9_inhibitor: bool,
    pub cyp2d6_inhibitor: bool,
    pub cyp3a4_inhibitor: bool,
    pub log_kp_cm_s: f64,
}

/// Outcome of a druglikeness / leadlikeness rule filter.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RuleResult {
    pub pass: bool,
    pub violations: usize,
    pub details: Vec<String>,
}

/// Druglikeness section (5 pharmaceutical filters + Abbott Bioavailability Score).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Druglikeness {
    pub lipinski: RuleResult,
    pub ghose: RuleResult,
    pub veber: RuleResult,
    pub egan: RuleResult,
    pub muegge: RuleResult,
    pub bioavailability_score: f64,
}

/// Medicinal Chemistry alerts and scores.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MedicinalChemistry {
    pub pains_alerts: Vec<String>,
    pub brenk_alerts: Vec<String>,
    pub leadlikeness: RuleResult,
    pub synthetic_accessibility: f64,
}

/// BOILED-Egg coordinate and region classification.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoiledEggPoint {
    pub tpsa: f64,
    pub wlogp: f64,
    pub in_white: bool,
    pub in_yolk: bool,
    pub pgp_substrate: bool,
    pub in_plot_bounds: bool,
}

/// Six-axis Bioavailability Radar descriptors and normalized values `[0, 1]` where `[0.25, 0.85]`
/// highlights the optimal pink physicochemical zone.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BioavailabilityRadar {
    pub lipo_xlogp3: f64,
    pub size_mw: f64,
    pub polar_tpsa: f64,
    pub insolu_log_s: f64,
    pub insatu_fsp3: f64,
    pub flex_rot_bonds: usize,
    /// Normalized axis scores `[0.0, 1.0]` in order `[LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX]`.
    pub normalized_axes: [f64; 6],
    /// True if all 6 axes lie inside the SwissADME pink optimal zone.
    pub all_in_optimal_zone: bool,
}

/// Complete ADMET evaluation for one molecule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MoleculeEvaluation {
    pub name: String,
    pub input_smiles: String,
    pub canonical_smiles: String,
    pub physicochemical: PhysicochemicalProperties,
    pub lipophilicity: Lipophilicity,
    pub water_solubility: WaterSolubility,
    pub pharmacokinetics: Pharmacokinetics,
    pub druglikeness: Druglikeness,
    pub medicinal_chemistry: MedicinalChemistry,
    pub boiled_egg: BoiledEggPoint,
    pub bioavailability_radar: BioavailabilityRadar,
}
