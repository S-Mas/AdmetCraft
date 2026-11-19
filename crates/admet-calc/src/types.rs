use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

/// Gastrointestinal passive absorption class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AbsorptionLevel {
    High,
    Low,
}

impl AbsorptionLevel {
    /// Display label.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::High => "High",
            Self::Low => "Low",
        }
    }
}

/// Qualitative aqueous solubility class.
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

    /// Display label.
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

/// Physicochemical descriptors section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhysicochemicalProperties {
    pub formula: String,
    pub molecular_weight: f64,
    #[serde(default)]
    pub exact_mass: f64,
    pub heavy_atoms: usize,
    pub aromatic_heavy_atoms: usize,
    pub fraction_csp3: f64,
    pub rotatable_bonds: usize,
    pub h_bond_acceptors: usize,
    pub h_bond_donors: usize,
    pub molar_refractivity: f64,
    /// Ertl (2000) TPSA including N, O, S, and P (SwissADME convention).
    pub tpsa: f64,
    /// Ertl (2000) TPSA counting N and O only (RDKit / QED convention).
    #[serde(default)]
    pub tpsa_no_only: f64,
    /// Double Bond Equivalents (Rings + Double Bonds).
    #[serde(default)]
    pub dbe: f64,
    #[serde(default)]
    pub num_rings: usize,
    #[serde(default)]
    pub num_aromatic_rings: usize,
}

/// Dominant physiological ionization state at pH 7.4.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum IonizationClass {
    #[default]
    Neutral,
    Acid,
    Base,
    Zwitterion,
}

impl IonizationClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Neutral => "Neutral",
            Self::Acid => "Acid",
            Self::Base => "Base",
            Self::Zwitterion => "Zwitterion",
        }
    }
}

/// Physiological pKa, ionization class, and Henderson-Hasselbalch LogD at pH 7.4.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct IonizationProfile {
    pub strongest_acidic_pka: Option<f64>,
    pub strongest_basic_pka: Option<f64>,
    pub ionization_class: IonizationClass,
    pub estimated_charge_ph74: f64,
    pub logd_ph74: f64,
}

/// Qualitative hERG cardiotoxicity structural/physicochemical risk class (heuristic).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum HergRisk {
    #[default]
    Low,
    Moderate,
    High,
}

impl HergRisk {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Low => "Low",
            Self::Moderate => "Moderate",
            Self::High => "High",
        }
    }
}

/// Extended toxicity & reactive-metabolite structural safety screening panel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct SafetyProfile {
    pub herg_risk: HergRisk,
    pub herg_reasons: Vec<String>,
    pub ames_alerts: Vec<String>,
    pub dili_alerts: Vec<String>,
    pub safety_alerts_are_structural_heuristics: bool,
}

/// Lipophilicity predictions: 3 distinct models (`XLOGP3`-heuristic, `WLOGP`, `MLOGP`),
/// their arithmetic mean (`consensus_log_p`), plus 2 clearly labeled experimental linear surrogates.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lipophilicity {
    /// Wildman & Crippen (1999) atomic logP.
    pub wlogp: f64,
    /// Atom-additive logP heuristic inspired by Cheng et al. (2007).
    pub xlogp3: f64,
    /// Moriguchi et al. (1992) topological logP.
    pub mlogp: f64,
    /// Experimental solvation-inspired linear heuristic (NOT SwissADME iLOGP).
    pub solvation_logp_est: f64,
    /// Experimental topological linear heuristic (NOT Silicos-IT).
    pub topological_logp_est: f64,
    /// Arithmetic mean of the 3 primary models (`wlogp`, `xlogp3`, `mlogp`).
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

/// Water solubility section (ESOL, Ali, and an experimental topological heuristic).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterSolubility {
    pub esol: WaterSolubilityModel,
    pub ali: WaterSolubilityModel,
    pub topological_est: WaterSolubilityModel,
}

/// Pharmacokinetics (ADME) predictions.
/// Note: GI/BBB use the BOILED-Egg ellipses and log Kp uses Potts & Guy (1992);
/// P-gp substrate and CYP450 flags are experimental rule-based heuristics (not trained SVMs).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pharmacokinetics {
    pub gi_absorption: AbsorptionLevel,
    pub bbb_permeant: bool,
    /// Experimental heuristic P-gp substrate flag.
    pub pgp_substrate: bool,
    /// Experimental heuristic CYP1A2 inhibitor flag.
    pub cyp1a2_inhibitor: bool,
    /// Experimental heuristic CYP2C19 inhibitor flag.
    pub cyp2c19_inhibitor: bool,
    /// Experimental heuristic CYP2C9 inhibitor flag.
    pub cyp2c9_inhibitor: bool,
    /// Experimental heuristic CYP2D6 inhibitor flag.
    pub cyp2d6_inhibitor: bool,
    /// Experimental heuristic CYP3A4 inhibitor flag.
    pub cyp3a4_inhibitor: bool,
    pub log_kp_cm_s: f64,
    pub pgp_cyp_are_experimental_heuristics: bool,
}

/// Outcome of a druglikeness / leadlikeness rule filter.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
pub struct RuleResult {
    pub pass: bool,
    pub violations: usize,
    pub details: Vec<String>,
}

/// Druglikeness section (SwissADME 5 filters + Bickerton QED + Pfizer CNS MPO + Industrial rules).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Druglikeness {
    pub lipinski: RuleResult,
    pub ghose: RuleResult,
    pub veber: RuleResult,
    pub egan: RuleResult,
    pub muegge: RuleResult,
    pub bioavailability_score: f64,
    /// Bickerton et al. (2012, *Nature Chemistry*) Quantitative Estimate of Drug-likeness `[0.0, 1.0]`.
    #[serde(default)]
    pub qed: f64,
    /// Wager et al. (2010) Pfizer CNS Multiparameter Optimization score `[0.0, 6.0]`.
    #[serde(default)]
    pub cns_mpo: f64,
    /// Hughes et al. (2008) Pfizer 3/75 in vivo toxicity rule (`pass = !(LogP > 3 && TPSA < 75)`).
    #[serde(default)]
    pub pfizer_3_75: RuleResult,
    /// Gleeson (2008) GSK 4/400 favorable ADMET space (`MW <= 400 && LogP <= 4`).
    #[serde(default)]
    pub gsk_4_400: RuleResult,
    /// Johnson et al. (2009) Golden Triangle (`200 <= MW <= 500 && -2 <= LogD7.4 <= 5`).
    #[serde(default)]
    pub golden_triangle: RuleResult,
    /// Congreve et al. (2003) Fragment Rule of Three (`MW <= 300, LogP <= 3, HBD <= 3, HBA <= 3, RB <= 3`).
    #[serde(default)]
    pub ro3_fragment: RuleResult,
    /// Doak et al. (2014) Beyond Rule of 5 (`MW > 500` or `TPSA > 140` or `HBD > 5` or `RB > 10`).
    #[serde(default)]
    pub beyond_ro5: bool,
}

/// Medicinal Chemistry alerts and scores.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MedicinalChemistry {
    pub pains_alerts: Vec<String>,
    /// Number of curated PAINS SMARTS checked (subset of Baell & Holloway ~480 catalog).
    pub pains_catalog_size: usize,
    pub brenk_alerts: Vec<String>,
    /// Number of curated Brenk SMARTS checked (subset of Brenk et al. ~105 catalog).
    pub brenk_catalog_size: usize,
    pub leadlikeness: RuleResult,
    /// Heuristic Synthetic Accessibility proxy in `[1.0, 10.0]` (does not include PubChem fragment table).
    pub synthetic_accessibility: f64,
    pub sa_is_heuristic_proxy: bool,
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

/// Shared geometry sampler for BOILED-Egg ellipses across PDF, SVG, and GUI renderers.
pub struct BoiledEggGeometry;

impl BoiledEggGeometry {
    /// Sample points in `(tpsa, wlogp)` coordinates along a rotated ellipse.
    pub fn sample_ellipse(
        center: (f64, f64),
        axes: (f64, f64),
        angle_rad: f64,
        steps: usize,
    ) -> Vec<(f64, f64)> {
        let cos_a = angle_rad.cos();
        let sin_a = angle_rad.sin();
        (0..steps)
            .map(|i| {
                let t = 2.0 * PI * (i as f64) / (steps as f64);
                let lx = axes.0 * t.cos();
                let ly = axes.1 * t.sin();
                (
                    center.0 + lx * cos_a - ly * sin_a,
                    center.1 + lx * sin_a + ly * cos_a,
                )
            })
            .collect()
    }
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
    /// True if all 6 axes lie inside the optimal physicochemical zone.
    pub all_in_optimal_zone: bool,
}

impl BioavailabilityRadar {
    /// Axis labels in clockwise order starting at 12 o'clock.
    pub const AXIS_LABELS: [&'static str; 6] =
        ["LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"];

    /// Compute normalized unit-circle `(x, y)` coordinates (`x` right, `y` up) for a hexagon at `radius`.
    pub fn unit_hexagon(radius: f64) -> [(f64, f64); 6] {
        let mut pts = [(0.0, 0.0); 6];
        for (i, slot) in pts.iter_mut().enumerate() {
            let ang = PI / 2.0 - (i as f64) * PI / 3.0;
            *slot = (radius * ang.cos(), radius * ang.sin());
        }
        pts
    }

    /// Compute normalized `(x, y)` coordinates (`x` right, `y` up) for this molecule's radar polygon.
    pub fn unit_polygon(&self) -> [(f64, f64); 6] {
        let mut pts = [(0.0, 0.0); 6];
        for (i, (&val, slot)) in self.normalized_axes.iter().zip(pts.iter_mut()).enumerate() {
            let ang = PI / 2.0 - (i as f64) * PI / 3.0;
            let r = val.clamp(0.05, 1.0);
            *slot = (r * ang.cos(), r * ang.sin());
        }
        pts
    }
}

/// Complete ADMET evaluation for one molecule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MoleculeEvaluation {
    pub name: String,
    pub input_smiles: String,
    pub canonical_smiles: String,
    #[serde(default)]
    pub warnings: Vec<String>,
    pub physicochemical: PhysicochemicalProperties,
    pub lipophilicity: Lipophilicity,
    pub water_solubility: WaterSolubility,
    pub pharmacokinetics: Pharmacokinetics,
    pub druglikeness: Druglikeness,
    pub medicinal_chemistry: MedicinalChemistry,
    #[serde(default)]
    pub ionization: IonizationProfile,
    #[serde(default)]
    pub safety: SafetyProfile,
    pub boiled_egg: BoiledEggPoint,
    pub bioavailability_radar: BioavailabilityRadar,
}

/// Side-by-side ADMET delta between two molecules (e.g. parent vs analogue / bioisostere).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdmetComparison {
    pub name_a: String,
    pub name_b: String,
    pub smiles_a: String,
    pub smiles_b: String,
    pub delta_mw: f64,
    pub delta_tpsa: f64,
    pub delta_wlogp: f64,
    pub delta_consensus_logp: f64,
    pub delta_logd_ph74: f64,
    pub delta_esol_log_s: f64,
    pub delta_qed: f64,
    pub delta_cns_mpo: f64,
    pub delta_sa: f64,
    pub alerts_added: Vec<String>,
    pub alerts_removed: Vec<String>,
    pub eval_a: MoleculeEvaluation,
    pub eval_b: MoleculeEvaluation,
}
