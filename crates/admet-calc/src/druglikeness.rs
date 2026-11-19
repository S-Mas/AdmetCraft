use chem_core::Molecule;

use crate::pharmacokinetics::has_carboxylic_or_sulfonic_acid;
use crate::physicochem::{count_n_or_o, count_nh_or_oh, round2};
use crate::types::{
    Druglikeness, IonizationProfile, Lipophilicity, PhysicochemicalProperties, RuleResult,
};

/// Evaluate Lipinski, Ghose, Veber, Egan, Muegge, Abbott Bioavailability Score,
/// Bickerton (2012) QED, Pfizer CNS MPO (Wager 2010), Pfizer 3/75, GSK 4/400,
/// Golden Triangle, Rule of Three (Ro3), and Beyond Rule of 5 (bRo5).
pub fn compute_druglikeness(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
    ion: &IonizationProfile,
    brenk_alert_count: usize,
) -> Druglikeness {
    let lipinski = eval_lipinski(mol, phys, lipo);
    let ghose = eval_ghose(mol, phys, lipo);
    let veber = eval_veber(phys);
    let egan = eval_egan(phys, lipo);
    let muegge = eval_muegge(mol, phys, lipo);
    let bioavailability_score = compute_abbott_bioavailability(mol, phys, &lipinski);
    let qed = compute_qed(phys, lipo, brenk_alert_count);
    let cns_mpo = compute_cns_mpo(phys, lipo, ion);
    let pfizer_3_75 = eval_pfizer_3_75(phys, lipo);
    let gsk_4_400 = eval_gsk_4_400(phys, lipo);
    let golden_triangle = eval_golden_triangle(phys, ion);
    let ro3_fragment = eval_ro3_fragment(phys, lipo);
    let beyond_ro5 = phys.molecular_weight > 500.0
        || phys.tpsa > 140.0
        || phys.h_bond_donors > 5
        || phys.rotatable_bonds > 10;

    Druglikeness {
        lipinski,
        ghose,
        veber,
        egan,
        muegge,
        bioavailability_score,
        qed,
        cns_mpo,
        pfizer_3_75,
        gsk_4_400,
        golden_triangle,
        ro3_fragment,
        beyond_ro5,
    }
}

/// Exact Bickerton et al. (2012, *Nature Chemistry* 4, 90-98) Quantitative Estimate of Drug-likeness
/// (`QED_w`, weighted geometric mean of 8 asymmetric double-sigmoidal desirability curves).
pub fn compute_qed(
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
    structural_alerts: usize,
) -> f64 {
    #[allow(clippy::too_many_arguments)]
    fn ads(x: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, dmax: f64) -> f64 {
        let s1 = 1.0 / (1.0 + (-(x - c + d * 0.5) / e).exp());
        let s2 = 1.0 - 1.0 / (1.0 + (-(x - c - d * 0.5) / f).exp());
        ((a + b * s1 * s2) / dmax).clamp(1e-6, 1.0)
    }

    let d_mw = ads(
        phys.molecular_weight,
        2.817_065_973,
        392.575_495_3,
        290.748_976_4,
        2.419_764_353,
        49.223_256_77,
        65.370_517_07,
        104.980_556_1,
    );
    let d_alogp = ads(
        lipo.wlogp,
        3.172_690_585,
        137.862_475_1,
        2.534_937_431,
        4.581_497_897,
        0.822_739_154,
        0.576_295_591,
        131.318_660_4,
    );
    let d_hba = ads(
        phys.h_bond_acceptors as f64,
        2.948_620_388,
        160.460_597_2,
        3.615_294_657,
        4.435_986_202,
        0.290_141_953,
        1.300_669_958,
        148.776_304_6,
    );
    let d_hbd = ads(
        phys.h_bond_donors as f64,
        1.618_662_227,
        1_010.051_101,
        0.985_094_388,
        1e-9,
        0.713_820_843,
        0.920_922_555,
        258.163_261_6,
    );
    let d_psa = ads(
        phys.tpsa_no_only,
        1.876_861_559,
        125.223_265_7,
        62.907_735_54,
        87.833_666_14,
        12.019_998_24,
        28.513_247_32,
        104.568_616_7,
    );
    let d_rotb = ads(
        phys.rotatable_bonds as f64,
        0.01,
        272.412_142_7,
        2.558_379_97,
        1.565_547_684,
        1.271_567_166,
        2.758_063_707,
        105.442_040_3,
    );
    let d_arom = ads(
        phys.num_aromatic_rings as f64,
        3.217_788_97,
        957.737_410_8,
        2.274_627_939,
        1e-9,
        1.317_690_384,
        0.375_760_881,
        312.337_261,
    );
    let d_alerts = ads(
        structural_alerts as f64,
        0.01,
        1_199.094_025,
        -0.090_028_83,
        1e-9,
        0.185_904_477,
        0.875_193_782,
        417.725_314,
    );

    let weights = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95];
    let ds = [d_mw, d_alogp, d_hba, d_hbd, d_psa, d_rotb, d_arom, d_alerts];
    let w_sum: f64 = weights.iter().sum();
    let log_sum: f64 = weights.iter().zip(ds.iter()).map(|(w, d)| w * d.ln()).sum();
    let raw = (log_sum / w_sum).exp().clamp(0.0, 1.0);
    (raw * 1000.0).round() / 1000.0
}

/// Wager et al. (2010, *ACS Chem. Neurosci.* 1, 435-449) Pfizer CNS MPO score `[0.00, 6.00]`.
pub fn compute_cns_mpo(
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
    ion: &IonizationProfile,
) -> f64 {
    let mono = |v: f64, lo: f64, hi: f64| -> f64 {
        if v <= lo {
            1.0
        } else if v >= hi {
            0.0
        } else {
            (hi - v) / (hi - lo)
        }
    };
    let t_clogp = mono(lipo.consensus_log_p, 3.0, 5.0);
    let t_clogd = mono(ion.logd_ph74, 2.0, 4.0);
    let t_mw = mono(phys.molecular_weight, 360.0, 500.0);
    let t_tpsa = if phys.tpsa <= 20.0 || phys.tpsa >= 120.0 {
        0.0
    } else if phys.tpsa < 40.0 {
        (phys.tpsa - 20.0) / 20.0
    } else if phys.tpsa <= 90.0 {
        1.0
    } else {
        (120.0 - phys.tpsa) / 30.0
    };
    let t_hbd = if phys.h_bond_donors == 0 {
        1.0
    } else {
        ((3.5 - phys.h_bond_donors as f64) / 3.0).clamp(0.0, 1.0)
    };
    let t_pka = ion
        .strongest_basic_pka
        .map_or(1.0, |pka| mono(pka, 8.0, 10.0));
    round2((t_clogp + t_clogd + t_mw + t_tpsa + t_hbd + t_pka).clamp(0.0, 6.0))
}

fn eval_pfizer_3_75(phys: &PhysicochemicalProperties, lipo: &Lipophilicity) -> RuleResult {
    let mut details = Vec::new();
    if lipo.consensus_log_p > 3.0 && phys.tpsa < 75.0 {
        details.push(format!(
            "cLogP>3 ({:.2}) & TPSA<75 ({:.1}) — elevated in vivo toxicity risk",
            lipo.consensus_log_p, phys.tpsa
        ));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_gsk_4_400(phys: &PhysicochemicalProperties, lipo: &Lipophilicity) -> RuleResult {
    let mut details = Vec::new();
    if phys.molecular_weight > 400.0 {
        details.push(format!("MW>400 ({:.1})", phys.molecular_weight));
    }
    if lipo.consensus_log_p > 4.0 {
        details.push(format!("cLogP>4 ({:.2})", lipo.consensus_log_p));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_golden_triangle(phys: &PhysicochemicalProperties, ion: &IonizationProfile) -> RuleResult {
    let mut details = Vec::new();
    if phys.molecular_weight < 200.0 || phys.molecular_weight > 500.0 {
        details.push(format!(
            "MW outside [200, 500] ({:.1})",
            phys.molecular_weight
        ));
    }
    if ion.logd_ph74 < -2.0 || ion.logd_ph74 > 5.0 {
        details.push(format!("LogD7.4 outside [-2, 5] ({:.2})", ion.logd_ph74));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_ro3_fragment(phys: &PhysicochemicalProperties, lipo: &Lipophilicity) -> RuleResult {
    let mut details = Vec::new();
    if phys.molecular_weight > 300.0 {
        details.push(format!("MW>300 ({:.1})", phys.molecular_weight));
    }
    if lipo.consensus_log_p > 3.0 {
        details.push(format!("cLogP>3 ({:.2})", lipo.consensus_log_p));
    }
    if phys.h_bond_donors > 3 {
        details.push(format!("HBD>3 ({})", phys.h_bond_donors));
    }
    if phys.h_bond_acceptors > 3 {
        details.push(format!("HBA>3 ({})", phys.h_bond_acceptors));
    }
    if phys.rotatable_bonds > 3 {
        details.push(format!("Rotors>3 ({})", phys.rotatable_bonds));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_lipinski(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> RuleResult {
    let n_or_o = count_n_or_o(mol);
    let nh_or_oh = count_nh_or_oh(mol);
    let mut details = Vec::new();
    if phys.molecular_weight > 500.0 {
        details.push(format!("MW>500 ({:.2})", phys.molecular_weight));
    }
    if lipo.mlogp > 4.15 {
        details.push(format!("MLOGP>4.15 ({:.2})", lipo.mlogp));
    }
    if n_or_o > 10 {
        details.push(format!("NorO>10 ({n_or_o})"));
    }
    if nh_or_oh > 5 {
        details.push(format!("NHorOH>5 ({nh_or_oh})"));
    }
    let violations = details.len();
    RuleResult {
        pass: violations <= 1,
        violations,
        details,
    }
}

fn eval_ghose(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> RuleResult {
    let total_atoms = mol.atoms.len()
        + mol
            .atoms
            .iter()
            .map(|a| usize::from(a.total_h()))
            .sum::<usize>();
    let mut details = Vec::new();
    if phys.molecular_weight < 160.0 || phys.molecular_weight > 480.0 {
        details.push(format!(
            "MW outside [160, 480] ({:.2})",
            phys.molecular_weight
        ));
    }
    if lipo.wlogp < -0.4 || lipo.wlogp > 5.6 {
        details.push(format!("WLOGP outside [-0.4, 5.6] ({:.2})", lipo.wlogp));
    }
    if phys.molar_refractivity < 40.0 || phys.molar_refractivity > 130.0 {
        details.push(format!(
            "MR outside [40, 130] ({:.2})",
            phys.molar_refractivity
        ));
    }
    if !(20..=70).contains(&total_atoms) {
        details.push(format!("#atoms outside [20, 70] ({total_atoms})"));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_veber(phys: &PhysicochemicalProperties) -> RuleResult {
    let mut details = Vec::new();
    if phys.rotatable_bonds > 10 {
        details.push(format!("Rotors>10 ({})", phys.rotatable_bonds));
    }
    if phys.tpsa > 140.0 {
        details.push(format!("TPSA>140 ({:.2})", phys.tpsa));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_egan(phys: &PhysicochemicalProperties, lipo: &Lipophilicity) -> RuleResult {
    let mut details = Vec::new();
    if lipo.wlogp > 5.88 {
        details.push(format!("WLOGP>5.88 ({:.2})", lipo.wlogp));
    }
    if phys.tpsa > 131.6 {
        details.push(format!("TPSA>131.6 ({:.2})", phys.tpsa));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

fn eval_muegge(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> RuleResult {
    let num_carbons = mol.atoms.iter().filter(|a| a.atomic_number == 6).count();
    let num_hetero = mol
        .atoms
        .iter()
        .filter(|a| !matches!(a.atomic_number, 1 | 6))
        .count();
    let num_rings = mol.rings.len();
    let mut details = Vec::new();
    if phys.molecular_weight < 200.0 || phys.molecular_weight > 600.0 {
        details.push(format!(
            "MW outside [200, 600] ({:.2})",
            phys.molecular_weight
        ));
    }
    if lipo.xlogp3 < -2.0 || lipo.xlogp3 > 5.0 {
        details.push(format!("XLOGP3 outside [-2, 5] ({:.2})", lipo.xlogp3));
    }
    if phys.tpsa > 150.0 {
        details.push(format!("TPSA>150 ({:.2})", phys.tpsa));
    }
    if num_rings > 7 {
        details.push(format!("Num. rings>7 ({num_rings})"));
    }
    if num_carbons <= 4 {
        details.push(format!("Num. carbons<=4 ({num_carbons})"));
    }
    if num_hetero <= 1 {
        details.push(format!("Num. heteroatoms<=1 ({num_hetero})"));
    }
    if phys.rotatable_bonds > 15 {
        details.push(format!(
            "Num. rotatable bonds>15 ({})",
            phys.rotatable_bonds
        ));
    }
    if phys.h_bond_acceptors > 10 {
        details.push(format!("HBA>10 ({})", phys.h_bond_acceptors));
    }
    if phys.h_bond_donors > 5 {
        details.push(format!("HBD>5 ({})", phys.h_bond_donors));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

/// Martin (2005) *J. Med. Chem.* 48, 3164-3170 Abbott Bioavailability Score.
fn compute_abbott_bioavailability(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipinski: &RuleResult,
) -> f64 {
    let net_charge: i32 = mol.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();
    let is_anion = net_charge < 0 || has_carboxylic_or_sulfonic_acid(mol);
    if is_anion {
        if phys.tpsa <= 75.0 {
            0.85
        } else if phys.tpsa <= 150.0 {
            0.56
        } else {
            0.11
        }
    } else if lipinski.pass {
        0.55
    } else {
        0.17
    }
}
