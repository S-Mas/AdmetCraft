//! Layer 1 (`admet-calc`): Pure-Rust parallel ADMET, physicochemical, lipophilicity,
//! solubility, pharmacokinetics (BOILED-Egg), druglikeness, and medicinal chemistry calculators.

#![forbid(unsafe_code)]

pub mod druglikeness;
pub mod ionization;
pub mod lipophilicity;
pub mod medchem;
pub mod pharmacokinetics;
pub mod physicochem;
pub mod solubility;
pub mod types;

use std::collections::BTreeSet;

use chem_core::{ChemError, Molecule, parse_smiles};
use rayon::prelude::*;
use thiserror::Error;

pub use types::*;

/// Errors returned by `admet-calc`.
#[derive(Debug, Error)]
pub enum AdmetError {
    #[error("chemistry parsing error: {0}")]
    Chem(#[from] ChemError),
}

/// Evaluate a parsed `Molecule` and produce the full ADMET profile.
pub fn evaluate_molecule(mol: &Molecule) -> MoleculeEvaluation {
    let name = mol.name.clone().unwrap_or_else(|| "Molecule 1".to_string());
    let mut warnings = mol.warnings.clone();
    let canonical_smiles = match mol.try_to_canonical_smiles() {
        Ok(s) => s,
        Err(e) => {
            warnings.push(format!("Canonical SMILES omitted: {e}"));
            String::new()
        }
    };
    let input_smiles = mol
        .input_smiles
        .clone()
        .unwrap_or_else(|| canonical_smiles.clone());

    let physicochemical = physicochem::compute_physicochemical(mol);
    let lipophilicity = lipophilicity::compute_lipophilicity(mol);
    let ionization = ionization::compute_ionization(mol, &lipophilicity);
    let water_solubility =
        solubility::compute_water_solubility(mol, &physicochemical, &lipophilicity);
    let (pharmacokinetics, boiled_egg) =
        pharmacokinetics::compute_pharmacokinetics(mol, &physicochemical, &lipophilicity);
    let medicinal_chemistry =
        medchem::compute_medicinal_chemistry(mol, &physicochemical, &lipophilicity);
    let druglikeness = druglikeness::compute_druglikeness(
        mol,
        &physicochemical,
        &lipophilicity,
        &ionization,
        medicinal_chemistry.brenk_alerts.len(),
    );
    let safety =
        medchem::compute_safety_profile(mol, &physicochemical, &lipophilicity, &ionization);
    let bioavailability_radar = compute_radar(&physicochemical, &lipophilicity, &water_solubility);

    MoleculeEvaluation {
        name,
        input_smiles,
        canonical_smiles,
        warnings,
        physicochemical,
        lipophilicity,
        water_solubility,
        pharmacokinetics,
        druglikeness,
        medicinal_chemistry,
        ionization,
        safety,
        boiled_egg,
        bioavailability_radar,
        qsar: None,
        metabolism_sites: Vec::new(),
        tpp_mpo: None,
        tox21_summary: None,
        off_targets: None,
        metabolite_guard: None,
    }
}

/// Compare two evaluated molecules side by side (delta ADMET, QED, CNS MPO, LogD7.4, and alert diff).
pub fn compare_evaluations(a: MoleculeEvaluation, b: MoleculeEvaluation) -> AdmetComparison {
    let r2 = physicochem::round2;
    let r3 = |v: f64| (v * 1000.0).round() / 1000.0;
    let collect_alerts = |ev: &MoleculeEvaluation| -> BTreeSet<String> {
        let mut s = BTreeSet::new();
        for x in &ev.medicinal_chemistry.pains_alerts {
            s.insert(format!("PAINS:{x}"));
        }
        for x in &ev.medicinal_chemistry.brenk_alerts {
            s.insert(format!("Brenk:{x}"));
        }
        for x in &ev.safety.ames_alerts {
            s.insert(format!("AMES:{x}"));
        }
        for x in &ev.safety.dili_alerts {
            s.insert(format!("DILI:{x}"));
        }
        s
    };
    let set_a = collect_alerts(&a);
    let set_b = collect_alerts(&b);
    let alerts_added = set_b.difference(&set_a).cloned().collect();
    let alerts_removed = set_a.difference(&set_b).cloned().collect();

    AdmetComparison {
        name_a: a.name.clone(),
        name_b: b.name.clone(),
        smiles_a: a.canonical_smiles.clone(),
        smiles_b: b.canonical_smiles.clone(),
        delta_mw: r2(b.physicochemical.molecular_weight - a.physicochemical.molecular_weight),
        delta_tpsa: r2(b.physicochemical.tpsa - a.physicochemical.tpsa),
        delta_wlogp: r2(b.lipophilicity.wlogp - a.lipophilicity.wlogp),
        delta_consensus_logp: r2(b.lipophilicity.consensus_log_p - a.lipophilicity.consensus_log_p),
        delta_logd_ph74: r2(b.ionization.logd_ph74 - a.ionization.logd_ph74),
        delta_esol_log_s: r2(b.water_solubility.esol.log_s - a.water_solubility.esol.log_s),
        delta_qed: r3(b.druglikeness.qed - a.druglikeness.qed),
        delta_cns_mpo: r2(b.druglikeness.cns_mpo - a.druglikeness.cns_mpo),
        delta_sa: r2(b.medicinal_chemistry.synthetic_accessibility
            - a.medicinal_chemistry.synthetic_accessibility),
        alerts_added,
        alerts_removed,
        eval_a: a,
        eval_b: b,
    }
}

/// Evaluate and compare two SMILES strings side by side.
pub fn compare_smiles(
    smiles_a: &str,
    name_a: Option<&str>,
    smiles_b: &str,
    name_b: Option<&str>,
) -> Result<AdmetComparison, AdmetError> {
    let a = evaluate_smiles(smiles_a, name_a)?;
    let b = evaluate_smiles(smiles_b, name_b)?;
    Ok(compare_evaluations(a, b))
}

/// Parse a SMILES string and evaluate all ADMET properties.
pub fn evaluate_smiles(smiles: &str, name: Option<&str>) -> Result<MoleculeEvaluation, AdmetError> {
    let mut mol = parse_smiles(smiles)?;
    if let Some(n) = name {
        mol.name = Some(n.to_string());
    }
    Ok(evaluate_molecule(&mol))
}

/// Evaluate a batch of molecules in parallel across all CPU cores using `rayon`.
pub fn evaluate_batch(molecules: &[Molecule]) -> Vec<MoleculeEvaluation> {
    molecules.par_iter().map(evaluate_molecule).collect()
}

/// Compute 6-axis Bioavailability Radar values (SwissADME pink zone:
/// LIPO: -0.7..5.0, SIZE: 150..500, POLAR: 20..130, INSOLU: -6..0, INSATU: >=0.25, FLEX: <=9).
fn compute_radar(
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
    sol: &WaterSolubility,
) -> BioavailabilityRadar {
    let xlogp3 = lipo.xlogp3;
    let mw = phys.molecular_weight;
    let tpsa = phys.tpsa;
    let log_s = sol.esol.log_s;
    let fsp3 = phys.fraction_csp3;
    let rb = phys.rotatable_bonds;

    // Normalize so optimal range maps to [0.25, 0.85] (matching SwissADME pink radar ring!)
    let norm_lipo = map_band(xlogp3, -0.7, 5.0, -2.5, 7.0);
    let norm_size = map_band(mw, 150.0, 500.0, 50.0, 750.0);
    let norm_polar = map_band(tpsa, 20.0, 130.0, 0.0, 200.0);
    // INSOLU axis plots -log S (0 = soluble, 6 = limit of pink zone)
    let norm_insolu = map_band(-log_s, 0.0, 6.0, -1.0, 10.0);
    // INSATU: Fraction Csp3 >= 0.25 is optimal (map 0.25..1.0 into 0.25..0.85)
    let norm_insatu = if fsp3 < 0.25 {
        (fsp3 / 0.25) * 0.25
    } else {
        0.25 + ((fsp3 - 0.25) / 0.75) * 0.60
    }
    .clamp(0.0, 1.0);
    // FLEX: rotatable bonds <= 9 is optimal (0..9 -> 0.25..0.85, >9 -> 0.85..1.0)
    let norm_flex = if rb <= 9 {
        0.25 + (rb as f64 / 9.0) * 0.60
    } else {
        0.85 + (((rb - 9) as f64) / 8.0) * 0.15
    }
    .clamp(0.0, 1.0);

    let all_in_optimal_zone = (-0.7..=5.0).contains(&xlogp3)
        && (150.0..=500.0).contains(&mw)
        && (20.0..=130.0).contains(&tpsa)
        && (-6.0..=0.0).contains(&log_s)
        && fsp3 >= 0.25
        && rb <= 9;

    BioavailabilityRadar {
        lipo_xlogp3: xlogp3,
        size_mw: mw,
        polar_tpsa: tpsa,
        insolu_log_s: log_s,
        insatu_fsp3: fsp3,
        flex_rot_bonds: rb,
        normalized_axes: [
            norm_lipo,
            norm_size,
            norm_polar,
            norm_insolu,
            norm_insatu,
            norm_flex,
        ],
        all_in_optimal_zone,
    }
}

fn map_band(val: f64, opt_lo: f64, opt_hi: f64, min_v: f64, max_v: f64) -> f64 {
    if val < opt_lo {
        let t = (val - min_v) / (opt_lo - min_v).max(1e-6);
        (t * 0.25).clamp(0.0, 0.25)
    } else if val <= opt_hi {
        let t = (val - opt_lo) / (opt_hi - opt_lo).max(1e-6);
        (0.25 + t * 0.60).clamp(0.25, 0.85)
    } else {
        let t = (val - opt_hi) / (max_v - opt_hi).max(1e-6);
        (0.85 + t * 0.15).clamp(0.85, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    struct RefFixture {
        molecules: Vec<RefMol>,
    }

    #[derive(Debug, Deserialize)]
    struct RefMol {
        name: String,
        smiles: String,
        formula: String,
        mw: f64,
        heavy_atoms: usize,
        aromatic_heavy_atoms: usize,
        fsp3: f64,
        rotatable_bonds: usize,
        hbd: usize,
        tpsa: f64,
        wlogp_range: [f64; 2],
        esol_range: [f64; 2],
        gi_absorption: String,
        bbb_permeant: Option<bool>,
        lipinski_yes: bool,
        lipinski_max_violations: usize,
        veber_yes: bool,
    }

    #[test]
    fn reference_15_drugs_match_swissadme_corpus() -> Result<(), Box<dyn std::error::Error>> {
        let raw = std::fs::read_to_string("../../tests/fixtures/reference_molecules.json")?;
        let fixture: RefFixture = serde_json::from_str(&raw)?;
        assert_eq!(fixture.molecules.len(), 15);

        let mut failures = Vec::new();
        for r in &fixture.molecules {
            let eval = evaluate_smiles(&r.smiles, Some(&r.name))?;
            eprintln!(
                "{:<14} {:<12} MW={:>6.2} HA={:>2} Ar={:>2} Fsp3={:.2} RB={:>2} HBD={} TPSA={:>6.2} WLOGP={:>5.2} ESOL={:>5.2} GI={} BBB={} Lip={}({}) Veb={}",
                r.name,
                eval.physicochemical.formula,
                eval.physicochemical.molecular_weight,
                eval.physicochemical.heavy_atoms,
                eval.physicochemical.aromatic_heavy_atoms,
                eval.physicochemical.fraction_csp3,
                eval.physicochemical.rotatable_bonds,
                eval.physicochemical.h_bond_donors,
                eval.physicochemical.tpsa,
                eval.lipophilicity.wlogp,
                eval.water_solubility.esol.log_s,
                eval.pharmacokinetics.gi_absorption.as_str(),
                eval.pharmacokinetics.bbb_permeant,
                eval.druglikeness.lipinski.pass,
                eval.druglikeness.lipinski.violations,
                eval.druglikeness.veber.pass,
            );
            if eval.physicochemical.formula != r.formula {
                failures.push(format!(
                    "{}: formula {} != {}",
                    r.name, eval.physicochemical.formula, r.formula
                ));
            }
            if (eval.physicochemical.molecular_weight - r.mw).abs() >= 0.1 {
                failures.push(format!(
                    "{}: MW {} != {}",
                    r.name, eval.physicochemical.molecular_weight, r.mw
                ));
            }
            if eval.physicochemical.heavy_atoms != r.heavy_atoms {
                failures.push(format!(
                    "{}: HA {} != {}",
                    r.name, eval.physicochemical.heavy_atoms, r.heavy_atoms
                ));
            }
            if eval.physicochemical.aromatic_heavy_atoms != r.aromatic_heavy_atoms {
                failures.push(format!(
                    "{}: ArHA {} != {}",
                    r.name, eval.physicochemical.aromatic_heavy_atoms, r.aromatic_heavy_atoms
                ));
            }
            if (eval.physicochemical.fraction_csp3 - r.fsp3).abs() > 0.02 {
                failures.push(format!(
                    "{}: Fsp3 {} != {}",
                    r.name, eval.physicochemical.fraction_csp3, r.fsp3
                ));
            }
            if eval.physicochemical.rotatable_bonds != r.rotatable_bonds {
                failures.push(format!(
                    "{}: RB {} != {}",
                    r.name, eval.physicochemical.rotatable_bonds, r.rotatable_bonds
                ));
            }
            if eval.physicochemical.h_bond_donors != r.hbd {
                failures.push(format!(
                    "{}: HBD {} != {}",
                    r.name, eval.physicochemical.h_bond_donors, r.hbd
                ));
            }
            if (eval.physicochemical.tpsa - r.tpsa).abs() > 0.25 {
                failures.push(format!(
                    "{}: TPSA {} != {}",
                    r.name, eval.physicochemical.tpsa, r.tpsa
                ));
            }
            if eval.lipophilicity.wlogp < r.wlogp_range[0]
                || eval.lipophilicity.wlogp > r.wlogp_range[1]
            {
                failures.push(format!(
                    "{}: WLOGP {} not in {:?}",
                    r.name, eval.lipophilicity.wlogp, r.wlogp_range
                ));
            }
            if eval.water_solubility.esol.log_s < r.esol_range[0]
                || eval.water_solubility.esol.log_s > r.esol_range[1]
            {
                failures.push(format!(
                    "{}: ESOL {} not in {:?}",
                    r.name, eval.water_solubility.esol.log_s, r.esol_range
                ));
            }
            if eval.pharmacokinetics.gi_absorption.as_str() != r.gi_absorption {
                failures.push(format!(
                    "{}: GI {} != {}",
                    r.name,
                    eval.pharmacokinetics.gi_absorption.as_str(),
                    r.gi_absorption
                ));
            }
            if let Some(bbb) = r.bbb_permeant
                && eval.pharmacokinetics.bbb_permeant != bbb
            {
                failures.push(format!(
                    "{}: BBB {} != {}",
                    r.name, eval.pharmacokinetics.bbb_permeant, bbb
                ));
            }
            if eval.druglikeness.lipinski.pass != r.lipinski_yes {
                failures.push(format!(
                    "{}: Lipinski {} != {}",
                    r.name, eval.druglikeness.lipinski.pass, r.lipinski_yes
                ));
            }
            if eval.druglikeness.lipinski.violations > r.lipinski_max_violations {
                failures.push(format!(
                    "{}: Lipinski viol {} > {}",
                    r.name, eval.druglikeness.lipinski.violations, r.lipinski_max_violations
                ));
            }
            if eval.druglikeness.veber.pass != r.veber_yes {
                failures.push(format!(
                    "{}: Veber {} != {}",
                    r.name, eval.druglikeness.veber.pass, r.veber_yes
                ));
            }
        }
        assert!(
            failures.is_empty(),
            "Reference mismatches:\n{}",
            failures.join("\n")
        );
        Ok(())
    }

    #[test]
    fn pains_detects_rhodanine_and_quinone() -> Result<(), Box<dyn std::error::Error>> {
        let quinone = evaluate_smiles("O=C1C=CC(=O)C=C1", Some("Benzoquinone"))?;
        assert!(!quinone.medicinal_chemistry.pains_alerts.is_empty());
        assert!(!quinone.medicinal_chemistry.brenk_alerts.is_empty());

        let aspirin = evaluate_smiles("CC(=O)Oc1ccccc1C(=O)O", Some("Aspirin"))?;
        assert!(aspirin.medicinal_chemistry.pains_alerts.is_empty());
        assert_eq!(
            aspirin.medicinal_chemistry.pains_catalog_size,
            medchem::PAINS_CATALOG_SIZE
        );
        assert_eq!(
            aspirin.medicinal_chemistry.brenk_catalog_size,
            medchem::BRENK_CATALOG_SIZE
        );
        assert!(aspirin.medicinal_chemistry.sa_is_heuristic_proxy);

        // Verify QED, CNS MPO, Ionization/LogD7.4, Industrial rules & Safety panel
        let ibu = evaluate_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1", Some("Ibuprofen"))?;
        assert!(
            (ibu.druglikeness.qed - 0.822).abs() < 0.06,
            "ibu QED={}",
            ibu.druglikeness.qed
        );
        assert_eq!(ibu.ionization.ionization_class, IonizationClass::Acid);
        assert!(ibu.ionization.logd_ph74 < ibu.lipophilicity.consensus_log_p - 1.5);
        assert_eq!(ibu.safety.herg_risk, HergRisk::Low);

        let diaz = evaluate_smiles("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", Some("Diazepam"))?;
        assert!(
            diaz.druglikeness.cns_mpo >= 4.0,
            "diazepam CNS MPO={}",
            diaz.druglikeness.cns_mpo
        );
        assert!(diaz.medicinal_chemistry.brenk_alerts.is_empty());

        // Verify [Br] SMARTS parsing fix (Benzyl bromide -> alkyl_halide & AMES alert)
        let bbr = evaluate_smiles("BrCc1ccccc1", Some("Benzyl_bromide"))?;
        assert!(
            bbr.medicinal_chemistry
                .brenk_alerts
                .iter()
                .any(|a| a == "alkyl_halide"),
            "Benzyl bromide missed alkyl_halide: {:?}",
            bbr.medicinal_chemistry.brenk_alerts
        );
        assert!(
            bbr.safety
                .ames_alerts
                .iter()
                .any(|a| a == "alkyl_halide_mutagen"),
            "Benzyl bromide missed AMES alert: {:?}",
            bbr.safety.ames_alerts
        );

        let cmp = compare_smiles(
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
            Some("Ibuprofen"),
            "CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1",
            Some("Ibu_Tetrazole"),
        )?;
        assert!(cmp.delta_tpsa > 15.0);
        Ok(())
    }

    #[test]
    fn sa_score_distinguishes_ethanol_from_taxol() -> Result<(), Box<dyn std::error::Error>> {
        let ethanol = evaluate_smiles("CCO", Some("Ethanol"))?;
        assert!(
            ethanol.medicinal_chemistry.synthetic_accessibility <= 1.2,
            "Ethanol SA too high: {}",
            ethanol.medicinal_chemistry.synthetic_accessibility
        );
        let taxol = evaluate_smiles(
            "CC(=O)OC1C(=O)C2(C)C(O)CC3OCC3(OC(C)=O)C2C(OC(=O)c2ccccc2)C2(O)CC(OC(=O)C(O)C(NC(=O)c3ccccc3)c3ccccc3)C(C)=C1C2(C)C",
            Some("Taxol"),
        )?;
        assert!(
            taxol.medicinal_chemistry.synthetic_accessibility > 5.0
                && taxol.medicinal_chemistry.synthetic_accessibility < 9.8,
            "Taxol SA unexpected: {}",
            taxol.medicinal_chemistry.synthetic_accessibility
        );
        Ok(())
    }

    #[test]
    fn experimental_cyp_and_pgp_heuristics_run_cleanly() -> Result<(), Box<dyn std::error::Error>> {
        let atorvastatin = evaluate_smiles(
            "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
            Some("Atorvastatin"),
        )?;
        assert!(
            atorvastatin
                .pharmacokinetics
                .pgp_cyp_are_experimental_heuristics
        );
        assert!(atorvastatin.pharmacokinetics.pgp_substrate);
        assert!(atorvastatin.pharmacokinetics.cyp3a4_inhibitor);

        let fluoxetine = evaluate_smiles("CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", Some("Fluoxetine"))?;
        assert!(fluoxetine.pharmacokinetics.cyp2d6_inhibitor);
        Ok(())
    }

    #[test]
    fn warnings_propagate_for_salts_and_stereo() -> Result<(), Box<dyn std::error::Error>> {
        let salt = evaluate_smiles("CC(=O)[O-].[Na+]", Some("SodiumAcetate"))?;
        assert!(salt.warnings.iter().any(|w| w.contains("Desalted")));
        assert!((salt.physicochemical.molecular_weight - 59.04).abs() < 0.1);

        let chiral = evaluate_smiles("C[C@@H](C(=O)O)N", Some("L-Alanine"))?;
        assert!(
            chiral
                .warnings
                .iter()
                .any(|w| w.contains("Stereochemistry"))
        );
        Ok(())
    }
}
