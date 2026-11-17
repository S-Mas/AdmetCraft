use admet_calc::MoleculeEvaluation;

use crate::ReportError;

/// Export evaluations to CSV (roundtrip-compatible with `chem_core::parse_csv_batch`).
pub fn export_csv_string(evals: &[MoleculeEvaluation]) -> Result<String, ReportError> {
    let mut wtr = csv::Writer::from_writer(Vec::new());
    wtr.write_record([
        "Molecule",
        "Canonical SMILES",
        "Formula",
        "MW",
        "#Heavy atoms",
        "#Aromatic heavy atoms",
        "Fraction Csp3",
        "#Rotatable bonds",
        "#H-bond acceptors",
        "#H-bond donors",
        "MR",
        "TPSA",
        "WLOGP",
        "XLOGP3 (heur.)",
        "MLOGP",
        "Solvation LogP (exp.)",
        "Topological LogP (exp.)",
        "Consensus Log P",
        "ESOL Log S",
        "ESOL Solubility (mg/ml)",
        "ESOL Solubility (mol/l)",
        "ESOL Class",
        "Ali Log S",
        "Ali Solubility (mg/ml)",
        "Ali Solubility (mol/l)",
        "Ali Class",
        "Topological LogS (exp.)",
        "Topological Solubility (mg/ml)",
        "Topological Solubility (mol/l)",
        "Topological Class",
        "GI absorption",
        "BBB permeant",
        "Pgp substrate (exp.)",
        "CYP1A2 inhibitor (exp.)",
        "CYP2C19 inhibitor (exp.)",
        "CYP2C9 inhibitor (exp.)",
        "CYP2D6 inhibitor (exp.)",
        "CYP3A4 inhibitor (exp.)",
        "log Kp (cm/s)",
        "Lipinski #violations",
        "Ghose #violations",
        "Veber #violations",
        "Egan #violations",
        "Muegge #violations",
        "Bioavailability Score",
        "PAINS #alerts (26 checked)",
        "PAINS alerts",
        "Brenk #alerts (25 checked)",
        "Brenk alerts",
        "Leadlikeness #violations",
        "Synthetic Accessibility (proxy)",
        "Warnings",
    ])?;

    for e in evals {
        let p = &e.physicochemical;
        let l = &e.lipophilicity;
        let w = &e.water_solubility;
        let pk = &e.pharmacokinetics;
        let d = &e.druglikeness;
        let mc = &e.medicinal_chemistry;
        let pains_str = mc.pains_alerts.join("; ");
        let brenk_str = mc.brenk_alerts.join("; ");
        let warn_str = e.warnings.join("; ");
        wtr.write_record(&[
            sanitize_csv_cell(&e.name),
            e.canonical_smiles.clone(),
            p.formula.clone(),
            format!("{:.2}", p.molecular_weight),
            p.heavy_atoms.to_string(),
            p.aromatic_heavy_atoms.to_string(),
            format!("{:.2}", p.fraction_csp3),
            p.rotatable_bonds.to_string(),
            p.h_bond_acceptors.to_string(),
            p.h_bond_donors.to_string(),
            format!("{:.2}", p.molar_refractivity),
            format!("{:.2}", p.tpsa),
            format!("{:.2}", l.wlogp),
            format!("{:.2}", l.xlogp3),
            format!("{:.2}", l.mlogp),
            format!("{:.2}", l.solvation_logp_est),
            format!("{:.2}", l.topological_logp_est),
            format!("{:.2}", l.consensus_log_p),
            format!("{:.2}", w.esol.log_s),
            format!("{:.2e}", w.esol.solubility_mg_ml),
            format!("{:.2e}", w.esol.solubility_mol_l),
            w.esol.class.as_str().to_string(),
            format!("{:.2}", w.ali.log_s),
            format!("{:.2e}", w.ali.solubility_mg_ml),
            format!("{:.2e}", w.ali.solubility_mol_l),
            w.ali.class.as_str().to_string(),
            format!("{:.2}", w.topological_est.log_s),
            format!("{:.2e}", w.topological_est.solubility_mg_ml),
            format!("{:.2e}", w.topological_est.solubility_mol_l),
            w.topological_est.class.as_str().to_string(),
            pk.gi_absorption.as_str().to_string(),
            yes_no(pk.bbb_permeant).to_string(),
            yes_no(pk.pgp_substrate).to_string(),
            yes_no(pk.cyp1a2_inhibitor).to_string(),
            yes_no(pk.cyp2c19_inhibitor).to_string(),
            yes_no(pk.cyp2c9_inhibitor).to_string(),
            yes_no(pk.cyp2d6_inhibitor).to_string(),
            yes_no(pk.cyp3a4_inhibitor).to_string(),
            format!("{:.2}", pk.log_kp_cm_s),
            d.lipinski.violations.to_string(),
            d.ghose.violations.to_string(),
            d.veber.violations.to_string(),
            d.egan.violations.to_string(),
            d.muegge.violations.to_string(),
            format!("{:.2}", d.bioavailability_score),
            mc.pains_alerts.len().to_string(),
            sanitize_csv_cell(&pains_str),
            mc.brenk_alerts.len().to_string(),
            sanitize_csv_cell(&brenk_str),
            mc.leadlikeness.violations.to_string(),
            format!("{:.2}", mc.synthetic_accessibility),
            sanitize_csv_cell(&warn_str),
        ])?;
    }

    let bytes = wtr
        .into_inner()
        .map_err(|e| ReportError::Format(e.to_string()))?;
    String::from_utf8(bytes).map_err(|e| ReportError::Format(e.to_string()))
}

/// Prevent spreadsheet formula injection (`=`, `+`, `-`, `@`, `\t`, `\r`) in free-text CSV cells.
fn sanitize_csv_cell(val: &str) -> String {
    let trimmed = val.trim_start_matches(' ');
    if trimmed.starts_with(['=', '+', '-', '@', '\t', '\r']) {
        format!("'{val}")
    } else {
        val.to_string()
    }
}

/// Export evaluations to pretty-printed JSON.
pub fn export_json_string(evals: &[MoleculeEvaluation]) -> Result<String, ReportError> {
    Ok(serde_json::to_string_pretty(evals)?)
}

fn yes_no(b: bool) -> &'static str {
    if b { "Yes" } else { "No" }
}
