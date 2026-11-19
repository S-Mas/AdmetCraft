use admet_calc::MoleculeEvaluation;

use crate::ReportError;

/// Export evaluations to CSV (roundtrip-compatible with `chem_core::parse_csv_batch`).
pub fn export_csv_string(evals: &[MoleculeEvaluation]) -> Result<String, ReportError> {
    let mut wtr = csv::Writer::from_writer(Vec::new());
    let pains_n = admet_calc::medchem::PAINS_CATALOG_SIZE;
    let brenk_n = admet_calc::medchem::BRENK_CATALOG_SIZE;
    let pains_hdr = format!("PAINS #alerts ({pains_n} checked)");
    let brenk_hdr = format!("Brenk #alerts ({brenk_n} checked)");
    wtr.write_record([
        "Molecule",
        "Canonical SMILES",
        "Formula",
        "MW",
        "Exact Mass",
        "#Heavy atoms",
        "#Aromatic heavy atoms",
        "Fraction Csp3",
        "#Rotatable bonds",
        "#H-bond acceptors",
        "#H-bond donors",
        "MR",
        "TPSA",
        "TPSA (N/O only)",
        "DBE",
        "WLOGP",
        "XLOGP3 (heur.)",
        "MLOGP",
        "Solvation LogP (exp.)",
        "Topological LogP (exp.)",
        "Consensus Log P",
        "LogD 7.4",
        "Ionization Class (pH 7.4)",
        "Acidic pKa (est.)",
        "Basic pKa (est.)",
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
        "QED (Bickerton 2012)",
        "CNS MPO (Pfizer)",
        "Pfizer 3/75 Pass",
        "GSK 4/400 Pass",
        "Golden Triangle Pass",
        "Ro3 Fragment Pass",
        "bRo5",
        &pains_hdr,
        "PAINS alerts",
        &brenk_hdr,
        "Brenk alerts",
        "Leadlikeness #violations",
        "Synthetic Accessibility (proxy)",
        "hERG Risk (heur.)",
        "AMES Alerts (heur.)",
        "DILI Alerts (heur.)",
        "Warnings",
    ])?;

    for e in evals {
        let p = &e.physicochemical;
        let l = &e.lipophilicity;
        let ion = &e.ionization;
        let w = &e.water_solubility;
        let pk = &e.pharmacokinetics;
        let d = &e.druglikeness;
        let mc = &e.medicinal_chemistry;
        let sf = &e.safety;
        let pains_str = mc.pains_alerts.join("; ");
        let brenk_str = mc.brenk_alerts.join("; ");
        let ames_str = sf.ames_alerts.join("; ");
        let dili_str = sf.dili_alerts.join("; ");
        let warn_str = e.warnings.join("; ");
        wtr.write_record(&[
            sanitize_csv_cell(&e.name),
            e.canonical_smiles.clone(),
            p.formula.clone(),
            format!("{:.2}", p.molecular_weight),
            format!("{:.4}", p.exact_mass),
            p.heavy_atoms.to_string(),
            p.aromatic_heavy_atoms.to_string(),
            format!("{:.2}", p.fraction_csp3),
            p.rotatable_bonds.to_string(),
            p.h_bond_acceptors.to_string(),
            p.h_bond_donors.to_string(),
            format!("{:.2}", p.molar_refractivity),
            format!("{:.2}", p.tpsa),
            format!("{:.2}", p.tpsa_no_only),
            format!("{:.1}", p.dbe),
            format!("{:.2}", l.wlogp),
            format!("{:.2}", l.xlogp3),
            format!("{:.2}", l.mlogp),
            format!("{:.2}", l.solvation_logp_est),
            format!("{:.2}", l.topological_logp_est),
            format!("{:.2}", l.consensus_log_p),
            format!("{:.2}", ion.logd_ph74),
            ion.ionization_class.as_str().to_string(),
            ion.strongest_acidic_pka
                .map_or_else(String::new, |v| format!("{v:.2}")),
            ion.strongest_basic_pka
                .map_or_else(String::new, |v| format!("{v:.2}")),
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
            format!("{:.3}", d.qed),
            format!("{:.2}", d.cns_mpo),
            yes_no(d.pfizer_3_75.pass).to_string(),
            yes_no(d.gsk_4_400.pass).to_string(),
            yes_no(d.golden_triangle.pass).to_string(),
            yes_no(d.ro3_fragment.pass).to_string(),
            yes_no(d.beyond_ro5).to_string(),
            mc.pains_alerts.len().to_string(),
            sanitize_csv_cell(&pains_str),
            mc.brenk_alerts.len().to_string(),
            sanitize_csv_cell(&brenk_str),
            mc.leadlikeness.violations.to_string(),
            format!("{:.2}", mc.synthetic_accessibility),
            sf.herg_risk.as_str().to_string(),
            sanitize_csv_cell(&ames_str),
            sanitize_csv_cell(&dili_str),
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
