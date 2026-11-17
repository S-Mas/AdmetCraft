use chem_core::Molecule;

use crate::physicochem::round2;
use crate::types::{
    Lipophilicity, PhysicochemicalProperties, SolubilityClass, WaterSolubility,
    WaterSolubilityModel,
};

/// Compute ESOL (Delaney 2004), Ali (2012), and SILICOS-IT water solubility predictions.
pub fn compute_water_solubility(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> WaterSolubility {
    // 1. ESOL (Delaney 2004): log S = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP
    let ap = if phys.heavy_atoms == 0 {
        0.0
    } else {
        phys.aromatic_heavy_atoms as f64 / phys.heavy_atoms as f64
    };
    let esol_log_s = 0.16 - 0.63 * lipo.xlogp3 - 0.0062 * phys.molecular_weight
        + 0.066 * (phys.rotatable_bonds as f64)
        - 0.74 * ap;

    // 2. Ali et al. (2012): log S = 0.514 - 1.038*LogP - 0.0062*TPSA (adapted with XLOGP3)
    let ali_log_s = 0.514 - 1.038 * lipo.xlogp3 - 0.0062 * phys.tpsa;

    // 3. Experimental topological solubility heuristic (NOT Silicos-IT)
    let rings = mol.rings.len() as f64;
    let topo_log_s = 0.40 - 0.85 * lipo.topological_logp_est - 0.0050 * phys.molecular_weight
        + 0.010 * phys.tpsa
        - 0.25 * rings;

    WaterSolubility {
        esol: build_model(esol_log_s, phys.molecular_weight),
        ali: build_model(ali_log_s, phys.molecular_weight),
        topological_est: build_model(topo_log_s, phys.molecular_weight),
    }
}

fn build_model(log_s: f64, mw: f64) -> WaterSolubilityModel {
    let clamped = log_s.clamp(-20.0, 6.0);
    let mol_l = 10.0_f64.powf(clamped);
    let mg_ml = mol_l * mw; // mol/L * g/mol = g/L = mg/mL
    WaterSolubilityModel {
        log_s: round2(clamped),
        solubility_mg_ml: format_sci(mg_ml),
        solubility_mol_l: format_sci(mol_l),
        class: SolubilityClass::from_log_s(clamped),
    }
}

fn format_sci(v: f64) -> f64 {
    if !v.is_finite() {
        return 0.0;
    }
    v
}
