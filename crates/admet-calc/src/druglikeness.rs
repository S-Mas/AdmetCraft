use chem_core::Molecule;

use crate::pharmacokinetics::has_carboxylic_or_sulfonic_acid;
use crate::physicochem::{count_n_or_o, count_nh_or_oh};
use crate::types::{Druglikeness, Lipophilicity, PhysicochemicalProperties, RuleResult};

/// Evaluate Lipinski, Ghose, Veber, Egan, Muegge, and Abbott Bioavailability Score.
pub fn compute_druglikeness(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> Druglikeness {
    let lipinski = eval_lipinski(mol, phys, lipo);
    let ghose = eval_ghose(mol, phys, lipo);
    let veber = eval_veber(phys);
    let egan = eval_egan(phys, lipo);
    let muegge = eval_muegge(mol, phys, lipo);
    let bioavailability_score = compute_abbott_bioavailability(mol, phys, &lipinski);

    Druglikeness {
        lipinski,
        ghose,
        veber,
        egan,
        muegge,
        bioavailability_score,
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
