use std::collections::BTreeSet;
use std::sync::OnceLock;

use chem_core::{Molecule, SmartsPattern};

use crate::physicochem::round2;
use crate::types::{Lipophilicity, MedicinalChemistry, PhysicochemicalProperties, RuleResult};

struct CompiledRule {
    label: &'static str,
    pattern: SmartsPattern,
}

static PAINS_RULES: OnceLock<Vec<CompiledRule>> = OnceLock::new();
static BRENK_RULES: OnceLock<Vec<CompiledRule>> = OnceLock::new();

fn compile_rules(defs: &[(&'static str, &'static str)]) -> Vec<CompiledRule> {
    defs.iter()
        .filter_map(|&(label, smarts)| {
            SmartsPattern::parse(smarts)
                .ok()
                .map(|pattern| CompiledRule { label, pattern })
        })
        .collect()
}

/// Curated PAINS A, B, and C SMARTS alerts (Baell & Holloway 2010).
fn pains_rules() -> &'static [CompiledRule] {
    PAINS_RULES.get_or_init(|| {
        compile_rules(&[
            // PAINS Family A
            ("rhodanine_A", "[#6]1(=[O,S])[#7][#6](=[O,S])[#16][#6]1=[#6]"),
            ("quinone_A", "O=C1C=CC(=O)C=C1"),
            ("quinone_12_A", "O=C1C(=O)C=CC=C1"),
            ("catechol_A", "c1ccc(O)c(O)c1"),
            ("hydroquinone_A", "c1cc(O)ccc1O"),
            ("azo_A", "[#6][#7]=[#7][#6]"),
            ("ene_rhodanine_A", "S=C1NC(=O)C(=C)S1"),
            ("curcumin_ene_one_A", "O=C([#6]=[#6])[#6]=[#6]"),
            ("mannich_A", "c1ccc(O)c(CN([#6])[#6])c1"),
            ("anil_di_alk_A", "c1cc(N([#6])[#6])ccc1N=C"),
            ("hydroxyphenyl_hydrazone_A", "c1ccc(O)c(C=NN)c1"),
            ("toxoflavin_A", "O=c1[nH]c(=O)n2ncnc2n1"),
            // PAINS Family B
            ("ene_five_het_B", "[#6]1(=[#6])[#6](=O)[#7][#6](=[O,S])[#16,#7,#8]1"),
            ("imine_one_isatin_B", "O=C1Nc2ccccc2C1=N"),
            ("thioamide_ene_B", "[#6]=[#6][#6](=S)[#7]"),
            ("het_pyridinium_B", "c1cc[n+]([#6])cc1"),
            ("amino_acridine_B", "c1ccc2nc3ccccc3c(N)c2c1"),
            ("indole_3yl_alk_B", "c1ccc2[nH]cc(C=C)c2c1"),
            ("pyrrole_A_B", "c1cc([#6]=O)[nH]c1"),
            ("cyanopyridone_B", "O=c1[nH]ccc(C#N)c1"),
            // PAINS Family C
            ("quinone_methide_C", "O=C1C=CC(=C)C=C1"),
            ("aminothiazole_C", "c1sc(N)nc1"),
            ("alkylidene_barbiturate_C", "O=C1NC(=O)NC(=O)C1=C"),
            ("hydantoin_ene_C", "O=C1NC(=O)NC1=C"),
            ("maleimide_C", "O=C1C=CC(=O)N1"),
            ("thiosemicarbazone_C", "[#7][#7]=[#6][#16]"),
        ])
    })
}

/// Curated Brenk structural alerts (Brenk et al. 2008).
fn brenk_rules() -> &'static [CompiledRule] {
    BRENK_RULES.get_or_init(|| {
        compile_rules(&[
            ("aldehyde", "[CX3H1](=O)[#6]"),
            ("alkyl_halide", "[CX4][Cl,Br,I]"),
            ("acyl_halide", "[CX3](=O)[F,Cl,Br,I]"),
            ("sulfonyl_halide", "[SX4](=O)(=O)[F,Cl,Br,I]"),
            ("michael_acceptor", "[CX3]=[CX3][CX3](=O)[#6,#8,#7]"),
            ("epoxide", "C1OC1"),
            ("aziridine", "C1NC1"),
            ("peroxide", "[OX2][OX2]"),
            ("disulfide", "[SX2][SX2]"),
            ("thiol", "[SX2H1]"),
            ("hydrazine", "[NX3][NX3]"),
            ("hydroxamic_acid", "[CX3](=O)[NX3][OX2H1]"),
            ("nitro_group", "[NX3](=O)=O"),
            ("azide", "[NX2]=[NX2]=[NX1]"),
            ("isocyanate", "[NX2]=[CX2]=[OX1]"),
            ("isothiocyanate", "[NX2]=[CX2]=[SX1]"),
            ("thiourea", "[NX3][CX3](=S)[NX3]"),
            ("beta_keto_ester_or_anhydride", "O=C([#6])OC(=O)[#6]"),
            ("polyene_4", "[CX3]=[CX3][CX3]=[CX3][CX3]=[CX3][CX3]=[CX3]"),
            ("triple_bond_terminal", "[CX2H1]#[CX2]"),
            ("phosphor_ester", "[PX4](=O)([OX2])[OX2]"),
            ("thioester", "[CX3](=O)[SX2][#6]"),
            ("imine", "[CX3]=[NX2]"),
            ("alpha_halo_ketone", "[CX3](=O)[CX4][Cl,Br,I]"),
            ("polycyclic_aromatic_3_fused", "c1ccc2c(c1)ccc3ccccc23"),
        ])
    })
}

/// Compute PAINS alerts, Brenk alerts, Leadlikeness, and Synthetic Accessibility (1..=10).
pub fn compute_medicinal_chemistry(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> MedicinalChemistry {
    let pains_alerts: Vec<String> = pains_rules()
        .iter()
        .filter(|r| r.pattern.matches(mol))
        .map(|r| r.label.to_string())
        .collect();

    let brenk_alerts: Vec<String> = brenk_rules()
        .iter()
        .filter(|r| r.pattern.matches(mol))
        .map(|r| r.label.to_string())
        .collect();

    let leadlikeness = eval_leadlikeness(phys, lipo);
    let synthetic_accessibility = compute_sa_score(mol);

    MedicinalChemistry {
        pains_alerts,
        brenk_alerts,
        leadlikeness,
        synthetic_accessibility,
    }
}

fn eval_leadlikeness(phys: &PhysicochemicalProperties, lipo: &Lipophilicity) -> RuleResult {
    let mut details = Vec::new();
    if phys.molecular_weight < 250.0 {
        details.push(format!("MW<250 ({:.2})", phys.molecular_weight));
    } else if phys.molecular_weight > 350.0 {
        details.push(format!("MW>350 ({:.2})", phys.molecular_weight));
    }
    if lipo.xlogp3 > 3.5 {
        details.push(format!("XLOGP3>3.5 ({:.2})", lipo.xlogp3));
    }
    if phys.rotatable_bonds > 7 {
        details.push(format!("Rotors>7 ({})", phys.rotatable_bonds));
    }
    let violations = details.len();
    RuleResult {
        pass: violations == 0,
        violations,
        details,
    }
}

/// Ertl & Schuffenhauer (2009) *J. Cheminf.* 1:8 Synthetic Accessibility Score (1.0 easy .. 10.0 hard).
pub fn compute_sa_score(mol: &Molecule) -> f64 {
    let heavy = mol.heavy_atom_count().max(1);
    let n_atoms = heavy as f64;

    // 1. Fragment commonality proxy from circular 1-bond environments (ECFP2-like)
    let mut env_signatures: BTreeSet<(u8, usize, bool, u8)> = BTreeSet::new();
    let mut uncommon_count = 0usize;
    for (idx, a) in mol.atoms.iter().enumerate() {
        if a.atomic_number == 1 {
            continue;
        }
        let deg = mol.heavy_degree(idx);
        env_signatures.insert((a.atomic_number, deg, a.aromatic, a.total_h()));
        if !matches!(a.atomic_number, 6 | 7 | 8 | 9 | 16 | 17) || deg >= 4 && a.in_ring {
            uncommon_count += 1;
        }
    }
    let diversity_ratio = env_signatures.len() as f64 / n_atoms;

    // 2. Topological complexity penalties (Ertl & Schuffenhauer 2009)
    let size_penalty = if n_atoms > 12.0 {
        (n_atoms.powf(1.005) - n_atoms) + 0.04 * (n_atoms - 12.0)
    } else {
        0.0
    };

    // Chiral centers (explicit tetrahedral tags + stereogenic sp3 C with >=3 distinct heavy nbrs)
    let stereo_centers = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|&(idx, a)| {
            a.chiral != chem_core::ChiralTag::None
                || (a.atomic_number == 6
                    && !a.aromatic
                    && mol.heavy_degree(idx) >= 3
                    && a.in_ring
                    && mol.rings.iter().filter(|r| r.contains(&idx)).count() >= 2)
        })
        .count() as f64;
    let stereo_penalty = (stereo_centers + 1.0).ln() * 0.65;

    // Ring bridgehead & spiro complexity (atoms shared by >= 2 rings)
    let mut ring_shared_atoms = 0usize;
    let mut spiro_atoms = 0usize;
    for idx in 0..mol.atoms.len() {
        let r_cnt = mol.rings.iter().filter(|r| r.contains(&idx)).count();
        if r_cnt >= 2 {
            ring_shared_atoms += 1;
            if mol.heavy_degree(idx) >= 4 {
                spiro_atoms += 1;
            }
        }
    }
    let bridge_penalty = if ring_shared_atoms >= 3 {
        0.55 * ((ring_shared_atoms - 1) as f64).ln()
    } else {
        0.15 * (ring_shared_atoms as f64)
    };
    let spiro_penalty = 0.50 * (spiro_atoms as f64);

    // Macrocycle penalty (rings > 8 atoms)
    let macro_penalty = mol
        .rings
        .iter()
        .filter(|r| r.len() > 8)
        .count() as f64
        * 0.8;

    // Small strained ring penalty (3- or 4-membered rings like cyclopropyl, beta-lactam)
    let strained_penalty = mol
        .rings
        .iter()
        .filter(|r| r.len() <= 4)
        .count() as f64
        * 0.55;

    let raw = 1.05
        + 0.95 * diversity_ratio
        + 0.25 * (uncommon_count as f64)
        + size_penalty
        + stereo_penalty
        + bridge_penalty
        + spiro_penalty
        + macro_penalty
        + strained_penalty;

    round2(raw.clamp(1.0, 10.0))
}
