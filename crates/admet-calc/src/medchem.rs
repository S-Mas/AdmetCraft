use std::collections::BTreeSet;
use std::sync::OnceLock;

use chem_core::{Molecule, SmartsPattern};

use crate::physicochem::round2;
use crate::types::{
    HergRisk, IonizationProfile, Lipophilicity, MedicinalChemistry, PhysicochemicalProperties,
    RuleResult, SafetyProfile,
};

struct CompiledRule {
    label: &'static str,
    pattern: SmartsPattern,
}

static PAINS_RULES: OnceLock<Vec<CompiledRule>> = OnceLock::new();
static BRENK_RULES: OnceLock<Vec<CompiledRule>> = OnceLock::new();
static AMES_RULES: OnceLock<Vec<CompiledRule>> = OnceLock::new();
static DILI_RULES: OnceLock<Vec<CompiledRule>> = OnceLock::new();

fn compile_rules(defs: &[(&'static str, &'static str)]) -> Vec<CompiledRule> {
    defs.iter()
        .filter_map(|&(label, smarts)| {
            SmartsPattern::parse(smarts)
                .ok()
                .map(|pattern| CompiledRule { label, pattern })
        })
        .collect()
}

/// Curated PAINS A, B, and C SMARTS alerts (Baell & Holloway 2010) calibrated against RDKit FilterCatalog
/// to avoid false positives on approved clinical drugs (e.g. free catechol `[OX2H1]` vs aryl ether).
fn pains_rules() -> &'static [CompiledRule] {
    PAINS_RULES.get_or_init(|| {
        compile_rules(&[
            // PAINS Family A
            (
                "rhodanine_A",
                "[#6]1(=[O,S])[#7][#6](=[O,S])[#16][#6]1=[#6]",
            ),
            ("quinone_A", "O=C1C=CC(=O)C=C1"),
            ("quinone_12_A", "O=C1C(=O)C=CC=C1"),
            ("catechol_A", "c1ccc([OX2H1])c([OX2H1])c1"),
            ("anil_no_alk_A", "[NX3H2]c1ccc([N+](=O)[O-])cc1"),
            ("azo_A", "[#6][#7]=[#7][#6]"),
            ("azo_azide_A", "[#7]=[#7+]=[#7-]"),
            ("ene_rhodanine_A", "S=C1NC(=O)C(=C)S1"),
            (
                "curcumin_ene_one_A",
                "[#6;!R]=[#6;!R][#6](=O)[#6;!R]=[#6;!R]",
            ),
            ("mannich_A", "c1ccc([OX2H1])c(CN([#6])[#6])c1"),
            ("anil_di_alk_A", "c1cc(N([#6])[#6])ccc1[NX2]=[#6]"),
            (
                "hydroxyphenyl_hydrazone_A",
                "c1ccc([OX2H1])c([CX3]=[NX2][NX3])c1",
            ),
            ("toxoflavin_A", "O=c1[nH]c(=O)n2ncnc2n1"),
            (
                "anil_di_alk_B_para",
                "c1cc([NX3]([#6])[#6])ccc1[NX3]([#6])[#6]",
            ),
            ("indole_3yl_alk_A", "c1ccc2[nH]cc([CX3;!R]=[CX3;!R])c2c1"),
            ("pyrrole_A", "c1cc([CX3](=O)[#6,#8])[nH]c1[CX3](=O)"),
            ("thiosemicarbazone_A", "[NX3][NX2]=[CX3][SX2,SX1]"),
            // PAINS Family B
            (
                "ene_five_het_B",
                "[#6]1(=[#6;!R])[#6](=O)[#7][#6](=[O,S])[#16,#7,#8]1",
            ),
            ("thioamide_ene_B", "[#6;!R]=[#6;!R][#6](=S)[#7]"),
            ("het_pyridinium_B", "c1cc[n+]([#6])cc1"),
            ("amino_acridine_B", "c1ccc2nc3ccccc3c([NX3H2])c2c1"),
            ("cyanopyridone_B", "O=c1[nH]ccc(C#N)c1"),
            ("benzofurazan_B", "c1ccc2nonc2c1"),
            ("ene_cyano_B", "[CX3;!R]=[CX3;!R](C#N)C#N"),
            ("meldrum_ene_B", "O=C1OC(C)(C)OC(=O)C1=C"),
            ("naphthoquinone_B", "O=C1C=CC(=O)c2ccccc12"),
            ("anthraquinone_B", "O=C1c2ccccc2C(=O)c2ccccc12"),
            ("phenothiazine_B", "c1ccc2c(c1)[nH]c1ccccc1s2"),
            ("aminopyrimidine_nitroso_B", "O=Nc1c(N)ncnc1N"),
            ("nitro_furan_B", "c1cc([N+](=O)[O-])oc1"),
            // PAINS Family C
            ("quinone_methide_C", "O=C1C=CC(=C)C=C1"),
            ("aminothiazole_ene_C", "c1sc(N)nc1[CX3]=[CX3]"),
            ("alkylidene_barbiturate_C", "O=C1NC(=O)NC(=O)C1=C"),
            ("hydantoin_ene_C", "O=C1NC(=O)NC1=C"),
            ("ene_six_het_A", "O=C1C(=[#6;!R])C(=O)NC(=O)N1"),
            ("isothiazolone_C", "O=c1sncc1"),
            ("ene_nitro_C", "[CX3;!R]=[CX3;!R][N+](=O)[O-]"),
            ("sulfonamide_azo_C", "S(=O)(=O)c1ccc(N=N)cc1"),
            ("acyl_hydrazone_C", "[CX3](=O)[NX3][NX2]=[CX3]"),
            ("enaminone_C", "[NX3][CX3;!R]=[CX3;!R][CX3](=O)"),
        ])
    })
}

/// Expanded Brenk structural alerts (Brenk et al. 2008, *ChemMedChem* 3, 435-444)
/// calibrated against RDKit's BRENK FilterCatalog (55 high-precision SMARTS).
fn brenk_rules() -> &'static [CompiledRule] {
    BRENK_RULES.get_or_init(|| {
        compile_rules(&[
            ("aldehyde", "[CX3H1](=O)[#6]"),
            ("alkyl_halide", "[CX4][Cl,Br,I]"),
            ("acyl_halide", "[CX3](=O)[F,Cl,Br,I]"),
            ("sulfonyl_halide", "[SX4](=O)(=O)[F,Cl,Br,I]"),
            ("michael_acceptor", "[CX3;!R]=[CX3][CX3](=O)[#6,#8,#7]"),
            ("epoxide", "C1OC1"),
            ("aziridine", "C1NC1"),
            ("peroxide", "[OX2][OX2]"),
            ("disulfide", "[SX2][SX2]"),
            ("thiol", "[SX2H1]"),
            ("hydrazine", "[NX3][NX3]"),
            ("hydroxamic_acid", "[CX3](=O)[NX3][OX2H1]"),
            ("nitro_group", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
            ("azide", "[$([NX2]=[NX2]=[NX1]),$([NX2]=[NX2+]=[NX1-])]"),
            ("isocyanate", "[NX2]=[CX2]=[OX1]"),
            ("isothiocyanate", "[NX2]=[CX2]=[SX1]"),
            ("thiourea", "[NX3][CX3](=S)[NX3]"),
            ("anhydride", "O=C([#6])OC(=O)[#6]"),
            ("polyene_4", "[CX3]=[CX3][CX3]=[CX3][CX3]=[CX3][CX3]=[CX3]"),
            ("triple_bond_terminal", "[CX2H1]#[CX2]"),
            ("phosphor_ester", "[PX4](=O)([OX2])[OX2]"),
            ("thioester", "[CX3](=O)[SX2][#6]"),
            ("acyclic_imine", "[CX3]=[NX2;!R]"),
            ("alpha_halo_ketone", "[CX3](=O)[CX4][Cl,Br,I]"),
            ("polycyclic_aromatic_3_fused", "c1ccc2c(c1)ccc3ccccc23"),
            // High-value RDKit Brenk parity additions
            ("phenol_ester", "[OX2]([c])[CX3](=O)[#6]"),
            (
                "hydroquinone_or_aminophenol",
                "[OX2H1]c1ccc([OX2H1,NX3])cc1",
            ),
            ("catechol_brenk", "[OX2H1]c1ccccc1[OX2H1]"),
            ("aniline_free", "[NX3H2]c1ccccc1"),
            (
                "isolated_alkene",
                "[CX3;H1,H2;!$(C=CC=O)]=[CX3;H1,H2;!$(C=CC=O)]",
            ),
            ("coumarin", "O=c1ccc2ccccc2o1"),
            ("sulfonate_ester", "[SX4](=O)(=O)[OX2][#6]"),
            ("sulfate_ester", "[OX2][SX4](=O)(=O)[OX2]"),
            ("oxime", "[CX3]=[NX2][OX2H1]"),
            ("nitroso", "[NX2]=O"),
            ("diazo", "[#6]=[NX2+]=[NX1-]"),
            ("thiocarbonyl", "[CX3]=S"),
            ("boronic_acid", "[BX3]([OX2])[OX2]"),
            ("quat_ammonium", "[NX4+]([#6])([#6])([#6])[#6]"),
            ("azo_brenk", "[#7]=[#7]"),
            ("acyl_hydrazine", "[CX3](=O)[NX3][NX3]"),
            ("halo_alkene", "[CX3]=[CX3][Cl,Br,I]"),
            ("gem_dihalide", "[CX4]([Cl,Br,I])[Cl,Br,I]"),
            (
                "polyhalogenated_arom",
                "c1c([Cl,Br,I])c([Cl,Br,I])c([Cl,Br,I])cc1",
            ),
            ("enamine", "[CX3;!R]=[CX3;!R][NX3]"),
            ("lawesson_ps", "[PX4]=S"),
            ("cyanohydrin", "[CX4]([OX2H1])C#N"),
            ("nitrite_nitrate", "[OX2][NX3+](=O)[O-]"),
            ("chloramidine", "[CX3](=[NX2])[Cl,Br]"),
            ("carbodiimide", "[NX2]=[CX2]=[NX2]"),
            ("ketene", "[CX3]=[CX2]=O"),
            ("vinyl_sulfone", "[CX3;!R]=[CX3;!R][SX4](=O)(=O)"),
            ("maleimide_brenk", "O=C1C=CC(=O)N1"),
            ("quinone_brenk", "O=C1C=CC(=O)C=C1"),
            ("stilbene", "c1ccccc1[CX3]=[CX3]c2ccccc2"),
        ])
    })
}

/// Benigni-Bossa / Kazius structural alert rules for AMES mutagenicity screening (heuristic).
fn ames_rules() -> &'static [CompiledRule] {
    AMES_RULES.get_or_init(|| {
        compile_rules(&[
            ("aromatic_nitro", "c[N+](=O)[O-]"),
            ("aromatic_amine_primary", "[NX3H2]c1ccccc1"),
            ("aromatic_hydroxylamine", "[NX3]([OX2H1])c1ccccc1"),
            ("nitroso_mutagen", "[NX2]=O"),
            ("epoxide_mutagen", "C1OC1"),
            ("aziridine_mutagen", "C1NC1"),
            ("alkyl_halide_mutagen", "[CX4][Cl,Br,I]"),
            ("hydrazine_mutagen", "[NX3][NX3]"),
            ("azo_aromatic", "c[NX2]=[NX2]c"),
            ("azide_mutagen", "[NX2]=[NX2+]=[NX1-]"),
            (
                "michael_acceptor_mutagen",
                "[CX3;!R]=[CX3;!R][CX3](=O)[#6,#8]",
            ),
            ("sulfonate_ester_alkylator", "[SX4](=O)(=O)[OX2][CX4]"),
            ("pah_3_rings", "c1ccc2c(c1)ccc3ccccc23"),
            ("acyl_halide_mutagen", "[CX3](=O)[Cl,Br]"),
            ("alpha_halo_carbonyl", "[CX3](=O)[CX4][Cl,Br,I]"),
        ])
    })
}

/// Structural bioactivation motifs linked to DILI / reactive metabolite formation (heuristic).
fn dili_rules() -> &'static [CompiledRule] {
    DILI_RULES.get_or_init(|| {
        compile_rules(&[
            ("furan_bioactivation", "c1ccoc1"),
            ("thiophene_bioactivation", "c1ccsc1"),
            ("aniline_acetanilide", "[NX3;H1,H2]c1ccccc1"),
            ("para_aminophenol_napqi", "[OX2H1]c1ccc([NX3])cc1"),
            ("catechol_oquinone", "[OX2H1]c1ccccc1[OX2H1]"),
            ("hydroquinone_pquinone", "[OX2H1]c1ccc([OX2H1])cc1"),
            ("quinone_reactive", "O=C1C=CC(=O)C=C1"),
            ("nitroaromatic_dili", "c[N+](=O)[O-]"),
            ("hydrazine_dili", "[NX3][NX3]"),
            ("thiourea_dili", "[NX3][CX3](=S)[NX3]"),
            ("terminal_alkyne_mbi", "[CX2H1]#[CX2]"),
            ("epoxide_dili", "C1OC1"),
            ("michael_acceptor_dili", "[CX3;!R]=[CX3;!R][CX3](=O)"),
            ("rhodanine_dili", "O=C1NC(=S)SC1"),
        ])
    })
}

/// Number of curated PAINS SMARTS filters implemented.
pub const PAINS_CATALOG_SIZE: usize = 40;
/// Number of curated Brenk SMARTS filters implemented (55 SMARTS + `>_2_ester_groups` count rule).
pub const BRENK_CATALOG_SIZE: usize = 56;

/// Compute PAINS alerts, Brenk alerts, Leadlikeness, and heuristic Synthetic Accessibility (1..=10).
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

    let mut brenk_alerts: Vec<String> = brenk_rules()
        .iter()
        .filter(|r| r.pattern.matches(mol))
        .map(|r| r.label.to_string())
        .collect();

    // RDKit Brenk count rule: >2 ester groups (e.g. Taxol / polyesters)
    static ESTER_PAT: OnceLock<Option<SmartsPattern>> = OnceLock::new();
    if let Some(pat) = ESTER_PAT.get_or_init(|| SmartsPattern::parse("[CX3](=O)[OX2][#6]").ok())
        && pat.count_unique_matches(mol) > 2
    {
        brenk_alerts.push(">_2_ester_groups".to_string());
    }

    let leadlikeness = eval_leadlikeness(phys, lipo);
    let synthetic_accessibility = compute_sa_score(mol);

    MedicinalChemistry {
        pains_alerts,
        pains_catalog_size: pains_rules().len(),
        brenk_alerts,
        brenk_catalog_size: brenk_rules().len() + 1,
        leadlikeness,
        synthetic_accessibility,
        sa_is_heuristic_proxy: true,
    }
}

/// Compute hERG cardiotoxicity risk (descriptor/pharmacophore heuristic) plus AMES & DILI structural alerts.
pub fn compute_safety_profile(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
    ion: &IonizationProfile,
) -> SafetyProfile {
    let mut herg_score = 0u8;
    let mut herg_reasons = Vec::new();

    let is_basic = ion.strongest_basic_pka.is_some_and(|pka| pka >= 7.5);
    let is_acidic = ion.strongest_acidic_pka.is_some_and(|pka| pka <= 5.5);

    if is_basic {
        herg_score += 2;
        herg_reasons.push(format!(
            "Basic ionizable center (pKa≈{:.1})",
            ion.strongest_basic_pka.unwrap_or(8.0)
        ));
    }
    if lipo.consensus_log_p >= 3.5 {
        herg_score += 2;
        herg_reasons.push(format!(
            "High lipophilicity (cLogP={:.2})",
            lipo.consensus_log_p
        ));
    } else if lipo.consensus_log_p >= 2.5 {
        herg_score += 1;
    }
    if phys.num_aromatic_rings >= 2 {
        herg_score += 1;
        herg_reasons.push(format!(
            "Multiple aromatic rings ({}) for pi-stacking",
            phys.num_aromatic_rings
        ));
    }
    if phys.molecular_weight >= 350.0 && phys.tpsa < 85.0 {
        herg_score += 1;
        herg_reasons.push(format!(
            "MW={:.0} with low polar surface (TPSA={:.1})",
            phys.molecular_weight, phys.tpsa
        ));
    }
    if is_acidic && !is_basic {
        // Acidic anions rarely block the hERG K+ channel pore
        herg_score = herg_score.saturating_sub(2);
    }

    let herg_risk = if herg_score >= 5 {
        HergRisk::High
    } else if herg_score >= 3 {
        HergRisk::Moderate
    } else {
        HergRisk::Low
    };

    let ames_alerts: Vec<String> = ames_rules()
        .iter()
        .filter(|r| r.pattern.matches(mol))
        .map(|r| r.label.to_string())
        .collect();

    let dili_alerts: Vec<String> = dili_rules()
        .iter()
        .filter(|r| r.pattern.matches(mol))
        .map(|r| r.label.to_string())
        .collect();

    SafetyProfile {
        herg_risk,
        herg_reasons,
        ames_alerts,
        dili_alerts,
        safety_alerts_are_structural_heuristics: true,
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

/// Heuristic Synthetic Accessibility proxy (1.0 easy .. 10.0 hard).
/// Note: calibrated topological/fragment-diversity proxy; does not embed the 13M PubChem fragment
/// frequency table from Ertl & Schuffenhauer (2009).
pub fn compute_sa_score(mol: &Molecule) -> f64 {
    let heavy = mol.heavy_atom_count().max(1);
    let n_atoms = heavy as f64;

    // 1. Local environment diversity above simple small-molecule baseline
    let mut env_signatures: BTreeSet<(u8, usize, bool, u8)> = BTreeSet::new();
    let mut uncommon_count = 0usize;
    for (idx, a) in mol.atoms.iter().enumerate() {
        if a.atomic_number == 1 {
            continue;
        }
        let deg = mol.heavy_degree(idx);
        env_signatures.insert((a.atomic_number, deg, a.aromatic, a.total_h()));
        if !matches!(a.atomic_number, 6 | 7 | 8 | 9 | 16 | 17) || (deg >= 4 && a.in_ring) {
            uncommon_count += 1;
        }
    }
    let excess_envs = (env_signatures.len() as f64 - 4.0).max(0.0);
    let complexity_ratio = excess_envs / n_atoms;

    // 2. Topological complexity penalties
    let size_penalty = if n_atoms > 12.0 {
        0.042 * (n_atoms - 12.0)
    } else {
        0.0
    };

    // Chiral centers (explicit tetrahedral tags + stereogenic ring junctions)
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
    let stereo_penalty = (stereo_centers + 1.0).ln() * 0.55;

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
    let spiro_penalty = 0.45 * (spiro_atoms as f64);

    // Macrocycle & strained 3/4-membered ring penalties
    let macro_penalty = mol.rings.iter().filter(|r| r.len() > 8).count() as f64 * 0.65;
    let strained_penalty = mol.rings.iter().filter(|r| r.len() <= 4).count() as f64 * 0.45;

    let raw = 1.0
        + 0.85 * complexity_ratio
        + 0.20 * (uncommon_count as f64)
        + size_penalty
        + stereo_penalty
        + bridge_penalty
        + spiro_penalty
        + macro_penalty
        + strained_penalty;

    round2(raw.clamp(1.0, 10.0))
}
