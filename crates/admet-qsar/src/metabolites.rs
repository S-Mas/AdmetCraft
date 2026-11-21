#![forbid(unsafe_code)]

//! Pure-Rust **Metabolite Guard**: deterministic Phase I & Phase II hepatic/plasma
//! biotransformation engine operating directly on `chem_core::Molecule` graphs from predicted
//! Sites of Metabolism (`som.rs`), followed by automatic re-screening for **reactive
//! bioactivation** and **toxicity escalation** (DILI / AMES / hERG / Tox21 / structural alerts).

use std::collections::{BTreeSet, HashSet, VecDeque};

use admet_calc::{
    MetabolismPhase, MetabolismSite, MetaboliteGuardReport, MoleculeEvaluation,
    PredictedMetabolite, evaluate_molecule,
};
use chem_core::{Atom, BondOrder, BondStereo, Hybridization, Molecule, parse_smiles};

use crate::fingerprint::extract_mol_features;
use crate::pack::AdmetModelPack;

/// Maximum distinct metabolites re-evaluated per parent molecule by default.
pub const DEFAULT_MAX_METABOLITES: usize = 10;

struct RawMetaboliteCandidate {
    mol: Molecule,
    smiles: String,
    phase: MetabolismPhase,
    reaction: String,
    enzyme: String,
    parent_som_atom: usize,
    lability_score: f64,
    reactive_bioactivation: bool,
    bioactivation_reason: Option<String>,
}

/// Predicts major Phase I & Phase II metabolites from `mol` and re-screens each generated
/// metabolite for reactive bioactivation, newly introduced structural alerts, and QSAR toxicity
/// escalation relative to `parent_eval`.
#[must_use]
pub fn run_metabolite_guard(
    mol: &Molecule,
    parent_eval: &MoleculeEvaluation,
    sites: &[MetabolismSite],
    pack: Option<&AdmetModelPack>,
    max_metabolites: usize,
) -> MetaboliteGuardReport {
    let cap = max_metabolites.clamp(1, 24);
    let mut seen_smiles: HashSet<String> = HashSet::new();
    seen_smiles.insert(parent_eval.canonical_smiles.clone());

    let mut raw_cands: Vec<RawMetaboliteCandidate> = Vec::new();

    for site in sites {
        let idx = site.atom_index;
        if idx >= mol.atoms.len() {
            continue;
        }
        let rxn = site.reaction.as_str();

        if rxn.contains("Benzylic") || rxn.contains("Allylic") || rxn.contains("Tertiary Aliphatic")
        {
            if let Some((smi, met_mol)) = attach_hydroxyl(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: format!("{} -> Alcohol", site.reaction),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: false,
                    bioactivation_reason: None,
                });
            }
        } else if rxn.contains("Aromatic C-H Hydroxylation") {
            if let Some((smi, met_mol)) = attach_hydroxyl(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                let quinone_precursor = has_quinone_or_quinone_imine_precursor(&met_mol);
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "Aromatic Hydroxylation -> Phenol".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: quinone_precursor,
                    bioactivation_reason: if quinone_precursor {
                        Some("Catechol / para-aminophenol quinone(-imine) bioactivation precursor (DILI risk)".to_string())
                    } else {
                        None
                    },
                });
            }
        } else if rxn.contains("O-Dealkylation") {
            if let Some(o_idx) = find_neighbor_by_z(mol, idx, 8)
                && let Some((smi, met_mol)) = cleave_bond_keeping_side(mol, idx, o_idx, o_idx)
                && seen_smiles.insert(smi.clone())
            {
                let quinone_precursor = has_quinone_or_quinone_imine_precursor(&met_mol);
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "O-Dealkylation -> Phenol / Alcohol".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: quinone_precursor,
                    bioactivation_reason: if quinone_precursor {
                        Some(
                            "O-dealkylation unmasks catechol / para-aminophenol quinone precursor"
                                .to_string(),
                        )
                    } else {
                        None
                    },
                });
            }
        } else if rxn.contains("N-Dealkylation") {
            if let Some(n_idx) = find_neighbor_by_z(mol, idx, 7)
                && let Some((smi, met_mol)) = cleave_bond_keeping_side(mol, idx, n_idx, n_idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "N-Dealkylation -> Dealkylated Amine".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: false,
                    bioactivation_reason: None,
                });
            }
        } else if rxn.contains("Heterocycle Epoxidation") {
            if let Some((smi, met_mol)) = epoxidize_five_ring_heterocycle(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "Furan/Thiophene 2,3-Epoxidation".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: true,
                    bioactivation_reason: Some(
                        "Reactive furan/thiophene epoxide (cis-enedial / covalent DILI bioactivation)"
                            .to_string(),
                    ),
                });
            }
        } else if rxn.contains("Terminal Alkene Epoxidation") {
            if let Some((smi, met_mol)) = epoxidize_alkene(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "Alkene Epoxidation -> Oxirane".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: true,
                    bioactivation_reason: Some(
                        "Reactive electrophilic epoxide (GSH trapping / DNA & protein adduct risk)"
                            .to_string(),
                    ),
                });
            }
        } else if rxn.contains("Aniline N-Hydroxylation") {
            if let Some((smi, met_mol)) = attach_hydroxyl(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "Aniline N-Hydroxylation -> Ar-NHOH".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: true,
                    bioactivation_reason: Some(
                        "Aromatic hydroxylamine / nitroso bioactivation (AMES mutagenicity & DILI)"
                            .to_string(),
                    ),
                });
            }
        } else if rxn.contains("Thioether S-Oxidation") {
            if let Some((smi, met_mol)) = oxidize_sulfur(mol, idx, 1)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "Thioether S-Oxidation -> Sulfoxide".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: false,
                    bioactivation_reason: None,
                });
            }
            if let Some((smi, met_mol)) = oxidize_sulfur(mol, idx, 2)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseI,
                    reaction: "Thioether S-Oxidation -> Sulfone".to_string(),
                    enzyme: site.isoform.clone(),
                    parent_som_atom: idx,
                    lability_score: (site.lability_score * 0.90).min(1.0),
                    reactive_bioactivation: false,
                    bioactivation_reason: None,
                });
            }
        } else if rxn.contains("Ester Hydrolysis") {
            if let Some((acid_res, alc_res)) = hydrolyze_ester(mol, idx) {
                if let Some((smi, met_mol)) = acid_res
                    && seen_smiles.insert(smi.clone())
                {
                    raw_cands.push(RawMetaboliteCandidate {
                        mol: met_mol,
                        smiles: smi,
                        phase: MetabolismPhase::PhaseI,
                        reaction: "Ester Hydrolysis -> Carboxylic Acid".to_string(),
                        enzyme: site.isoform.clone(),
                        parent_som_atom: idx,
                        lability_score: site.lability_score,
                        reactive_bioactivation: false,
                        bioactivation_reason: None,
                    });
                }
                if let Some((smi, met_mol)) = alc_res
                    && seen_smiles.insert(smi.clone())
                {
                    raw_cands.push(RawMetaboliteCandidate {
                        mol: met_mol,
                        smiles: smi,
                        phase: MetabolismPhase::PhaseI,
                        reaction: "Ester Hydrolysis -> Alcohol / Phenol".to_string(),
                        enzyme: site.isoform.clone(),
                        parent_som_atom: idx,
                        lability_score: (site.lability_score * 0.95).min(1.0),
                        reactive_bioactivation: false,
                        bioactivation_reason: None,
                    });
                }
            }
        } else if rxn.contains("Phenolic O-Glucuronidation") {
            if let Some((smi, met_mol)) = attach_glucuronide(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseII,
                    reaction: "Phenolic O-Glucuronidation (UGT)".to_string(),
                    enzyme: "UGT1A/2B (Phase II)".to_string(),
                    parent_som_atom: idx,
                    lability_score: site.lability_score,
                    reactive_bioactivation: false,
                    bioactivation_reason: None,
                });
            }
            if let Some((smi, met_mol)) = attach_sulfate(mol, idx)
                && seen_smiles.insert(smi.clone())
            {
                raw_cands.push(RawMetaboliteCandidate {
                    mol: met_mol,
                    smiles: smi,
                    phase: MetabolismPhase::PhaseII,
                    reaction: "Phenolic O-Sulfation (SULT)".to_string(),
                    enzyme: "SULT1A1 (Phase II)".to_string(),
                    parent_som_atom: idx,
                    lability_score: (site.lability_score * 0.92).min(1.0),
                    reactive_bioactivation: false,
                    bioactivation_reason: None,
                });
            }
        } else if rxn.contains("Acyl-Glucuronidation")
            && let Some(oh_idx) = find_acid_oh_oxygen(mol, idx)
            && let Some((smi, met_mol)) = attach_glucuronide(mol, oh_idx)
            && seen_smiles.insert(smi.clone())
        {
            raw_cands.push(RawMetaboliteCandidate {
                mol: met_mol,
                smiles: smi,
                phase: MetabolismPhase::PhaseII,
                reaction: "1-beta-O-Acyl-Glucuronidation (UGT)".to_string(),
                enzyme: site.isoform.clone(),
                parent_som_atom: idx,
                lability_score: site.lability_score,
                reactive_bioactivation: true,
                bioactivation_reason: Some(
                    "Reactive 1-beta-O-acyl-glucuronide (acyl migration / hepatic protein adduct DILI risk)"
                        .to_string(),
                ),
            });
        }
    }

    // Direct check for para-aminophenol / acetaminophen-like quinone-imine oxidation (NAPQI pathway)
    if let Some((smi, met_mol, atom_idx)) = oxidize_para_aminophenol_to_quinone_imine(mol)
        && seen_smiles.insert(smi.clone())
    {
        raw_cands.push(RawMetaboliteCandidate {
            mol: met_mol,
            smiles: smi,
            phase: MetabolismPhase::PhaseI,
            reaction: "2e- Oxidation -> Quinone-Imine (NAPQI-type)".to_string(),
            enzyme: "CYP2E1 / CYP1A2 / CYP3A4".to_string(),
            parent_som_atom: atom_idx,
            lability_score: 0.85,
            reactive_bioactivation: true,
            bioactivation_reason: Some(
                "Electrophilic quinone-imine bioactivation (GSH depletion & acute hepatotoxicity / DILI)"
                    .to_string(),
            ),
        });
    }

    raw_cands.sort_by(|a, b| {
        b.reactive_bioactivation
            .cmp(&a.reactive_bioactivation)
            .then_with(|| b.lability_score.total_cmp(&a.lability_score))
            .then_with(|| a.smiles.cmp(&b.smiles))
    });
    raw_cands.truncate(cap);

    let parent_alerts = collect_alert_set(parent_eval);
    let parent_dili = qsar_prob(parent_eval, "dili");
    let parent_ames = qsar_prob(parent_eval, "ames");
    let parent_herg = qsar_prob(parent_eval, "herg_karim");

    let mut metabolites = Vec::with_capacity(raw_cands.len());
    let mut bioact_cnt = 0usize;
    let mut esc_cnt = 0usize;

    for (rank_idx, cand) in raw_cands.into_iter().enumerate() {
        let mut met_eval = evaluate_molecule(&cand.mol);
        if let Some(mp) = pack {
            let feat = extract_mol_features(&cand.mol, &met_eval);
            met_eval.qsar = Some(mp.predict_all(&feat, "metabolite-guard"));
        }

        let met_alerts = collect_alert_set(&met_eval);
        let new_alerts: Vec<String> = met_alerts.difference(&parent_alerts).cloned().collect();

        let dili_p = qsar_prob(&met_eval, "dili");
        let ames_p = qsar_prob(&met_eval, "ames");
        let herg_p = qsar_prob(&met_eval, "herg_karim");
        let max_tox21_p = max_tox21_prob(&met_eval);

        let qsar_jump = jumped_high(dili_p, parent_dili)
            || jumped_high(ames_p, parent_ames)
            || jumped_high(herg_p, parent_herg);

        let toxicity_escalation =
            cand.reactive_bioactivation || !new_alerts.is_empty() || qsar_jump;

        if cand.reactive_bioactivation {
            bioact_cnt += 1;
        }
        if toxicity_escalation {
            esc_cnt += 1;
        }

        metabolites.push(PredictedMetabolite {
            rank: rank_idx + 1,
            smiles: cand.smiles,
            phase: cand.phase,
            reaction: cand.reaction,
            enzyme: cand.enzyme,
            parent_som_atom: cand.parent_som_atom,
            lability_score: cand.lability_score,
            reactive_bioactivation: cand.reactive_bioactivation,
            bioactivation_reason: cand.bioactivation_reason,
            dili_prob: dili_p.map(round3),
            ames_prob: ames_p.map(round3),
            herg_prob: herg_p.map(round3),
            max_tox21_prob: max_tox21_p.map(round3),
            new_alerts,
            toxicity_escalation,
        });
    }

    MetaboliteGuardReport {
        total_generated: metabolites.len(),
        bioactivation_alerts: bioact_cnt,
        escalation_count: esc_cnt,
        metabolites,
    }
}

#[inline]
fn round3(v: f64) -> f64 {
    (v * 1000.0).round() / 1000.0
}

fn jumped_high(met: Option<f64>, parent: Option<f64>) -> bool {
    match (met, parent) {
        (Some(m), Some(p)) => m >= 0.65 && (m - p) >= 0.18,
        _ => false,
    }
}

fn qsar_prob(eval: &MoleculeEvaluation, id: &str) -> Option<f64> {
    eval.qsar.as_ref().and_then(|q| {
        q.endpoints
            .iter()
            .find(|e| e.endpoint_id == id)
            .map(|e| e.value)
    })
}

fn max_tox21_prob(eval: &MoleculeEvaluation) -> Option<f64> {
    let q = eval.qsar.as_ref()?;
    let mut mx: Option<f64> = None;
    for ep in &q.endpoints {
        if ep.endpoint_id.starts_with("tox21_") {
            mx = Some(mx.map_or(ep.value, |cur| cur.max(ep.value)));
        }
    }
    mx
}

fn collect_alert_set(eval: &MoleculeEvaluation) -> BTreeSet<String> {
    let mut set = BTreeSet::new();
    for a in &eval.medicinal_chemistry.pains_alerts {
        set.insert(format!("PAINS:{a}"));
    }
    for a in &eval.medicinal_chemistry.brenk_alerts {
        set.insert(format!("Brenk:{a}"));
    }
    for a in &eval.safety.ames_alerts {
        set.insert(format!("AMES:{a}"));
    }
    for a in &eval.safety.dili_alerts {
        set.insert(format!("DILI:{a}"));
    }
    set
}

fn finalize_metabolite(mut m: Molecule) -> Option<(String, Molecule)> {
    m.perceive().ok()?;
    let smi = m.try_to_canonical_smiles().ok()?;
    if smi.is_empty() {
        return None;
    }
    let reparsed = parse_smiles(&smi).ok()?;
    if reparsed.heavy_atom_count() < 4 {
        return None;
    }
    Some((reparsed.to_canonical_smiles(), reparsed))
}

fn reset_atom_h_flags(m: &mut Molecule, idx: usize) {
    if let Some(a) = m.atoms.get_mut(idx) {
        a.explicit_h = 0;
        a.from_bracket = false;
    }
}

fn attach_hydroxyl(mol: &Molecule, target_idx: usize) -> Option<(String, Molecule)> {
    let mut m = mol.clone();
    reset_atom_h_flags(&mut m, target_idx);
    let o_atom = Atom::new(0, "O").ok()?;
    let o_idx = m.add_atom(o_atom);
    m.add_bond(target_idx, o_idx, BondOrder::Single, BondStereo::None)
        .ok()?;
    finalize_metabolite(m)
}

fn find_neighbor_by_z(mol: &Molecule, idx: usize, z: u8) -> Option<usize> {
    mol.adj.get(idx)?.iter().find_map(|&(nb, _)| {
        mol.atoms
            .get(nb)
            .filter(|a| a.atomic_number == z)
            .map(|_| nb)
    })
}

fn find_acid_oh_oxygen(mol: &Molecule, c_idx: usize) -> Option<usize> {
    mol.adj.get(c_idx)?.iter().find_map(|&(nb, bid)| {
        let a = mol.atoms.get(nb)?;
        let b = mol.bonds.get(bid)?;
        if a.atomic_number == 8
            && b.order == BondOrder::Single
            && mol.adj.get(nb).is_some_and(|e| e.len() == 1)
        {
            Some(nb)
        } else {
            None
        }
    })
}

/// Cleaves covalent bond `(u, v)` and returns the connected component containing `keep_atom`.
fn cleave_bond_keeping_side(
    mol: &Molecule,
    u: usize,
    v: usize,
    keep_atom: usize,
) -> Option<(String, Molecule)> {
    let frag = extract_component_after_cleave(mol, u, v, keep_atom)?;
    finalize_metabolite(frag)
}

fn extract_component_after_cleave(
    mol: &Molecule,
    u: usize,
    v: usize,
    seed: usize,
) -> Option<Molecule> {
    let n = mol.atoms.len();
    if seed >= n {
        return None;
    }
    let mut visited = vec![false; n];
    let mut q = VecDeque::new();
    visited[seed] = true;
    q.push_back(seed);
    let mut comp = Vec::new();
    while let Some(cur) = q.pop_front() {
        comp.push(cur);
        for &(nb, _) in &mol.adj[cur] {
            if (cur == u && nb == v) || (cur == v && nb == u) {
                continue;
            }
            if !visited[nb] {
                visited[nb] = true;
                q.push_back(nb);
            }
        }
    }
    comp.sort_unstable();
    let mut old_to_new = vec![usize::MAX; n];
    let mut frag = Molecule::new();
    for &old_i in &comp {
        let mut a = mol.atoms[old_i].clone();
        if old_i == u || old_i == v {
            a.explicit_h = 0;
            a.from_bracket = false;
            if a.atomic_number == 8 && a.formal_charge < 0 {
                a.formal_charge = 0;
            }
        }
        let new_i = frag.add_atom(a);
        old_to_new[old_i] = new_i;
    }
    for b in &mol.bonds {
        if (b.u == u && b.v == v) || (b.u == v && b.v == u) {
            continue;
        }
        let nu = old_to_new[b.u];
        let nv = old_to_new[b.v];
        if nu != usize::MAX && nv != usize::MAX {
            frag.add_bond(nu, nv, b.order, b.stereo).ok()?;
        }
    }
    Some(frag)
}

fn epoxidize_alkene(mol: &Molecule, c_idx: usize) -> Option<(String, Molecule)> {
    let (other_c, bid) = mol.adj.get(c_idx)?.iter().find_map(|&(nb, bid)| {
        let nba = mol.atoms.get(nb)?;
        let b = mol.bonds.get(bid)?;
        if nba.atomic_number == 6 && b.order == BondOrder::Double {
            Some((nb, bid))
        } else {
            None
        }
    })?;
    let mut m = mol.clone();
    if let Some(b) = m.bonds.get_mut(bid) {
        b.order = BondOrder::Single;
    }
    reset_atom_h_flags(&mut m, c_idx);
    reset_atom_h_flags(&mut m, other_c);
    let o_idx = m.add_atom(Atom::new(0, "O").ok()?);
    m.add_bond(c_idx, o_idx, BondOrder::Single, BondStereo::None)
        .ok()?;
    m.add_bond(other_c, o_idx, BondOrder::Single, BondStereo::None)
        .ok()?;
    finalize_metabolite(m)
}

fn epoxidize_five_ring_heterocycle(mol: &Molecule, c_idx: usize) -> Option<(String, Molecule)> {
    let ring = mol.rings.iter().find(|r| {
        r.len() == 5
            && r.contains(&c_idx)
            && r.iter().any(|&ri| {
                mol.atoms
                    .get(ri)
                    .is_some_and(|a| a.atomic_number == 8 || a.atomic_number == 16)
            })
    })?;
    let het_pos = ring.iter().position(|&ri| {
        mol.atoms
            .get(ri)
            .is_some_and(|a| a.atomic_number == 8 || a.atomic_number == 16)
    })?;
    // Ordered ring starting at heteroatom: [X, C2, C3, C4, C5]
    let mut ord = [0usize; 5];
    for k in 0..5 {
        ord[k] = ring[(het_pos + k) % 5];
    }
    // Choose which double bond to epoxidize (C2-C3 if c_idx in {C2, C3}, else C4-C5)
    let (ep_a, ep_b, dbl_a, dbl_b) = if c_idx == ord[1] || c_idx == ord[2] {
        (ord[1], ord[2], ord[3], ord[4])
    } else {
        (ord[3], ord[4], ord[1], ord[2])
    };

    let mut m = mol.clone();
    for &ri in &ord {
        if let Some(a) = m.atoms.get_mut(ri) {
            a.aromatic = false;
            a.explicit_h = 0;
            a.from_bracket = false;
        }
    }
    for k in 0..5 {
        let u = ord[k];
        let v = ord[(k + 1) % 5];
        if let Some(bid) = m
            .adj
            .get(u)
            .and_then(|e| e.iter().find(|&&(nb, _)| nb == v).map(|&(_, b)| b))
            && let Some(b) = m.bonds.get_mut(bid)
        {
            b.order = if (u == dbl_a && v == dbl_b) || (u == dbl_b && v == dbl_a) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
        }
    }
    let o_idx = m.add_atom(Atom::new(0, "O").ok()?);
    m.add_bond(ep_a, o_idx, BondOrder::Single, BondStereo::None)
        .ok()?;
    m.add_bond(ep_b, o_idx, BondOrder::Single, BondStereo::None)
        .ok()?;
    finalize_metabolite(m)
}

fn oxidize_sulfur(mol: &Molecule, s_idx: usize, oxygens: usize) -> Option<(String, Molecule)> {
    let mut m = mol.clone();
    reset_atom_h_flags(&mut m, s_idx);
    for _ in 0..oxygens {
        let o_idx = m.add_atom(Atom::new(0, "O").ok()?);
        m.add_bond(s_idx, o_idx, BondOrder::Double, BondStereo::None)
            .ok()?;
    }
    finalize_metabolite(m)
}

type HydrolysisPair = (Option<(String, Molecule)>, Option<(String, Molecule)>);

fn hydrolyze_ester(mol: &Molecule, c_idx: usize) -> Option<HydrolysisPair> {
    let o_ester = mol.adj.get(c_idx)?.iter().find_map(|&(nb, bid)| {
        let a = mol.atoms.get(nb)?;
        let b = mol.bonds.get(bid)?;
        if a.atomic_number == 8
            && b.order == BondOrder::Single
            && mol.adj.get(nb).is_some_and(|e| e.len() >= 2)
        {
            Some(nb)
        } else {
            None
        }
    })?;

    // Acyl fragment -> add -OH back onto the carbonyl carbon
    let acid_side = extract_component_after_cleave(mol, c_idx, o_ester, c_idx).and_then(|mut f| {
        // Locate carbonyl carbon in the new component (it has a double-bonded O and degree <= 2)
        let carb = f.atoms.iter().enumerate().find_map(|(i, a)| {
            if a.atomic_number == 6
                && f.adj[i].iter().any(|&(nb, bid)| {
                    f.atoms.get(nb).is_some_and(|x| x.atomic_number == 8)
                        && f.bonds
                            .get(bid)
                            .is_some_and(|b| b.order == BondOrder::Double)
                })
                && !f.adj[i].iter().any(|&(nb, bid)| {
                    f.atoms.get(nb).is_some_and(|x| x.atomic_number == 8)
                        && f.bonds
                            .get(bid)
                            .is_some_and(|b| b.order == BondOrder::Single)
                })
            {
                Some(i)
            } else {
                None
            }
        })?;
        let oh = f.add_atom(Atom::new(0, "O").ok()?);
        f.add_bond(carb, oh, BondOrder::Single, BondStereo::None)
            .ok()?;
        finalize_metabolite(f)
    });

    // Alcohol/phenol leaving group fragment
    let alc_side = cleave_bond_keeping_side(mol, c_idx, o_ester, o_ester);
    Some((acid_side, alc_side))
}

/// Attaches a β-D-glucuronide pyranose ring (`-C1OC(C(=O)O)C(O)C(O)C1O`) onto oxygen `o_idx`.
fn attach_glucuronide(mol: &Molecule, o_idx: usize) -> Option<(String, Molecule)> {
    let mut m = mol.clone();
    reset_atom_h_flags(&mut m, o_idx);
    if let Some(oa) = m.atoms.get_mut(o_idx) {
        oa.formal_charge = 0;
    }
    // Build pyranose ring: C1 - O5 - C5(C6(=O)O) - C4(O) - C3(O) - C2(O) - C1
    let c1 = m.add_atom(Atom::new(0, "C").ok()?);
    let c2 = m.add_atom(Atom::new(0, "C").ok()?);
    let c3 = m.add_atom(Atom::new(0, "C").ok()?);
    let c4 = m.add_atom(Atom::new(0, "C").ok()?);
    let c5 = m.add_atom(Atom::new(0, "C").ok()?);
    let o5 = m.add_atom(Atom::new(0, "O").ok()?);
    let o2 = m.add_atom(Atom::new(0, "O").ok()?);
    let o3 = m.add_atom(Atom::new(0, "O").ok()?);
    let o4 = m.add_atom(Atom::new(0, "O").ok()?);
    let c6 = m.add_atom(Atom::new(0, "C").ok()?);
    let o6a = m.add_atom(Atom::new(0, "O").ok()?);
    let o6b = m.add_atom(Atom::new(0, "O").ok()?);

    m.add_bond(o_idx, c1, BondOrder::Single, BondStereo::None)
        .ok()?;
    for (u, v) in [
        (c1, c2),
        (c2, c3),
        (c3, c4),
        (c4, c5),
        (c5, o5),
        (o5, c1),
        (c2, o2),
        (c3, o3),
        (c4, o4),
        (c5, c6),
        (c6, o6b),
    ] {
        m.add_bond(u, v, BondOrder::Single, BondStereo::None).ok()?;
    }
    m.add_bond(c6, o6a, BondOrder::Double, BondStereo::None)
        .ok()?;
    finalize_metabolite(m)
}

/// Attaches a sulfate ester group (`-S(=O)(=O)O`) onto oxygen `o_idx`.
fn attach_sulfate(mol: &Molecule, o_idx: usize) -> Option<(String, Molecule)> {
    let mut m = mol.clone();
    reset_atom_h_flags(&mut m, o_idx);
    if let Some(oa) = m.atoms.get_mut(o_idx) {
        oa.formal_charge = 0;
    }
    let s = m.add_atom(Atom::new(0, "S").ok()?);
    let o1 = m.add_atom(Atom::new(0, "O").ok()?);
    let o2 = m.add_atom(Atom::new(0, "O").ok()?);
    let o3 = m.add_atom(Atom::new(0, "O").ok()?);
    m.add_bond(o_idx, s, BondOrder::Single, BondStereo::None)
        .ok()?;
    m.add_bond(s, o1, BondOrder::Double, BondStereo::None)
        .ok()?;
    m.add_bond(s, o2, BondOrder::Double, BondStereo::None)
        .ok()?;
    m.add_bond(s, o3, BondOrder::Single, BondStereo::None)
        .ok()?;
    finalize_metabolite(m)
}

/// Checks if `mol` contains a 6-membered aromatic ring with ortho (1,2) or para (1,4)
/// `-OH` / `-NH-` / `-NH2` substituents (catechol, hydroquinone, or p-aminophenol bioactivation).
fn has_quinone_or_quinone_imine_precursor(mol: &Molecule) -> bool {
    for ring in &mol.rings {
        if ring.len() != 6 {
            continue;
        }
        if !ring.iter().all(|&ri| {
            mol.atoms
                .get(ri)
                .is_some_and(|a| a.aromatic && a.atomic_number == 6)
        }) {
            continue;
        }
        let mut donor_positions = Vec::new();
        for (pos, &ri) in ring.iter().enumerate() {
            let has_donor = mol.adj[ri].iter().any(|&(nb, _)| {
                !ring.contains(&nb)
                    && mol.atoms.get(nb).is_some_and(|a| {
                        (a.atomic_number == 8 && a.total_h() >= 1)
                            || (a.atomic_number == 7 && a.total_h() >= 1 && !a.aromatic)
                    })
            });
            if has_donor {
                donor_positions.push(pos);
            }
        }
        for i in 0..donor_positions.len() {
            for j in (i + 1)..donor_positions.len() {
                let diff = donor_positions[j] - donor_positions[i];
                if diff == 1 || diff == 3 || diff == 5 {
                    return true;
                }
            }
        }
    }
    false
}

/// Detects a *para*-aminophenol or *para*-acetamidophenol motif (like Acetaminophen) and
/// generates its 2e- oxidized reactive **quinone-imine** metabolite (NAPQI-like).
fn oxidize_para_aminophenol_to_quinone_imine(mol: &Molecule) -> Option<(String, Molecule, usize)> {
    for ring in &mol.rings {
        if ring.len() != 6 {
            continue;
        }
        if !ring.iter().all(|&ri| {
            mol.atoms
                .get(ri)
                .is_some_and(|a| a.aromatic && a.atomic_number == 6)
        }) {
            continue;
        }
        for pos in 0..6 {
            let c_oh = ring[pos];
            let c_n = ring[(pos + 3) % 6];
            let o_idx = mol.adj[c_oh].iter().find_map(|&(nb, _)| {
                if ring.contains(&nb) {
                    return None;
                }
                let a = mol.atoms.get(nb)?;
                if a.atomic_number == 8 && a.total_h() >= 1 && mol.adj[nb].len() == 1 {
                    Some(nb)
                } else {
                    None
                }
            });
            let n_idx = mol.adj[c_n].iter().find_map(|&(nb, _)| {
                if ring.contains(&nb) {
                    return None;
                }
                let a = mol.atoms.get(nb)?;
                if a.atomic_number == 7
                    && !a.aromatic
                    && a.total_h() >= 1
                    && a.hybridization != Hybridization::Sp
                {
                    Some(nb)
                } else {
                    None
                }
            });
            if let (Some(oi), Some(ni)) = (o_idx, n_idx) {
                let mut m = mol.clone();
                // De-aromatize the 6-membered ring into cyclohexa-2,5-diene-1,4-imine-one
                let c1 = ring[pos];
                let c2 = ring[(pos + 1) % 6];
                let c3 = ring[(pos + 2) % 6];
                let c4 = ring[(pos + 3) % 6];
                let c5 = ring[(pos + 4) % 6];
                let c6 = ring[(pos + 5) % 6];
                for &ri in &[c1, c2, c3, c4, c5, c6, oi, ni] {
                    if let Some(a) = m.atoms.get_mut(ri) {
                        a.aromatic = false;
                        a.explicit_h = 0;
                        a.from_bracket = false;
                    }
                }
                for (u, v, ord) in [
                    (c1, oi, BondOrder::Double),
                    (c4, ni, BondOrder::Double),
                    (c1, c2, BondOrder::Single),
                    (c2, c3, BondOrder::Double),
                    (c3, c4, BondOrder::Single),
                    (c4, c5, BondOrder::Single),
                    (c5, c6, BondOrder::Double),
                    (c6, c1, BondOrder::Single),
                ] {
                    if let Some(bid) = m
                        .adj
                        .get(u)
                        .and_then(|e| e.iter().find(|&&(nb, _)| nb == v).map(|&(_, b)| b))
                        && let Some(b) = m.bonds.get_mut(bid)
                    {
                        b.order = ord;
                    }
                }
                if let Some((smi, met)) = finalize_metabolite(m) {
                    return Some((smi, met, ni));
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::som::predict_metabolism_sites;

    #[test]
    fn metabolite_guard_captures_ibuprofen_paracetamol_furan_and_aniline_bioactivation()
    -> Result<(), Box<dyn std::error::Error>> {
        // 1. Ibuprofen -> Benzylic hydroxylation + reactive 1-beta-O-acyl-glucuronide
        let ibu = parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")?;
        let ibu_ev = evaluate_molecule(&ibu);
        let ibu_sites = predict_metabolism_sites(&ibu);
        let ibu_rep = run_metabolite_guard(&ibu, &ibu_ev, &ibu_sites, None, 10);
        assert!(
            ibu_rep.total_generated >= 2,
            "expected >= 2 Ibuprofen metabolites, got {ibu_rep:?}"
        );
        assert!(
            ibu_rep
                .metabolites
                .iter()
                .any(|m| m.reaction.contains("Benzylic")),
            "expected benzylic hydroxylated metabolite for Ibuprofen"
        );
        assert!(
            ibu_rep
                .metabolites
                .iter()
                .any(|m| m.reaction.contains("Acyl-Glucuronidation") && m.reactive_bioactivation),
            "expected acyl-glucuronide bioactivation alert for Ibuprofen"
        );

        // 2. Paracetamol -> Phenolic O-glucuronide + O-sulfate + reactive NAPQI quinone-imine!
        let para = parse_smiles("CC(=O)Nc1ccc(O)cc1")?;
        let para_ev = evaluate_molecule(&para);
        let para_sites = predict_metabolism_sites(&para);
        let para_rep = run_metabolite_guard(&para, &para_ev, &para_sites, None, 10);
        assert!(
            para_rep
                .metabolites
                .iter()
                .any(|m| m.reaction.contains("Quinone-Imine") && m.reactive_bioactivation),
            "expected NAPQI quinone-imine bioactivation for Paracetamol, got {:?}",
            para_rep.metabolites
        );
        assert!(
            para_rep
                .metabolites
                .iter()
                .any(|m| m.phase == MetabolismPhase::PhaseII
                    && m.reaction.contains("O-Glucuronidation")),
            "expected Phase II O-glucuronide for Paracetamol"
        );

        // 3. Furan derivative -> reactive 2,3-epoxide bioactivation!
        let furan_drug = parse_smiles("O=C(Nc1ccccc1)c1ccco1")?;
        let fur_ev = evaluate_molecule(&furan_drug);
        let fur_sites = predict_metabolism_sites(&furan_drug);
        let fur_rep = run_metabolite_guard(&furan_drug, &fur_ev, &fur_sites, None, 10);
        assert!(
            fur_rep.bioactivation_alerts >= 1 && fur_rep.escalation_count >= 1,
            "expected furan epoxide bioactivation & toxicity escalation, got {fur_rep:?}"
        );

        // 4. Primary aniline -> reactive N-hydroxylamine Ar-NHOH
        let anil = parse_smiles("Nc1ccc(Cl)cc1")?;
        let anil_ev = evaluate_molecule(&anil);
        let anil_sites = predict_metabolism_sites(&anil);
        let anil_rep = run_metabolite_guard(&anil, &anil_ev, &anil_sites, None, 10);
        assert!(
            anil_rep
                .metabolites
                .iter()
                .any(|m| m.reaction.contains("N-Hydroxylation") && m.reactive_bioactivation),
            "expected aniline N-hydroxylation bioactivation, got {anil_rep:?}"
        );

        // 5. Propranolol -> O-dealkylation & N-dealkylation
        let prop = parse_smiles("CC(C)NCC(O)COc1cccc2ccccc12")?;
        let prop_ev = evaluate_molecule(&prop);
        let prop_sites = predict_metabolism_sites(&prop);
        let prop_rep = run_metabolite_guard(&prop, &prop_ev, &prop_sites, None, 10);
        assert!(
            prop_rep
                .metabolites
                .iter()
                .any(|m| m.reaction.contains("Dealkylation")),
            "expected O-/N-dealkylation metabolite for Propranolol, got {prop_rep:?}"
        );

        Ok(())
    }
}
