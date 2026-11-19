#![forbid(unsafe_code)]

use admet_calc::{MetabolismRisk, MetabolismSite};
use chem_core::{BondOrder, Hybridization, Molecule};

/// Predicts atomic Sites of Metabolism (CYP3A4 / CYP2D6 / CYP2C9 / CYP1A2 + Phase II UGT / Esterase)
/// with actionable medicinal-chemistry blocking recommendations.
#[must_use]
pub fn predict_metabolism_sites(mol: &Molecule) -> Vec<MetabolismSite> {
    let mut sites: Vec<MetabolismSite> = Vec::new();
    let has_basic_n = mol
        .atoms
        .iter()
        .enumerate()
        .any(|(i, a)| is_basic_amine(mol, i, a.atomic_number));

    for (idx, atom) in mol.atoms.iter().enumerate() {
        if atom.atomic_number <= 1 {
            continue;
        }
        let mut best: Option<MetabolismSite> = None;
        let mut consider = |isoform: &str, reaction: &str, score: f64, fix: &str| {
            let sc = (score.clamp(0.0, 1.0) * 100.0).round() / 100.0;
            if sc < 0.45 {
                return;
            }
            if best.as_ref().is_none_or(|b| sc > b.lability_score) {
                let risk = if sc >= 0.75 {
                    MetabolismRisk::High
                } else if sc >= 0.55 {
                    MetabolismRisk::Moderate
                } else {
                    MetabolismRisk::Low
                };
                best = Some(MetabolismSite {
                    atom_index: idx,
                    atom_symbol: atom.symbol.clone(),
                    isoform: isoform.to_string(),
                    reaction: reaction.to_string(),
                    lability_score: sc,
                    risk,
                    medchem_fix: fix.to_string(),
                });
            }
        };

        let h_cnt = atom.total_h();

        match atom.atomic_number {
            6 => {
                // Count fluorinated/halogenated neighbors (electronic/steric blocking)
                let f_nbrs = mol.adj[idx]
                    .iter()
                    .filter(|&&(nb, _)| mol.atoms.get(nb).is_some_and(|x| x.atomic_number == 9))
                    .count();
                let f_damp = 1.0 - 0.35 * (f_nbrs as f64).min(2.5);

                if !atom.aromatic && atom.hybridization == Hybridization::Sp3 && h_cnt >= 1 {
                    let mut adj_arom = false;
                    let mut adj_alkene = false;
                    let mut adj_basic_n = false;
                    let mut adj_other_n = false;
                    let mut adj_ether_o = false;
                    let mut heavy_deg = 0usize;

                    for &(nb, bid) in &mol.adj[idx] {
                        let Some(nba) = mol.atoms.get(nb) else {
                            continue;
                        };
                        if nba.atomic_number <= 1 {
                            continue;
                        }
                        heavy_deg += 1;
                        if nba.aromatic && nba.atomic_number == 6 {
                            adj_arom = true;
                        }
                        if nba.atomic_number == 6
                            && !nba.aromatic
                            && mol.adj[nb].iter().any(|&(nb2, bid2)| {
                                nb2 != idx
                                    && mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 6)
                                    && mol
                                        .bonds
                                        .get(bid2)
                                        .is_some_and(|b| b.order == BondOrder::Double)
                            })
                        {
                            adj_alkene = true;
                        }
                        if nba.atomic_number == 7 && !is_amide_nitrogen(mol, nb) {
                            if is_basic_amine(mol, nb, nba.atomic_number) {
                                adj_basic_n = true;
                            } else {
                                adj_other_n = true;
                            }
                        }
                        if nba.atomic_number == 8
                            && mol
                                .bonds
                                .get(bid)
                                .is_some_and(|b| b.order == BondOrder::Single)
                            && nba.total_h() == 0
                            && !is_ester_oxygen(mol, nb)
                        {
                            adj_ether_o = true;
                        }
                    }

                    if adj_basic_n {
                        let bonus = if h_cnt >= 2 { 0.05 } else { 0.0 };
                        consider(
                            "CYP2D6 / CYP3A4",
                            "N-Dealkylation (alpha-C-H)",
                            (0.86 + bonus) * f_damp,
                            "Replace N-alkyl with N-CD3, N-CF3, N-cyclopropyl, or cap as amide",
                        );
                    } else if adj_other_n {
                        consider(
                            "CYP3A4",
                            "N-Dealkylation (alpha-C-H)",
                            0.78 * f_damp,
                            "Replace N-Me/Et with N-cyclopropyl, CF3, or morpholine/azetidine",
                        );
                    }

                    if adj_ether_o {
                        let iso = if has_basic_n {
                            "CYP2D6 / CYP3A4"
                        } else {
                            "CYP1A2 / CYP3A4"
                        };
                        let sc = if h_cnt == 3 { 0.86 } else { 0.80 };
                        consider(
                            iso,
                            "O-Dealkylation (alpha-C-H)",
                            sc * f_damp,
                            "Replace -OCH3/-OR with -OCF3, -OCHF2, or halogen",
                        );
                    }

                    if adj_arom {
                        let sc = if h_cnt >= 2 { 0.88 } else { 0.84 };
                        consider(
                            "CYP3A4 / CYP2C9",
                            "Benzylic C-H Hydroxylation",
                            sc * f_damp,
                            "Block benzylic position with CF2, O, gem-dimethyl, or cyclopropyl",
                        );
                    } else if adj_alkene {
                        consider(
                            "CYP3A4",
                            "Allylic C-H Hydroxylation",
                            0.82 * f_damp,
                            "Saturate alkene or block allylic C-H with fluorine/cyclopropyl",
                        );
                    } else if heavy_deg == 3 && h_cnt == 1 && !atom.in_ring {
                        consider(
                            "CYP3A4",
                            "Tertiary Aliphatic C-H Hydroxylation",
                            0.66 * f_damp,
                            "Replace tertiary C-H with C-F, C-OH, or spiro-oxetane",
                        );
                    }
                } else if atom.aromatic && h_cnt >= 1 {
                    // Check furan/thiophene 5-membered bioactivation ring
                    let in_furan_or_thio = atom.ring_sizes.contains(&5)
                        && mol.rings.iter().any(|r| {
                            r.len() == 5
                                && r.contains(&idx)
                                && r.iter().any(|&ri| {
                                    mol.atoms.get(ri).is_some_and(|x| {
                                        x.atomic_number == 8 || x.atomic_number == 16
                                    })
                                })
                                && !r.iter().any(|&ri| {
                                    mol.atoms.get(ri).is_some_and(|x| x.atomic_number == 7)
                                })
                        });
                    if in_furan_or_thio {
                        consider(
                            "CYP3A4 / CYP1A2",
                            "Heterocycle Epoxidation (Bioactivation)",
                            0.87,
                            "Replace furan/thiophene with oxazole, thiazole, or phenyl ring",
                        );
                    } else if let Some((edg, ewg, has_ring_n)) = ring_electronic_profile(mol, idx) {
                        // Electron-rich aromatic C-H hydroxylation
                        let base = 0.56 + 0.14 * (edg as f64).min(2.0)
                            - 0.15 * (ewg as f64).min(2.0)
                            - if has_ring_n { 0.16 } else { 0.0 };
                        if base >= 0.58 {
                            let iso = if has_basic_n {
                                "CYP2D6 / CYP3A4"
                            } else {
                                "CYP1A2 / CYP2C9"
                            };
                            consider(
                                iso,
                                "Aromatic C-H Hydroxylation",
                                base.min(0.84),
                                "Block para/ortho C-H with F or Cl, or swap CH for ring N (pyridine)",
                            );
                        }
                    }
                } else if !atom.aromatic && atom.hybridization == Hybridization::Sp2 && h_cnt >= 2 {
                    // Terminal alkene =CH2 epoxidation
                    let is_alkene = mol.adj[idx].iter().any(|&(nb, bid)| {
                        mol.atoms.get(nb).is_some_and(|x| x.atomic_number == 6)
                            && mol
                                .bonds
                                .get(bid)
                                .is_some_and(|b| b.order == BondOrder::Double)
                    });
                    if is_alkene {
                        consider(
                            "CYP3A4 / CYP2E1",
                            "Terminal Alkene Epoxidation",
                            0.85,
                            "Hydrogenate double bond or replace vinyl with cyclopropyl/ethyl",
                        );
                    }
                } else if !atom.aromatic && is_ester_carbonyl(mol, idx) {
                    consider(
                        "Plasma Esterase (CES1/2)",
                        "Ester Hydrolysis",
                        0.79,
                        "Replace ester with amide, 1,2,4-oxadiazole, or sterically hindered ester",
                    );
                } else if !atom.aromatic && is_carboxylic_acid_carbon(mol, idx) {
                    consider(
                        "UGT (Phase II)",
                        "Acyl-Glucuronidation",
                        0.66,
                        "Consider tetrazole, acylsulfonamide, or alpha-steric substitution",
                    );
                }
            }
            8 => {
                // Phenolic OH -> Phase II UGT / SULT
                if h_cnt >= 1 {
                    let on_arom = mol.adj[idx]
                        .iter()
                        .any(|&(nb, _)| mol.atoms.get(nb).is_some_and(|x| x.aromatic));
                    if on_arom {
                        consider(
                            "UGT / SULT (Phase II)",
                            "Phenolic O-Glucuronidation / Sulfation",
                            0.83,
                            "Add ortho-F/Me shielding, or replace phenol with indazole/azaindole isostere",
                        );
                    }
                }
            }
            7 => {
                // Primary aniline N-hydroxylation / NAT acetylation
                if !atom.aromatic && h_cnt >= 2 {
                    let on_arom = mol.adj[idx]
                        .iter()
                        .any(|&(nb, _)| mol.atoms.get(nb).is_some_and(|x| x.aromatic));
                    if on_arom {
                        consider(
                            "CYP1A2 / NAT2",
                            "Aniline N-Hydroxylation / Bioactivation",
                            0.80,
                            "Replace free aniline with amino-pyridine, amide, or heteroaromatic isostere",
                        );
                    }
                }
            }
            16 => {
                // Divalent thioether S-oxidation
                if !atom.aromatic && h_cnt == 0 {
                    let n_o = mol.adj[idx]
                        .iter()
                        .filter(|&&(nb, _)| mol.atoms.get(nb).is_some_and(|x| x.atomic_number == 8))
                        .count();
                    let n_c = mol.adj[idx]
                        .iter()
                        .filter(|&&(nb, _)| mol.atoms.get(nb).is_some_and(|x| x.atomic_number == 6))
                        .count();
                    if n_o == 0 && n_c == 2 {
                        consider(
                            "CYP3A4 / FMO",
                            "Thioether S-Oxidation",
                            0.84,
                            "Pre-oxidize to sulfone (-SO2-) or replace thioether with -O- or -CH2-",
                        );
                    }
                }
            }
            _ => {}
        }

        if let Some(site) = best {
            sites.push(site);
        }
    }

    sites.sort_by(|a, b| {
        b.lability_score
            .total_cmp(&a.lability_score)
            .then_with(|| a.atom_index.cmp(&b.atom_index))
    });
    sites.truncate(10);
    sites
}

fn is_amide_nitrogen(mol: &Molecule, idx: usize) -> bool {
    mol.adj.get(idx).is_some_and(|nbrs| {
        nbrs.iter().any(|&(nb, _)| {
            mol.atoms.get(nb).is_some_and(|a| a.atomic_number == 6)
                && mol.adj.get(nb).is_some_and(|cnbrs| {
                    cnbrs.iter().any(|&(nb2, bid2)| {
                        mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 8)
                            && mol
                                .bonds
                                .get(bid2)
                                .is_some_and(|b| b.order == BondOrder::Double)
                    })
                })
        })
    })
}

fn is_basic_amine(mol: &Molecule, idx: usize, z: u8) -> bool {
    if z != 7 {
        return false;
    }
    let Some(a) = mol.atoms.get(idx) else {
        return false;
    };
    if a.aromatic || a.hybridization != Hybridization::Sp3 || is_amide_nitrogen(mol, idx) {
        return false;
    }
    // Exclude anilines and sulfonamides
    !mol.adj[idx].iter().any(|&(nb, _)| {
        mol.atoms
            .get(nb)
            .is_some_and(|x| x.aromatic || x.atomic_number == 16)
    })
}

fn is_ester_oxygen(mol: &Molecule, o_idx: usize) -> bool {
    mol.adj
        .get(o_idx)
        .is_some_and(|nbrs| nbrs.iter().any(|&(nb, _)| is_ester_carbonyl(mol, nb)))
}

fn is_ester_carbonyl(mol: &Molecule, c_idx: usize) -> bool {
    if mol.atoms.get(c_idx).is_none_or(|a| a.atomic_number != 6) {
        return false;
    }
    let mut has_double_o = false;
    let mut has_single_or = false;
    for &(nb, bid) in &mol.adj[c_idx] {
        let Some(nba) = mol.atoms.get(nb) else {
            continue;
        };
        if nba.atomic_number != 8 {
            continue;
        }
        let ord = mol.bonds.get(bid).map(|b| b.order);
        if ord == Some(BondOrder::Double) {
            has_double_o = true;
        } else if ord == Some(BondOrder::Single)
            && nba.total_h() == 0
            && mol.adj[nb].iter().any(|&(nb2, _)| {
                nb2 != c_idx && mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 6)
            })
        {
            has_single_or = true;
        }
    }
    has_double_o && has_single_or
}

fn is_carboxylic_acid_carbon(mol: &Molecule, c_idx: usize) -> bool {
    if mol.atoms.get(c_idx).is_none_or(|a| a.atomic_number != 6) {
        return false;
    }
    let mut has_double_o = false;
    let mut has_oh_or_om = false;
    for &(nb, bid) in &mol.adj[c_idx] {
        let Some(nba) = mol.atoms.get(nb) else {
            continue;
        };
        if nba.atomic_number != 8 {
            continue;
        }
        let ord = mol.bonds.get(bid).map(|b| b.order);
        if ord == Some(BondOrder::Double) {
            has_double_o = true;
        } else if ord == Some(BondOrder::Single)
            && (nba.total_h() >= 1 || nba.formal_charge < 0)
            && mol.adj[nb].len() == 1
        {
            has_oh_or_om = true;
        }
    }
    has_double_o && has_oh_or_om
}

/// Counts electron-donating groups (EDG), electron-withdrawing groups (EWG), and ring nitrogens
/// on the 6-membered aromatic ring containing `atom_idx`.
fn ring_electronic_profile(mol: &Molecule, atom_idx: usize) -> Option<(usize, usize, bool)> {
    let ring = mol
        .rings
        .iter()
        .find(|r| r.len() == 6 && r.contains(&atom_idx))?;
    let mut edg = 0usize;
    let mut ewg = 0usize;
    let mut has_ring_n = false;
    for &ri in ring {
        if mol.atoms.get(ri).is_some_and(|a| a.atomic_number == 7) {
            has_ring_n = true;
        }
        for &(nb, _) in &mol.adj[ri] {
            if ring.contains(&nb) {
                continue;
            }
            let Some(nba) = mol.atoms.get(nb) else {
                continue;
            };
            match nba.atomic_number {
                8 | 7 => {
                    // Check nitro/amide vs donor -OH/-OR/-NR2
                    let is_nitro = nba.atomic_number == 7
                        && mol.adj[nb].iter().any(|&(nb2, bid2)| {
                            mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 8)
                                && mol
                                    .bonds
                                    .get(bid2)
                                    .is_some_and(|b| b.order == BondOrder::Double)
                        });
                    if is_nitro {
                        ewg += 1;
                    } else {
                        edg += 1;
                    }
                }
                9 | 17 | 35 | 53 => ewg += 1,
                16 => ewg += 1,
                6 => {
                    // Check CF3 / CN / C=O electron-withdrawing substituents
                    let n_f = mol.adj[nb]
                        .iter()
                        .filter(|&&(nb2, _)| {
                            mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 9)
                        })
                        .count();
                    let has_multi_het = mol.adj[nb].iter().any(|&(nb2, bid2)| {
                        mol.atoms
                            .get(nb2)
                            .is_some_and(|x| x.atomic_number == 7 || x.atomic_number == 8)
                            && mol.bonds.get(bid2).is_some_and(|b| {
                                b.order == BondOrder::Double || b.order == BondOrder::Triple
                            })
                    });
                    if n_f >= 2 || has_multi_het {
                        ewg += 1;
                    }
                }
                _ => {}
            }
        }
    }
    Some((edg, ewg, has_ring_n))
}
