//! Physiological pKa estimation, dominant ionization state at pH 7.4, and
//! Henderson-Hasselbalch LogD(pH 7.4) calculation for CNS MPO & ADMET profiling.

use chem_core::{BondOrder, Hybridization, Molecule};

use crate::physicochem::round2;
use crate::types::{IonizationClass, IonizationProfile, Lipophilicity};

const PH: f64 = 7.4;

/// Estimate strongest acidic & basic pKa values, physiological charge at pH 7.4,
/// and Henderson-Hasselbalch `logD_7.4`.
pub fn compute_ionization(mol: &Molecule, lipo: &Lipophilicity) -> IonizationProfile {
    let mut min_acid_pka: Option<f64> = None;
    let mut max_base_pka: Option<f64> = None;
    let mut permanent_charge: f64 = 0.0;

    let mut update_acid = |pka: f64| {
        min_acid_pka = Some(min_acid_pka.map_or(pka, |cur| cur.min(pka)));
    };
    let mut update_base = |pka: f64| {
        max_base_pka = Some(max_base_pka.map_or(pka, |cur| cur.max(pka)));
    };

    for (idx, a) in mol.atoms.iter().enumerate() {
        // Permanent quaternary ammonium [N+] with 4 heavy bonds
        if a.atomic_number == 7 && a.formal_charge > 0 && mol.heavy_degree(idx) == 4 {
            permanent_charge += 1.0;
            continue;
        }

        // Acidic centers
        match a.atomic_number {
            8 => {
                // Carboxylic acid, sulfonic acid, phosphoric acid, phenol, hydroxamic acid
                let is_oh_or_om = a.total_h() > 0 || a.formal_charge < 0;
                if !is_oh_or_om {
                    continue;
                }
                for nbr in mol.neighbors(idx) {
                    let na = &mol.atoms[nbr];
                    if na.atomic_number == 6 {
                        let has_c_dbl_o = mol.adj[nbr].iter().any(|&(o2, bid)| {
                            o2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[o2].atomic_number == 8
                        });
                        if has_c_dbl_o {
                            // Check alpha-halogens or aromatic ring
                            let on_arom = mol
                                .neighbors(nbr)
                                .any(|x| x != idx && mol.atoms[x].aromatic);
                            let alpha_ewg = mol.neighbors(nbr).any(|x| {
                                x != idx
                                    && mol.atoms[x].atomic_number == 6
                                    && mol
                                        .neighbors(x)
                                        .any(|y| matches!(mol.atoms[y].atomic_number, 9 | 17 | 35))
                            });
                            let pka = if alpha_ewg {
                                2.8
                            } else if on_arom {
                                4.1
                            } else {
                                4.5
                            };
                            update_acid(pka);
                        } else if na.aromatic {
                            // Phenol (check ortho/para nitro or carbonyl)
                            let has_ewg = mol.neighbors(nbr).any(|x| {
                                x != idx
                                    && mol.neighbors(x).any(|y| {
                                        mol.atoms[y].atomic_number == 7
                                            && mol.atoms[y].formal_charge > 0
                                    })
                            });
                            update_acid(if has_ewg { 7.2 } else { 9.9 });
                        }
                    } else if na.atomic_number == 16 {
                        // Sulfonic acid R-SO3H
                        let dbl_o_cnt = mol.adj[nbr]
                            .iter()
                            .filter(|&&(o2, bid)| {
                                o2 != idx
                                    && mol.bonds[bid].order == BondOrder::Double
                                    && mol.atoms[o2].atomic_number == 8
                            })
                            .count();
                        if dbl_o_cnt >= 2 {
                            update_acid(1.0);
                        }
                    } else if na.atomic_number == 15 {
                        // Phosphate / phosphonate -P(=O)(OH)
                        let has_p_dbl_o = mol.adj[nbr].iter().any(|&(o2, bid)| {
                            o2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[o2].atomic_number == 8
                        });
                        if has_p_dbl_o {
                            update_acid(2.1);
                        }
                    }
                }
            }
            7 => {
                // Tetrazole acidic NH (5-ring with 4 nitrogens)
                if a.aromatic && a.total_h() > 0 && a.ring_sizes.contains(&5) {
                    let is_tetrazole = mol.rings.iter().any(|r| {
                        r.len() == 5
                            && r.contains(&idx)
                            && r.iter()
                                .filter(|&&rid| mol.atoms[rid].atomic_number == 7)
                                .count()
                                == 4
                    });
                    if is_tetrazole {
                        update_acid(4.9);
                        continue;
                    }
                }
                // Sulfonamide NH (weak acid ~9.8, acylsulfonamide ~4.8)
                let on_sulfonyl = mol.neighbors(idx).any(|nbr| {
                    mol.atoms[nbr].atomic_number == 16
                        && mol.adj[nbr]
                            .iter()
                            .filter(|&&(o, bid)| {
                                mol.bonds[bid].order == BondOrder::Double
                                    && mol.atoms[o].atomic_number == 8
                            })
                            .count()
                            >= 2
                });
                let on_carbonyl = mol.neighbors(idx).any(|nbr| {
                    mol.atoms[nbr].atomic_number == 6
                        && mol.adj[nbr].iter().any(|&(o, bid)| {
                            o != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[o].atomic_number == 8
                        })
                });
                if on_sulfonyl && a.total_h() > 0 {
                    update_acid(if on_carbonyl { 4.8 } else { 9.8 });
                    continue;
                }
                if on_carbonyl {
                    // Neutral amide / urea / lactam (not basic)
                    continue;
                }

                // Basic nitrogens
                // 1. Guanidine / Amidine C(=N)-N
                let on_amidine_c = mol.neighbors(idx).any(|c| {
                    mol.atoms[c].atomic_number == 6
                        && !mol.atoms[c].aromatic
                        && mol.adj[c].iter().any(|&(n2, bid)| {
                            n2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[n2].atomic_number == 7
                        })
                });
                if on_amidine_c {
                    update_base(12.2);
                    continue;
                }

                if a.aromatic {
                    // Imidazole (5-ring with 2 nitrogens, sp2 unprotonated N)
                    if a.total_h() == 0 && mol.heavy_degree(idx) == 2 {
                        let in_imidazole = mol.rings.iter().any(|r| {
                            r.len() == 5
                                && r.contains(&idx)
                                && r.iter()
                                    .filter(|&&rid| mol.atoms[rid].atomic_number == 7)
                                    .count()
                                    == 2
                        });
                        if in_imidazole {
                            update_base(6.9);
                        } else if a.ring_sizes.contains(&6) {
                            // Pyridine / pyrimidine
                            let n_in_ring = mol
                                .rings
                                .iter()
                                .filter(|r| r.len() == 6 && r.contains(&idx))
                                .map(|r| {
                                    r.iter()
                                        .filter(|&&rid| mol.atoms[rid].atomic_number == 7)
                                        .count()
                                })
                                .max()
                                .unwrap_or(1);
                            update_base(if n_in_ring >= 2 { 1.8 } else { 5.2 });
                        }
                    }
                    continue;
                }

                // Exclude nitrile, nitro, azide
                let has_multi = mol.adj[idx].iter().any(|&(_, bid)| {
                    matches!(mol.bonds[bid].order, BondOrder::Double | BondOrder::Triple)
                });
                if has_multi {
                    continue;
                }

                let adj_arom = mol.neighbors(idx).any(|n| mol.atoms[n].aromatic);
                if adj_arom {
                    update_base(4.6); // Aniline-like
                    continue;
                }

                // Aliphatic amine (primary / secondary / tertiary / piperazine / morpholine)
                let in_morpholine = mol.rings.iter().any(|r| {
                    r.len() == 6
                        && r.contains(&idx)
                        && r.iter().any(|&rid| mol.atoms[rid].atomic_number == 8)
                });
                let in_piperazine = mol.rings.iter().any(|r| {
                    r.len() == 6
                        && r.contains(&idx)
                        && r.iter()
                            .filter(|&&rid| mol.atoms[rid].atomic_number == 7)
                            .count()
                            == 2
                });
                let base_pka = if in_morpholine {
                    7.4
                } else if in_piperazine {
                    8.9
                } else if a.hybridization == Hybridization::Sp3 {
                    match mol.heavy_degree(idx) {
                        1 => 10.4,
                        2 => 10.1,
                        _ => 9.4,
                    }
                } else {
                    8.5
                };
                update_base(base_pka);
            }
            _ => {}
        }
    }

    // Henderson-Hasselbalch fractions ionized at pH 7.4
    let frac_deprot = min_acid_pka.map_or(0.0, |pka| 1.0 / (1.0 + 10.0_f64.powf(pka - PH)));
    let frac_prot = max_base_pka.map_or(0.0, |pka| 1.0 / (1.0 + 10.0_f64.powf(PH - pka)));
    let estimated_charge_ph74 = permanent_charge + frac_prot - frac_deprot;

    let is_acidic = min_acid_pka.is_some_and(|pka| pka < 6.5);
    let is_basic = max_base_pka.is_some_and(|pka| pka > 7.5) || permanent_charge > 0.5;
    let ionization_class = match (is_acidic, is_basic) {
        (true, true) => IonizationClass::Zwitterion,
        (true, false) => IonizationClass::Acid,
        (false, true) => IonizationClass::Base,
        (false, false) => IonizationClass::Neutral,
    };

    // Henderson-Hasselbalch LogD(7.4) with practical ion-pair floor (max drop capped at 3.5 log units)
    let mut delta_logd = 0.0_f64;
    if let Some(pka_a) = min_acid_pka {
        delta_logd += (1.0 + 10.0_f64.powf((PH - pka_a).min(4.0))).log10();
    }
    if let Some(pka_b) = max_base_pka {
        delta_logd += (1.0 + 10.0_f64.powf((pka_b - PH).min(4.0))).log10();
    }
    let logd_ph74 = lipo.consensus_log_p - delta_logd.min(3.5);

    IonizationProfile {
        strongest_acidic_pka: min_acid_pka.map(round2),
        strongest_basic_pka: max_base_pka.map(round2),
        ionization_class,
        estimated_charge_ph74: round2(estimated_charge_ph74),
        logd_ph74: round2(logd_ph74),
    }
}
