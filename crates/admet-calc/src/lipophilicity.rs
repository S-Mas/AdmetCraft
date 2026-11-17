use chem_core::{BondOrder, Hybridization, Molecule};

use crate::physicochem::{compute_tpsa_ertl, round2};
use crate::types::Lipophilicity;

/// Compute lipophilicity models (`WLOGP`, `XLOGP3`-heuristic, `MLOGP`), their 3-model
/// `consensus_log_p`, plus two experimental linear surrogates (`solvation_logp_est`,
/// `topological_logp_est`).
pub fn compute_lipophilicity(mol: &Molecule) -> Lipophilicity {
    let (wlogp, _mr) = compute_wlogp_and_mr(mol);
    let mlogp = compute_mlogp(mol);
    let xlogp3 = compute_xlogp3(mol);
    let solvation_logp_est = compute_solvation_logp_est(mol, xlogp3, wlogp);
    let topological_logp_est = compute_topological_logp_est(mol, xlogp3, wlogp);
    let consensus = (xlogp3 + wlogp + mlogp) / 3.0;

    Lipophilicity {
        wlogp: round2(wlogp),
        xlogp3: round2(xlogp3),
        mlogp: round2(mlogp),
        solvation_logp_est: round2(solvation_logp_est),
        topological_logp_est: round2(topological_logp_est),
        consensus_log_p: round2(consensus),
    }
}

/// Wildman & Crippen (1999) *J. Chem. Inf. Comput. Sci.* 39, 868–873
/// Atomic contributions for WLOGP and Molar Refractivity (MR).
pub fn compute_wlogp_and_mr(mol: &Molecule) -> (f64, f64) {
    let mut logp = 0.0;
    let mut mr = 0.0;

    for (idx, atom) in mol.atoms.iter().enumerate() {
        let z = atom.atomic_number;
        if z == 1 {
            // Standalone H atom
            let (lp, m) = crippen_h_contrib(mol, idx);
            logp += lp;
            mr += m;
            continue;
        }
        let (atom_lp, atom_mr) = crippen_heavy_contrib(mol, idx);
        logp += atom_lp;
        mr += atom_mr;

        // Add contributions for attached implicit + explicit hydrogens
        let h_cnt = f64::from(atom.total_h());
        if h_cnt > 0.0 {
            let (hlp, hmr) = crippen_attached_h_contrib(mol, idx);
            logp += h_cnt * hlp;
            mr += h_cnt * hmr;
        }
    }

    (logp, mr)
}

fn has_hetero_neighbor(mol: &Molecule, idx: usize) -> bool {
    mol.neighbors(idx)
        .any(|nbr| !matches!(mol.atoms[nbr].atomic_number, 1 | 6))
}

fn count_bonded_hetero(mol: &Molecule, idx: usize) -> usize {
    mol.neighbors(idx)
        .filter(|&nbr| !matches!(mol.atoms[nbr].atomic_number, 1 | 6))
        .count()
}

fn is_carbonyl_or_hetero_multi(mol: &Molecule, c_idx: usize) -> bool {
    mol.adj[c_idx].iter().any(|&(nbr, bid)| {
        let bo = mol.bonds[bid].order;
        let nz = mol.atoms[nbr].atomic_number;
        matches!(bo, BondOrder::Double | BondOrder::Triple) && !matches!(nz, 1 | 6)
    })
}

fn crippen_h_contrib(mol: &Molecule, h_idx: usize) -> (f64, f64) {
    if let Some(nbr) = mol.neighbors(h_idx).next() {
        crippen_attached_h_contrib(mol, nbr)
    } else {
        (0.1230, 1.057)
    }
}

fn crippen_attached_h_contrib(mol: &Molecule, parent_idx: usize) -> (f64, f64) {
    let p = &mol.atoms[parent_idx];
    match p.atomic_number {
        6 => {
            if is_carbonyl_or_hetero_multi(mol, parent_idx) || p.hybridization == Hybridization::Sp
            {
                (0.2667, 1.313) // H4: attached to sp C or C=X
            } else {
                (0.1230, 1.057) // H1: aliphatic or aromatic C-H
            }
        }
        7 => (0.2142, 0.9627), // H3: attached to N
        8 => {
            // H2: alcohol/phenol vs carboxylic acid
            let on_acid = mol.neighbors(parent_idx).any(|c| {
                mol.atoms[c].atomic_number == 6
                    && mol.adj[c].iter().any(|&(o2, bid)| {
                        o2 != parent_idx
                            && mol.bonds[bid].order == BondOrder::Double
                            && mol.atoms[o2].atomic_number == 8
                    })
            });
            if on_acid {
                (0.2980, 0.8939)
            } else {
                (-0.2677, 0.8939)
            }
        }
        16 => (0.1125, 1.057), // HS
        _ => (0.1230, 1.057),
    }
}

fn crippen_heavy_contrib(mol: &Molecule, idx: usize) -> (f64, f64) {
    let a = &mol.atoms[idx];
    let deg = mol.heavy_degree(idx);
    let hetero_nbrs = count_bonded_hetero(mol, idx);
    let mut double_c = false;
    let mut double_hetero = false;
    let mut triple_any = false;
    let mut arom_bonds = 0usize;
    for &(nbr, bid) in &mol.adj[idx] {
        let bo = mol.bonds[bid].order;
        let nz = mol.atoms[nbr].atomic_number;
        match bo {
            BondOrder::Double => {
                if nz == 6 {
                    double_c = true;
                } else {
                    double_hetero = true;
                }
            }
            BondOrder::Triple => triple_any = true,
            BondOrder::Aromatic => arom_bonds += 1,
            BondOrder::Single => {}
        }
    }

    match a.atomic_number {
        // Carbon types (C1..C27 from Wildman & Crippen 1999)
        6 => {
            if a.aromatic {
                let has_hetero_sub = mol.neighbors(idx).any(|nbr| {
                    let na = &mol.atoms[nbr];
                    !matches!(na.atomic_number, 1 | 6)
                        && mol
                            .bond_between(idx, nbr)
                            .is_some_and(|b| b.order != BondOrder::Aromatic)
                });
                let has_arom_hetero_nbr = mol.neighbors(idx).any(|nbr| {
                    let na = &mol.atoms[nbr];
                    !matches!(na.atomic_number, 1 | 6)
                        && mol
                            .bond_between(idx, nbr)
                            .is_some_and(|b| b.order == BondOrder::Aromatic)
                });
                if double_hetero {
                    (0.0566, 3.350) // aromatic carbonyl C (e.g. caffeine, coumarin)
                } else if a.total_h() > 0 {
                    (0.1581, 3.350) // C18: aromatic CH
                } else if has_hetero_sub {
                    (0.1360, 3.743) // C20/C21: aromatic C attached to hetero substituent
                } else if has_arom_hetero_nbr {
                    (-0.0044, 3.465) // C19: aromatic C adjacent to aromatic hetero
                } else {
                    (0.2961, 3.509) // C21: aromatic C-C (substituted / bridgehead)
                }
            } else if triple_any || a.hybridization == Hybridization::Sp {
                if hetero_nbrs > 0 {
                    (0.0017, 3.413) // C10: sp C with hetero
                } else {
                    (0.3300, 3.875) // C8: sp C
                }
            } else if double_hetero {
                // C=O, C=N, C=S
                if hetero_nbrs >= 2 {
                    (-0.0244, 2.753) // ester/acid/amide/urea carbonyl C
                } else {
                    (-0.1830, 3.120) // ketone/aldehyde C=O
                }
            } else if double_c || a.hybridization == Hybridization::Sp2 {
                if hetero_nbrs > 0 {
                    (0.0017, 3.413) // C7/C10: =C-X
                } else if a.total_h() > 0 {
                    (0.1551, 3.509) // C5/C6: =CH- / =CH2
                } else {
                    (0.2640, 3.882) // C6: =C<
                }
            } else {
                // sp3 carbon
                if hetero_nbrs == 0 {
                    match deg {
                        0 | 1 => (0.1441, 2.503), // C1: primary CH3
                        2 => (0.1441, 2.503),     // C1: secondary CH2
                        3 => (-0.2051, 2.433),    // C2: tertiary CH
                        _ => (-0.2051, 2.433),    // C2: quaternary C
                    }
                } else {
                    // sp3 C attached to heteroatom (C3/C4)
                    let attached_to_f_or_cl = mol
                        .neighbors(idx)
                        .any(|n| matches!(mol.atoms[n].atomic_number, 9 | 17 | 35 | 53));
                    if attached_to_f_or_cl && hetero_nbrs >= 2 {
                        (0.2050, 2.680) // CX2 / CX3 (e.g. CF3)
                    } else if deg <= 2 {
                        (-0.2035, 2.753) // C3: CH3-X, CH2-X
                    } else {
                        (-0.2051, 2.433) // C4: CH-X, C-X
                    }
                }
            }
        }
        // Nitrogen types (N1..N14)
        7 => {
            if a.formal_charge > 0 {
                if double_hetero {
                    (-0.3396, 2.800) // nitro N+
                } else {
                    (-0.6600, 2.100) // quaternary N+
                }
            } else if a.aromatic {
                (-0.3239, 2.202) // N11/N12: aromatic nitrogen
            } else if triple_any {
                (-0.5666, 2.936) // N9: nitrile #N
            } else if double_c || double_hetero {
                if a.total_h() > 0 {
                    (-0.4806, 2.500) // imine =NH
                } else {
                    (-0.3187, 3.000) // imine =N-
                }
            } else {
                // sp3 amine or amide N
                let is_amide = mol.neighbors(idx).any(|c| {
                    matches!(mol.atoms[c].atomic_number, 6 | 16)
                        && mol.adj[c].iter().any(|&(o2, bid)| {
                            o2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[o2].atomic_number == 8
                        })
                });
                let is_aromatic_adj = mol.neighbors(idx).any(|nbr| mol.atoms[nbr].aromatic);
                if is_amide {
                    if a.total_h() >= 1 {
                        (-0.5060, 2.400) // amide NH
                    } else {
                        (-0.4458, 2.830) // tertiary amide N
                    }
                } else if is_aromatic_adj {
                    (-0.4458, 3.000) // aniline N
                } else {
                    match deg {
                        0 | 1 => (-1.0190, 2.262), // N1: R-NH2
                        2 => (-0.7096, 2.173),     // N2: R2NH
                        _ => (-0.3187, 2.836),     // N7: R3N
                    }
                }
            }
        }
        // Oxygen types (O1..O12)
        8 => {
            if a.formal_charge < 0 {
                (-0.7790, 1.380) // O-
            } else if a.aromatic {
                (0.1552, 1.080) // O1: furan o
            } else if double_c || double_hetero {
                (-0.2476, 0.790) // O9: carbonyl/sulfoxide/sulfone =O
            } else if a.total_h() > 0 {
                let on_arom = mol.neighbors(idx).any(|n| mol.atoms[n].aromatic);
                if on_arom {
                    (-0.2893, 1.150) // O2: phenolic -OH
                } else {
                    (-0.2893, 1.080) // O2: alcoholic/carboxylic -OH
                }
            } else {
                // -O- ether / ester
                let on_arom = mol.neighbors(idx).any(|n| mol.atoms[n].aromatic);
                let on_carbonyl = mol.neighbors(idx).any(|c| {
                    mol.atoms[c].atomic_number == 6
                        && mol.adj[c].iter().any(|&(o2, bid)| {
                            o2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[o2].atomic_number == 8
                        })
                });
                if on_arom && on_carbonyl {
                    (-0.0850, 1.180) // aryl ester -O-
                } else if on_arom {
                    (-0.1414, 1.180) // O4: Ar-O-R
                } else if on_carbonyl {
                    (-0.2000, 1.080) // alkyl ester -O-
                } else {
                    (-0.4195, 1.080) // O3: R-O-R
                }
            }
        }
        // Halogens
        9 => (0.4202, 1.108),  // F
        17 => (0.6895, 5.853), // Cl
        35 => (0.8456, 8.927), // Br
        53 => (0.8857, 14.02), // I
        // Sulfur
        16 => {
            if arom_bonds >= 2 {
                (0.6237, 7.591) // aromatic s
            } else if double_hetero {
                (-0.0024, 6.920) // S(=O) or S(=O)2
            } else {
                (0.6482, 7.591) // thioether / thiol
            }
        }
        // Phosphorus
        15 => (0.8612, 6.920),
        _ => (0.0, 2.0),
    }
}

/// Moriguchi et al. (1992, 1994) *Chem. Pharm. Bull.* 13-parameter MLOGP model.
pub fn compute_mlogp(mol: &Molecule) -> f64 {
    let mut cx = 0.0f64;
    let mut no = 0.0f64;
    let mut ub = 0.0f64;
    let mut has_qn = false;
    let mut no2_count = 0.0f64;
    let mut halogens = 0usize;

    for (idx, a) in mol.atoms.iter().enumerate() {
        match a.atomic_number {
            6 => cx += 1.0,
            7 | 8 => no += 1.0,
            9 => {
                cx += 0.5;
                halogens += 1;
            }
            17 => {
                cx += 1.0;
                halogens += 1;
            }
            35 => {
                cx += 1.5;
                halogens += 1;
            }
            53 => {
                cx += 2.0;
                halogens += 1;
            }
            16 => cx += 1.0,
            15 => cx += 1.0,
            _ => {}
        }
        if a.atomic_number == 7 && a.formal_charge > 0 && mol.degree(idx) == 4 {
            has_qn = true;
        }
        // Nitro group check: N bonded to 2 oxygens with at least one =O
        if a.atomic_number == 7 && mol.degree(idx) == 3 {
            let o_nbrs: Vec<usize> = mol
                .neighbors(idx)
                .filter(|&n| mol.atoms[n].atomic_number == 8 && mol.degree(n) == 1)
                .collect();
            if o_nbrs.len() == 2 {
                no2_count += 1.0;
            }
        }
    }

    // Unsaturated bonds (UB): double=1, triple=2, aromatic=0.5 per bond
    for b in &mol.bonds {
        match b.order {
            BondOrder::Double => ub += 1.0,
            BondOrder::Triple => ub += 2.0,
            BondOrder::Aromatic => ub += 0.5,
            BondOrder::Single => {}
        }
    }

    // Proximity effect (PRX): pairs of N/O separated by 1 or 2 bonds
    let mut prx = 0.0f64;
    for (i, a) in mol.atoms.iter().enumerate() {
        if !matches!(a.atomic_number, 7 | 8) {
            continue;
        }
        for (j, b) in mol.atoms.iter().enumerate().skip(i + 1) {
            if !matches!(b.atomic_number, 7 | 8) {
                continue;
            }
            if mol.bond_between(i, j).is_some() {
                prx += 2.0;
            } else {
                // Check 2-bond separation (e.g. O=C-O or O=C-N or N-C=N)
                let dist2 = mol
                    .neighbors(i)
                    .any(|mid| mol.bond_between(mid, j).is_some());
                if dist2 {
                    prx += 1.0;
                }
            }
        }
    }
    prx = prx.min(6.0);

    // POL: number of aromatic polar substituents (e.g. OH/NH2/COOH attached to aromatic ring)
    let mut pol = 0.0f64;
    for (idx, a) in mol.atoms.iter().enumerate() {
        if a.aromatic && has_hetero_neighbor(mol, idx) {
            pol += 0.5;
        }
    }
    pol = pol.min(4.0);

    // RNG: ring indicator (1 if non-aromatic ring / >1 ring, Moriguchi RNG = 0 or 1 for alicyclic)
    let has_non_arom_ring = mol
        .rings
        .iter()
        .any(|r| r.iter().any(|&aid| !mol.atoms[aid].aromatic));
    let rng = if has_non_arom_ring { 1.0 } else { 0.0 };

    // ALK: alkane / hydrocarbon flag (only C, H, and optional halogens)
    let is_alk = no == 0.0 && halogens == 0;
    let alk = if is_alk { 1.0 } else { 0.0 };

    // AMP: amphoteric (has both basic amine and acidic COOH/SO3H)
    let has_acid = mol.atoms.iter().enumerate().any(|(i, a)| {
        a.atomic_number == 6
            && mol.adj[i].iter().any(|&(o1, b1)| {
                mol.bonds[b1].order == BondOrder::Double && mol.atoms[o1].atomic_number == 8
            })
            && mol.adj[i].iter().any(|&(o2, b2)| {
                mol.bonds[b2].order == BondOrder::Single
                    && mol.atoms[o2].atomic_number == 8
                    && mol.atoms[o2].total_h() > 0
            })
    });
    let has_basic_n = mol.atoms.iter().enumerate().any(|(i, a)| {
        a.atomic_number == 7 && !a.aromatic && a.formal_charge == 0 && !is_amide_n(mol, i)
    });
    let amp = if has_acid && has_basic_n { 1.0 } else { 0.0 };

    // HB: intramolecular H-bond indicator (e.g. salicylic acid / beta-diketone / ortho-substituted)
    let hb = if has_intramol_hbond(mol) { 1.0 } else { 0.0 };

    // QN, NCS, BLM
    let qn = if has_qn { 1.0 } else { 0.0 };
    let ncs = 0.0;
    let blm = if mol
        .rings
        .iter()
        .any(|r| r.len() == 4 && r.iter().any(|&a| mol.atoms[a].atomic_number == 7))
    {
        1.0
    } else {
        0.0
    };

    1.244 * cx.powf(0.6) - 1.017 * no.powf(0.9) + 0.406 * prx - 0.145 * ub.powf(0.8)
        + 0.511 * hb
        + 0.268 * pol
        - 2.215 * amp
        + 0.912 * alk
        - 0.392 * rng
        - 3.684 * qn
        + 0.474 * no2_count
        + 1.582 * ncs
        + 0.773 * blm
        - 1.041
}

fn is_amide_n(mol: &Molecule, n_idx: usize) -> bool {
    mol.neighbors(n_idx).any(|c| {
        matches!(mol.atoms[c].atomic_number, 6 | 16)
            && mol.adj[c].iter().any(|&(o, b)| {
                o != n_idx
                    && mol.bonds[b].order == BondOrder::Double
                    && mol.atoms[o].atomic_number == 8
            })
    })
}

fn has_intramol_hbond(mol: &Molecule) -> bool {
    // Detect 5- or 6-membered intramolecular H-bond donor-acceptor pair (e.g., salicylate, enol-ketone)
    for (i, a) in mol.atoms.iter().enumerate() {
        if !matches!(a.atomic_number, 7 | 8) || a.total_h() == 0 {
            continue;
        }
        for n1 in mol.neighbors(i) {
            for n2 in mol.neighbors(n1) {
                if n2 == i {
                    continue;
                }
                for n3 in mol.neighbors(n2) {
                    if n3 == n1 {
                        continue;
                    }
                    // 1,4-donor-acceptor (5-ring H-bond) or 1,5-donor-acceptor (6-ring H-bond)
                    if matches!(mol.atoms[n3].atomic_number, 7 | 8) {
                        return true;
                    }
                    for n4 in mol.neighbors(n3) {
                        if n4 != n2 && matches!(mol.atoms[n4].atomic_number, 7 | 8) {
                            return true;
                        }
                    }
                }
            }
        }
    }
    false
}

/// Cheng et al. (2007) *J. Chem. Inf. Model.* XLOGP3 atom-additive + correction-factor model.
pub fn compute_xlogp3(mol: &Molecule) -> f64 {
    let mut val = 0.0;
    let mut hydrophobic_c_count = 0usize;

    for (idx, a) in mol.atoms.iter().enumerate() {
        let deg = mol.heavy_degree(idx);
        let h = f64::from(a.total_h());
        let hetero_nbrs = count_bonded_hetero(mol, idx);
        match a.atomic_number {
            6 => {
                if a.aromatic {
                    if hetero_nbrs == 0 {
                        if h > 0.0 {
                            val += 0.337; // aromatic CH
                            hydrophobic_c_count += 1;
                        } else {
                            val += 0.296; // aromatic C
                            hydrophobic_c_count += 1;
                        }
                    } else {
                        val += 0.130; // aromatic C-X
                    }
                } else if a.hybridization == Hybridization::Sp {
                    val += if h > 0.0 { 0.25 } else { 0.20 };
                } else if a.hybridization == Hybridization::Sp2 {
                    if is_carbonyl_or_hetero_multi(mol, idx) {
                        if hetero_nbrs >= 2 {
                            val -= 0.04; // COOH / ester / amide C=O
                        } else {
                            val += 0.02; // ketone / aldehyde
                        }
                    } else {
                        val += if h > 0.0 { 0.28 } else { 0.35 };
                    }
                } else {
                    // sp3 C
                    if hetero_nbrs == 0 {
                        hydrophobic_c_count += 1;
                        val += match deg {
                            0 | 1 => 0.528, // CH3
                            2 => 0.358,     // CH2
                            3 => 0.120,     // CH
                            _ => -0.050,    // C
                        };
                    } else {
                        let hal_cnt = mol
                            .neighbors(idx)
                            .filter(|&n| matches!(mol.atoms[n].atomic_number, 9 | 17 | 35 | 53))
                            .count();
                        if hal_cnt >= 2 {
                            val += 0.35; // -CF3 / -CCl3 carbon
                        } else {
                            val += match deg {
                                0 | 1 => 0.220,
                                2 => -0.010,
                                _ => -0.220,
                            };
                        }
                    }
                }
            }
            7 => {
                if a.aromatic {
                    if a.total_h() > 0 || deg >= 3 {
                        val -= 0.250; // pyrrole-like
                    } else {
                        val -= 0.493; // pyridine-like
                    }
                } else if is_amide_n(mol, idx) {
                    val -= if a.total_h() > 0 { 0.650 } else { 0.420 };
                } else if a.formal_charge > 0 {
                    val -= 1.25;
                } else if mol.neighbors(idx).any(|n| mol.atoms[n].aromatic) {
                    val -= 0.350; // aniline-like
                } else {
                    val -= match deg {
                        0 | 1 => 0.890, // primary amine
                        2 => 0.650,     // secondary amine
                        _ => 0.380,     // tertiary amine
                    };
                }
            }
            8 => {
                if a.aromatic {
                    val -= 0.10; // furan
                } else if mol.adj[idx]
                    .iter()
                    .any(|&(_, b)| mol.bonds[b].order == BondOrder::Double)
                {
                    val -= 0.399; // =O
                } else if a.total_h() > 0 {
                    let on_acid = mol.neighbors(idx).any(|c| {
                        mol.atoms[c].atomic_number == 6
                            && mol.adj[c].iter().any(|&(o2, bid)| {
                                o2 != idx
                                    && mol.bonds[bid].order == BondOrder::Double
                                    && mol.atoms[o2].atomic_number == 8
                            })
                    });
                    let on_arom = mol.neighbors(idx).any(|n| mol.atoms[n].aromatic);
                    if on_acid {
                        val -= 0.32;
                    } else if on_arom {
                        val -= 0.04; // phenol OH
                    } else {
                        val -= 0.467; // alcohol OH
                    }
                } else {
                    val -= 0.180; // ether / ester -O-
                }
            }
            9 => val += 0.38,
            17 => val += 0.72,
            35 => val += 0.88,
            53 => val += 1.05,
            16 => {
                if a.aromatic {
                    val += 0.55;
                } else if mol.adj[idx]
                    .iter()
                    .any(|&(_, b)| mol.bonds[b].order == BondOrder::Double)
                {
                    val -= 0.95; // sulfoxide/sulfone S
                } else {
                    val += 0.42; // thioether
                }
            }
            15 => val += 0.30,
            _ => {}
        }
    }

    // Hydrophobic carbon pair correction (characteristic of XLOGP3)
    if hydrophobic_c_count >= 10 {
        val += 0.05 * ((hydrophobic_c_count - 9) as f64).min(8.0);
    }
    // Intramolecular H-bond bonus
    if has_intramol_hbond(mol) {
        val += 0.25;
    }

    val
}

/// Experimental solvation-inspired linear logP surrogate (NOT SwissADME iLOGP).
fn compute_solvation_logp_est(mol: &Molecule, xlogp3: f64, wlogp: f64) -> f64 {
    let tpsa = compute_tpsa_ertl(mol);
    let heavy = mol.heavy_atom_count() as f64;
    let arom = mol.aromatic_heavy_atom_count() as f64;
    let nonpolar_term = 0.145 * heavy + 0.055 * arom;
    let polar_desolv = 0.021 * tpsa;
    let phys_est = nonpolar_term - polar_desolv + 0.18;
    0.50 * phys_est + 0.25 * xlogp3 + 0.25 * wlogp
}

/// Experimental topological linear logP surrogate (NOT Silicos-IT).
fn compute_topological_logp_est(mol: &Molecule, xlogp3: f64, wlogp: f64) -> f64 {
    let tpsa = compute_tpsa_ertl(mol);
    let heavy = mol.heavy_atom_count() as f64;
    let n_hetero = mol
        .atoms
        .iter()
        .filter(|a| !matches!(a.atomic_number, 1 | 6))
        .count() as f64;
    let topo = 0.16 * (heavy - n_hetero) - 0.025 * tpsa + 0.10;
    0.35 * topo + 0.35 * xlogp3 + 0.30 * wlogp
}
