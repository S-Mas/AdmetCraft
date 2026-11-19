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
        6 => (0.1230, 1.0570),       // H1: aliphatic/aromatic/sp C-H
        7 => (0.2142, 0.9627),       // H3: N-H
        8 | 16 => (-0.2677, 1.3950), // H2: O-H / S-H (Wildman & Crippen 1999)
        _ => (0.1230, 1.0570),
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
    let mut arom_nbr_cnt = 0usize;
    for &(nbr, bid) in &mol.adj[idx] {
        if mol.atoms[nbr].atomic_number == 1 {
            continue;
        }
        if mol.atoms[nbr].aromatic {
            arom_nbr_cnt += 1;
        }
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
        // Carbon types (Wildman & Crippen 1999 C1..C27)
        6 => {
            if a.aromatic {
                let exo_nbr = mol.adj[idx].iter().copied().find(|&(nbr, bid)| {
                    mol.atoms[nbr].atomic_number != 1 && mol.bonds[bid].order != BondOrder::Aromatic
                });
                if let Some((ev, bid)) = exo_nbr {
                    if mol.bonds[bid].order == BondOrder::Double
                        && matches!(mol.atoms[ev].atomic_number, 8 | 16)
                    {
                        // Aromatic carbonyl c(=O) (caffeine, uracil, coumarin, pyridone)
                        (-0.8186, 4.0000)
                    } else {
                        match mol.atoms[ev].atomic_number {
                            6 => {
                                if mol.atoms[ev].aromatic {
                                    (0.2713, 3.9040)
                                } else {
                                    (0.1360, 3.5090)
                                }
                            }
                            8 => (0.5437, 3.8530),
                            7 => (0.4619, 4.0670),
                            16 => (0.2450, 3.5640),
                            9 => (0.0000, 3.2570),
                            17 => (0.2450, 3.5640),
                            35 => (0.1980, 3.1800),
                            53 => (0.0000, 3.1040),
                            _ => (0.1581, 3.3500),
                        }
                    }
                } else if deg >= 3 {
                    // Fused aromatic junction c(a)(a)a
                    if hetero_nbrs > 0 {
                        (0.2955, 3.9040)
                    } else {
                        (0.2713, 3.9040)
                    }
                } else {
                    // Unsubstituted aromatic cH (C18)
                    (0.1581, 3.3500)
                }
            } else if triple_any || a.hybridization == Hybridization::Sp {
                // sp alkyne / nitrile C
                (0.0017, 3.8880)
            } else if double_hetero {
                // Carbonyl / imine / thiocarbonyl sp2 C
                (-0.2783, 5.0070)
            } else if double_c || a.hybridization == Hybridization::Sp2 {
                if hetero_nbrs > 0 {
                    (-0.0300, 3.6000)
                } else {
                    (0.1551, 3.5130)
                }
            } else {
                // sp3 aliphatic C
                let hal_cnt = mol
                    .neighbors(idx)
                    .filter(|&n| matches!(mol.atoms[n].atomic_number, 9 | 17 | 35 | 53))
                    .count();
                if hal_cnt >= 2 {
                    (-0.0967, 2.5760)
                } else if arom_nbr_cnt >= 1 {
                    if deg == 1 {
                        (0.0845, 2.5030)
                    } else {
                        (-0.0516, 2.4880)
                    }
                } else if hetero_nbrs >= 2 {
                    (-0.2051, 2.7530)
                } else if hetero_nbrs == 1 {
                    let on_arom_het = mol.neighbors(idx).any(|n| {
                        !matches!(mol.atoms[n].atomic_number, 1 | 6) && mol.atoms[n].aromatic
                    });
                    if on_arom_het {
                        (-0.1444, 2.4120)
                    } else {
                        (-0.2035, 2.7530)
                    }
                } else if deg >= 3 {
                    // C2: tertiary CH / quaternary C
                    (0.0000, 2.4330)
                } else {
                    // C1: CH3 / CH2
                    (0.1441, 2.5030)
                }
            }
        }
        // Nitrogen types (N1..N14)
        7 => {
            if a.aromatic {
                (-0.3239, 2.2020)
            } else if a.formal_charge > 0 {
                let o_cnt = mol
                    .neighbors(idx)
                    .filter(|&n| mol.atoms[n].atomic_number == 8)
                    .count();
                let dbl_n_cnt = mol.adj[idx]
                    .iter()
                    .filter(|&&(n, bid)| {
                        mol.atoms[n].atomic_number == 7 && mol.bonds[bid].order == BondOrder::Double
                    })
                    .count();
                if o_cnt >= 2 {
                    (-0.3396, 0.2604) // nitro [N+]
                } else if dbl_n_cnt >= 2 {
                    (0.2887, 3.3590) // azide =[N+]=
                } else {
                    (-0.3396, 1.8000)
                }
            } else if a.formal_charge < 0 {
                (0.2887, 3.3590) // azide =[N-]
            } else if triple_any {
                (0.0151, 1.7250) // nitrile #N
            } else if double_c || double_hetero {
                let is_guanidine = mol.neighbors(idx).any(|c| {
                    mol.atoms[c].atomic_number == 6
                        && mol
                            .neighbors(c)
                            .filter(|&x| mol.atoms[x].atomic_number == 7)
                            .count()
                            >= 3
                });
                if is_guanidine {
                    (-0.3239, 2.2020)
                } else {
                    (0.1836, 2.4280) // imine / azo / azide N1
                }
            } else {
                let is_amide = is_amide_n(mol, idx);
                let is_arom_adj = mol.neighbors(idx).any(|nbr| mol.atoms[nbr].aromatic);
                if is_arom_adj {
                    match deg {
                        0 | 1 => (-1.0270, 2.8270),
                        2 => (-0.5188, 2.9500),
                        _ => (-0.4458, 2.8190),
                    }
                } else if is_amide {
                    match deg {
                        0 | 1 => (-1.0190, 2.2620),
                        #[allow(clippy::approx_constant)]
                        2 => (-0.5262, 2.7180),
                        _ => (-0.3187, 2.8190),
                    }
                } else {
                    match deg {
                        0 | 1 => (-1.0190, 2.2620),
                        2 => (-0.7096, 2.1730),
                        _ => (-0.3187, 2.8190),
                    }
                }
            }
        }
        // Oxygen types (O1..O12)
        8 => {
            if a.aromatic {
                (0.1552, 1.0800)
            } else if a.formal_charge < 0 {
                let adj_n_plus = mol
                    .neighbors(idx)
                    .any(|n| mol.atoms[n].atomic_number == 7 && mol.atoms[n].formal_charge > 0);
                if adj_n_plus {
                    (0.0335, 3.3670) // nitro [O-]
                } else {
                    (-0.3339, 0.7774)
                }
            } else if double_c || double_hetero {
                let adj_n_plus = mol
                    .neighbors(idx)
                    .any(|n| mol.atoms[n].atomic_number == 7 && mol.atoms[n].formal_charge > 0);
                let adj_s = mol.neighbors(idx).any(|n| mol.atoms[n].atomic_number == 16);
                let adj_p = mol.neighbors(idx).any(|n| mol.atoms[n].atomic_number == 15);
                let adj_arom_c = mol
                    .neighbors(idx)
                    .any(|n| mol.atoms[n].atomic_number == 6 && mol.atoms[n].aromatic);
                let adj_amide_or_acid = mol.neighbors(idx).any(|n| {
                    mol.atoms[n].atomic_number == 6
                        && mol
                            .neighbors(n)
                            .any(|u| u != idx && matches!(mol.atoms[u].atomic_number, 7 | 8))
                });
                if adj_n_plus {
                    (0.0335, 3.3670)
                } else if adj_s {
                    (-0.3339, 0.7774)
                } else if adj_p {
                    (-0.1188, 0.6865)
                } else if adj_arom_c {
                    (0.1788, 1.0000)
                } else if adj_amide_or_acid {
                    (0.1129, 0.7774)
                } else {
                    (-0.1526, 0.0000)
                }
            } else if a.total_h() > 0 {
                (-0.2893, 0.8238) // O2: alcohol / phenol / carboxylic OH / phosphate OH
            } else {
                let on_carbonyl = mol.neighbors(idx).any(|c| {
                    mol.atoms[c].atomic_number == 6
                        && mol.adj[c].iter().any(|&(o2, bid)| {
                            o2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[o2].atomic_number == 8
                        })
                });
                if on_carbonyl {
                    (-0.0684, 1.0850) // ester/acid -O-C(=O)
                } else {
                    (-0.4195, 1.1820) // aliphatic or aryl ether -O-
                }
            }
        }
        // Halogens
        9 => (0.4202, 1.1080),
        17 => (0.6895, 5.8530),
        35 => (0.8456, 8.9270),
        53 => (0.8857, 14.0200),
        // Sulfur
        16 => {
            if arom_bonds >= 2 {
                (0.6237, 6.6910)
            } else if double_hetero || deg >= 3 {
                (-0.0024, 7.3650)
            } else {
                (0.6482, 7.5910)
            }
        }
        // Phosphorus
        15 => (0.8612, 6.9200),
        5 => (-0.3808, 4.0410),
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
