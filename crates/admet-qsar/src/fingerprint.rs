#![forbid(unsafe_code)]

use admet_calc::{HergRisk, MoleculeEvaluation};
use chem_core::{BondOrder, ChiralTag, Hybridization, Molecule};
use serde::{Deserialize, Serialize};

/// Number of dense physicochemical + graph descriptors prepended before ECFP4 bits.
pub const NUM_DENSE_FEATURES: usize = 56;
/// Number of folded Morgan / ECFP4 fingerprint bins (1024 = 16 x u64 words).
pub const FP_BITS: usize = 1024;
/// Total feature vector dimension (`NUM_DENSE_FEATURES + FP_BITS`).
pub const TOTAL_FEATURES: usize = NUM_DENSE_FEATURES + FP_BITS;

/// Compact 1024-bit ECFP4 (Morgan radius 2) fingerprint with hardware POPCNT Tanimoto.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct Ecfp4BitVec {
    pub words: [u64; 16],
}

impl Ecfp4BitVec {
    #[inline]
    pub fn set_bit(&mut self, bit: usize) {
        let idx = (bit & 1023) >> 6;
        let off = bit & 63;
        self.words[idx] |= 1u64 << off;
    }

    #[inline]
    #[must_use]
    pub fn popcount(&self) -> u32 {
        let mut c = 0u32;
        for w in self.words {
            c += w.count_ones();
        }
        c
    }

    /// Exact Tanimoto (Jaccard) similarity in `[0.0, 1.0]` using `u64::count_ones()`.
    #[inline]
    #[must_use]
    pub fn tanimoto(&self, other: &Self) -> f64 {
        let mut inter = 0u32;
        let mut uni = 0u32;
        for i in 0..16 {
            inter += (self.words[i] & other.words[i]).count_ones();
            uni += (self.words[i] | other.words[i]).count_ones();
        }
        if uni == 0 {
            1.0
        } else {
            f64::from(inter) / f64::from(uni)
        }
    }
}

/// Complete featurization of a molecule for QSAR training/inference and AD checking.
#[derive(Debug, Clone)]
pub struct MolFeatures {
    pub dense: [f32; NUM_DENSE_FEATURES],
    pub fp: Ecfp4BitVec,
    pub fp_counts: [u8; FP_BITS],
    pub scaffold_hash: u64,
}

impl MolFeatures {
    /// Materialize the full 1,080-element feature vector (`56 dense + 1024 ECFP4 counts`).
    #[must_use]
    pub fn to_vec(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(TOTAL_FEATURES);
        out.extend_from_slice(&self.dense);
        for &c in &self.fp_counts {
            out.push(f32::from(c));
        }
        out
    }

    /// Read feature index `idx` directly without allocating.
    #[inline]
    #[must_use]
    pub fn get(&self, idx: usize) -> f32 {
        if idx < NUM_DENSE_FEATURES {
            self.dense[idx]
        } else if idx < TOTAL_FEATURES {
            f32::from(self.fp_counts[idx - NUM_DENSE_FEATURES])
        } else {
            0.0
        }
    }
}

#[inline]
fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

#[inline]
fn bond_code(order: BondOrder) -> u64 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

#[inline]
fn hyb_code(h: Hybridization) -> u64 {
    match h {
        Hybridization::Sp => 1,
        Hybridization::Sp2 => 2,
        Hybridization::Sp3 => 3,
        Hybridization::Other => 0,
    }
}

/// Computes circular Morgan / ECFP4 (radius 2) bitvector and capped-count array (`1024` bins).
#[must_use]
pub fn compute_ecfp4(mol: &Molecule) -> (Ecfp4BitVec, [u8; FP_BITS]) {
    let n = mol.atoms.len();
    let mut fp = Ecfp4BitVec::default();
    let mut counts = [0u8; FP_BITS];
    if n == 0 {
        return (fp, counts);
    }

    let mut inv = vec![0u64; n];
    let mut heavy_mask = vec![false; n];
    for (i, a) in mol.atoms.iter().enumerate() {
        if a.atomic_number <= 1 {
            continue;
        }
        heavy_mask[i] = true;
        let heavy_deg = mol.adj[i]
            .iter()
            .filter(|&&(nb, _)| mol.atoms.get(nb).is_some_and(|x| x.atomic_number > 1))
            .count() as u64;
        let min_ring = a.ring_sizes.iter().copied().min().unwrap_or(0) as u64;
        let mut h = u64::from(a.atomic_number);
        h = mix64(h ^ (heavy_deg << 8));
        h = mix64(h ^ (u64::from(a.total_h()) << 16));
        h = mix64(h ^ ((i64::from(a.formal_charge) + 8) as u64) << 24);
        h = mix64(h ^ (u64::from(a.aromatic) << 32) ^ (u64::from(a.in_ring) << 33));
        h = mix64(h ^ (min_ring << 36) ^ (hyb_code(a.hybridization) << 44));
        inv[i] = h;
        let bit = (h as usize) & (FP_BITS - 1);
        fp.set_bit(bit);
        counts[bit] = counts[bit].saturating_add(1).min(8);
    }

    // Radius 1 and Radius 2 Morgan circular updates
    let mut next_inv = vec![0u64; n];
    let mut nbr_pairs: Vec<(u64, u64)> = Vec::with_capacity(8);
    for radius in 1u64..=2u64 {
        for i in 0..n {
            if !heavy_mask[i] {
                continue;
            }
            nbr_pairs.clear();
            for &(nb, bid) in &mol.adj[i] {
                if !heavy_mask[nb] {
                    continue;
                }
                let bc = mol.bonds.get(bid).map_or(1, |b| bond_code(b.order));
                nbr_pairs.push((bc, inv[nb]));
            }
            nbr_pairs.sort_unstable();
            let mut h = mix64(inv[i] ^ (radius.wrapping_mul(0x517c_c1b7_2722_0a95)));
            for &(bc, nh) in &nbr_pairs {
                h = mix64(h ^ mix64(nh.wrapping_add(bc.wrapping_mul(0x9e37_79b9))));
            }
            next_inv[i] = h;
            let bit = (h as usize) & (FP_BITS - 1);
            fp.set_bit(bit);
            counts[bit] = counts[bit].saturating_add(1).min(8);
        }
        inv.copy_from_slice(&next_inv);
    }

    (fp, counts)
}

/// Computes a deterministic Bemis-Murcko (1996) molecular scaffold hash for zero-leakage
/// scaffold splitting (`train / cal / test`).
#[must_use]
pub fn compute_murcko_scaffold_hash(mol: &Molecule) -> u64 {
    let n = mol.atoms.len();
    if n == 0 {
        return 0;
    }
    let mut active: Vec<bool> = mol.atoms.iter().map(|a| a.atomic_number > 1).collect();

    // Iteratively prune non-ring degree-1 sidechains, preserving ring atoms, inter-ring linkers,
    // and exocyclic double-bonded heteroatoms directly attached to a ring/linker.
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..n {
            if !active[i] || mol.atoms[i].in_ring {
                continue;
            }
            let mut active_nbrs = 0usize;
            let mut single_nbr_bond: Option<(usize, BondOrder)> = None;
            for &(nb, bid) in &mol.adj[i] {
                if active[nb] {
                    active_nbrs += 1;
                    if let Some(b) = mol.bonds.get(bid) {
                        single_nbr_bond = Some((nb, b.order));
                    }
                }
            }
            if active_nbrs == 0 {
                active[i] = false;
                changed = true;
            } else if active_nbrs == 1 {
                let keep_exo_double = match single_nbr_bond {
                    Some((nb, BondOrder::Double)) => {
                        mol.atoms[i].atomic_number != 6 && mol.atoms[nb].in_ring
                    }
                    _ => false,
                };
                if !keep_exo_double {
                    active[i] = false;
                    changed = true;
                }
            }
        }
    }

    // If molecule has no rings (acyclic), hash the full heavy-atom graph so identical acyclic
    // frameworks still group together deterministically.
    if !active.iter().any(|&x| x) {
        for (i, a) in mol.atoms.iter().enumerate() {
            active[i] = a.atomic_number > 1;
        }
    }

    // 2-round Weisfeiler-Lehman canonical signature on the active Murcko scaffold subgraph
    let mut color = vec![0u64; n];
    for i in 0..n {
        if !active[i] {
            continue;
        }
        let a = &mol.atoms[i];
        let deg = mol.adj[i].iter().filter(|&&(nb, _)| active[nb]).count() as u64;
        color[i] = mix64(
            u64::from(a.atomic_number)
                ^ (deg << 8)
                ^ (u64::from(a.aromatic) << 16)
                ^ (u64::from(a.in_ring) << 17),
        );
    }
    let mut next_c = vec![0u64; n];
    let mut nbrs: Vec<(u64, u64)> = Vec::with_capacity(8);
    for round in 1u64..=2u64 {
        for i in 0..n {
            if !active[i] {
                continue;
            }
            nbrs.clear();
            for &(nb, bid) in &mol.adj[i] {
                if active[nb] {
                    let bc = mol.bonds.get(bid).map_or(1, |b| bond_code(b.order));
                    nbrs.push((bc, color[nb]));
                }
            }
            nbrs.sort_unstable();
            let mut h = mix64(color[i] ^ (round.wrapping_mul(0x2545_f491_4f6c_dd1d)));
            for &(bc, nc) in &nbrs {
                h = mix64(h ^ mix64(nc.wrapping_add(bc.wrapping_mul(0x9e37_79b9))));
            }
            next_c[i] = h;
        }
        color.copy_from_slice(&next_c);
    }

    let mut active_colors: Vec<u64> = (0..n).filter(|&i| active[i]).map(|i| color[i]).collect();
    active_colors.sort_unstable();
    let mut final_h = mix64(active_colors.len() as u64);
    for c in active_colors {
        final_h = mix64(final_h ^ c);
    }
    final_h
}

/// Extracts the 56 dense physicochemical/graph descriptors + 1024-bit ECFP4 + Murcko scaffold.
#[must_use]
pub fn extract_mol_features(mol: &Molecule, eval: &MoleculeEvaluation) -> MolFeatures {
    let (fp, fp_counts) = compute_ecfp4(mol);
    let scaffold_hash = compute_murcko_scaffold_hash(mol);

    let p = &eval.physicochemical;
    let l = &eval.lipophilicity;
    let s = &eval.water_solubility;
    let d = &eval.druglikeness;
    let m = &eval.medicinal_chemistry;
    let ion = &eval.ionization;
    let saf = &eval.safety;

    let mut n_c = 0.0f32;
    let mut n_n = 0.0f32;
    let mut n_o = 0.0f32;
    let mut n_s = 0.0f32;
    let mut n_p = 0.0f32;
    let mut n_f = 0.0f32;
    let mut n_cl = 0.0f32;
    let mut n_br = 0.0f32;
    let mut n_i = 0.0f32;
    let mut n_basic_n = 0.0f32;
    let mut n_amide_n = 0.0f32;
    let mut n_sulfonamide = 0.0f32;
    let mut n_arom_n = 0.0f32;
    let mut n_nitro = 0.0f32;
    let mut n_pos = 0.0f32;
    let mut n_neg = 0.0f32;
    let mut n_chiral = 0.0f32;
    let mut n_fused = 0.0f32;
    let mut n_sp3_c = 0.0f32;
    let mut n_small_ring = 0.0f32;

    for (idx, a) in mol.atoms.iter().enumerate() {
        if a.atomic_number <= 1 {
            continue;
        }
        match a.atomic_number {
            6 => {
                n_c += 1.0;
                if a.hybridization == Hybridization::Sp3 {
                    n_sp3_c += 1.0;
                }
            }
            7 => {
                n_n += 1.0;
                if a.aromatic {
                    n_arom_n += 1.0;
                }
                let mut adj_carbonyl = false;
                let mut adj_sulfonyl = false;
                let mut adj_o_double = 0usize;
                let mut adj_arom = false;
                for &(nb, bid) in &mol.adj[idx] {
                    if let Some(nba) = mol.atoms.get(nb) {
                        if nba.aromatic {
                            adj_arom = true;
                        }
                        if nba.atomic_number == 8
                            && mol
                                .bonds
                                .get(bid)
                                .is_some_and(|b| b.order == BondOrder::Double)
                        {
                            adj_o_double += 1;
                        }
                        if nba.atomic_number == 6 {
                            let has_co = mol.adj[nb].iter().any(|&(nb2, bid2)| {
                                mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 8)
                                    && mol
                                        .bonds
                                        .get(bid2)
                                        .is_some_and(|b| b.order == BondOrder::Double)
                            });
                            if has_co {
                                adj_carbonyl = true;
                            }
                        }
                        if nba.atomic_number == 16 {
                            let n_so = mol.adj[nb]
                                .iter()
                                .filter(|&&(nb2, bid2)| {
                                    mol.atoms.get(nb2).is_some_and(|x| x.atomic_number == 8)
                                        && mol
                                            .bonds
                                            .get(bid2)
                                            .is_some_and(|b| b.order == BondOrder::Double)
                                })
                                .count();
                            if n_so >= 2 {
                                adj_sulfonyl = true;
                            }
                        }
                    }
                }
                if adj_o_double >= 1 && (a.formal_charge > 0 || adj_o_double >= 2) {
                    n_nitro += 1.0;
                } else if adj_sulfonyl {
                    n_sulfonamide += 1.0;
                } else if adj_carbonyl {
                    n_amide_n += 1.0;
                } else if !a.aromatic && !adj_arom && a.hybridization == Hybridization::Sp3 {
                    n_basic_n += 1.0;
                }
            }
            8 => n_o += 1.0,
            16 => n_s += 1.0,
            15 => n_p += 1.0,
            9 => n_f += 1.0,
            17 => n_cl += 1.0,
            35 => n_br += 1.0,
            53 => n_i += 1.0,
            _ => {}
        }
        if a.formal_charge > 0 {
            n_pos += 1.0;
        } else if a.formal_charge < 0 {
            n_neg += 1.0;
        }
        if a.chiral != ChiralTag::None {
            n_chiral += 1.0;
        }
        if a.ring_sizes.len() >= 2 {
            n_fused += 1.0;
        }
        if a.ring_sizes.iter().any(|&rs| rs == 3 || rs == 4) {
            n_small_ring += 1.0;
        }
    }

    let ha = (p.heavy_atoms.max(1)) as f32;
    let frac_hetero = ((p.heavy_atoms as f32) - n_c).max(0.0) / ha;
    let herg_code = match saf.herg_risk {
        HergRisk::Low => 0.0f32,
        HergRisk::Moderate => 1.0f32,
        HergRisk::High => 2.0f32,
    };

    let dense = [
        p.molecular_weight as f32,
        p.exact_mass as f32,
        p.heavy_atoms as f32,
        p.aromatic_heavy_atoms as f32,
        p.fraction_csp3 as f32,
        p.rotatable_bonds as f32,
        p.h_bond_acceptors as f32,
        p.h_bond_donors as f32,
        p.molar_refractivity as f32,
        p.tpsa as f32,
        p.tpsa_no_only as f32,
        p.dbe as f32,
        p.num_rings as f32,
        p.num_aromatic_rings as f32,
        l.wlogp as f32,
        l.xlogp3 as f32,
        l.mlogp as f32,
        l.consensus_log_p as f32,
        ion.logd_ph74 as f32,
        ion.estimated_charge_ph74 as f32,
        ion.strongest_acidic_pka.unwrap_or(14.0) as f32,
        ion.strongest_basic_pka.unwrap_or(0.0) as f32,
        if ion.strongest_acidic_pka.is_some() {
            1.0
        } else {
            0.0
        },
        if ion.strongest_basic_pka.is_some() {
            1.0
        } else {
            0.0
        },
        s.esol.log_s as f32,
        s.ali.log_s as f32,
        d.qed as f32,
        d.cns_mpo as f32,
        m.synthetic_accessibility as f32,
        d.lipinski.violations as f32,
        m.pains_alerts.len() as f32,
        m.brenk_alerts.len() as f32,
        herg_code,
        saf.ames_alerts.len() as f32,
        saf.dili_alerts.len() as f32,
        n_c,
        n_n,
        n_o,
        n_s,
        n_p,
        n_f,
        n_cl,
        n_br,
        n_i,
        frac_hetero,
        n_basic_n,
        n_amide_n,
        n_sulfonamide,
        n_arom_n,
        n_nitro,
        n_pos,
        n_neg,
        n_chiral,
        n_fused,
        n_sp3_c,
        n_small_ring,
    ];

    MolFeatures {
        dense,
        fp,
        fp_counts,
        scaffold_hash,
    }
}
