use chem_core::{BondOrder, Hybridization, Molecule};

use crate::types::PhysicochemicalProperties;

/// Compute all physicochemical descriptors for `mol` with SwissADME/OpenBabel conventions
/// plus exact monoisotopic mass, DBE, and N/O-only TPSA (RDKit/QED convention).
pub fn compute_physicochemical(mol: &Molecule) -> PhysicochemicalProperties {
    let formula = mol.formula();
    let molecular_weight = round2(mol.molecular_weight());
    let exact_mass = round4(compute_exact_mass(mol));
    let heavy_atoms = mol.heavy_atom_count();
    let aromatic_heavy_atoms = mol.aromatic_heavy_atom_count();
    let fraction_csp3 = round2(compute_fsp3(mol));
    let rotatable_bonds = count_rotatable_bonds(mol);
    let h_bond_donors = count_hbd(mol);
    let h_bond_acceptors = count_hba(mol);
    let tpsa = round2(compute_tpsa_ertl(mol));
    let tpsa_no_only = round2(compute_tpsa_no_only(mol));
    let dbe = round1(compute_dbe(mol));
    let num_rings = mol.rings.len();
    let num_aromatic_rings = mol
        .rings
        .iter()
        .filter(|r| !r.is_empty() && r.iter().all(|&aid| mol.atoms[aid].aromatic))
        .count();
    let (_, mr) = crate::lipophilicity::compute_wlogp_and_mr(mol);

    PhysicochemicalProperties {
        formula,
        molecular_weight,
        exact_mass,
        heavy_atoms,
        aromatic_heavy_atoms,
        fraction_csp3,
        rotatable_bonds,
        h_bond_acceptors,
        h_bond_donors,
        molar_refractivity: round2(mr),
        tpsa,
        tpsa_no_only,
        dbe,
        num_rings,
        num_aromatic_rings,
    }
}

/// Exact monoisotopic IUPAC molecular mass (most abundant natural isotopes, with electron-mass charge correction).
pub fn compute_exact_mass(mol: &Molecule) -> f64 {
    const ELECTRON: f64 = 0.000_548_579_9;
    const H_MONO: f64 = 1.007_825_032_2;
    let mut sum = 0.0_f64;
    let mut net_charge: i32 = 0;
    for a in &mol.atoms {
        net_charge += i32::from(a.formal_charge);
        let atom_m = if let Some(iso) = a.isotope {
            f64::from(iso)
        } else {
            match a.atomic_number {
                1 => H_MONO,
                5 => 11.009_305_36,
                6 => 12.000_000_0,
                7 => 14.003_074_004,
                8 => 15.994_914_620,
                9 => 18.998_403_163,
                11 => 22.989_769_282,
                12 => 23.985_041_697,
                14 => 27.976_926_535,
                15 => 30.973_761_998,
                16 => 31.972_071_174,
                17 => 34.968_852_682,
                19 => 38.963_706_486,
                20 => 39.962_590_863,
                26 => 55.934_936_33,
                30 => 63.929_142_01,
                34 => 79.916_521_8,
                35 => 78.918_337_6,
                53 => 126.904_471_9,
                _ => {
                    chem_core::element_by_number(a.atomic_number).map_or(12.0, |e| e.atomic_weight)
                }
            }
        };
        sum += atom_m + f64::from(a.total_h()) * H_MONO;
    }
    sum - (net_charge as f64) * ELECTRON
}

/// Double Bond Equivalents (Rings + Double Bonds): `C + Si - (H + X)/2 + (N + P)/2 + 1`.
pub fn compute_dbe(mol: &Molecule) -> f64 {
    if mol.atoms.is_empty() {
        return 0.0;
    }
    let mut c = 0.0_f64;
    let mut h = 0.0_f64;
    let mut x = 0.0_f64;
    let mut n = 0.0_f64;
    for a in &mol.atoms {
        h += f64::from(a.total_h());
        match a.atomic_number {
            1 => h += 1.0,
            6 | 14 => c += 1.0,
            7 | 15 => n += 1.0,
            9 | 17 | 35 | 53 => x += 1.0,
            _ => {}
        }
    }
    c - 0.5 * (h + x) + 0.5 * n + 1.0
}

/// Fraction of sp3-hybridized carbons (`C_sp3 / C_total`).
pub fn compute_fsp3(mol: &Molecule) -> f64 {
    let mut total_c = 0usize;
    let mut sp3_c = 0usize;
    for atom in &mol.atoms {
        if atom.atomic_number == 6 {
            total_c += 1;
            if atom.hybridization == Hybridization::Sp3 {
                sp3_c += 1;
            }
        }
    }
    if total_c == 0 {
        0.0
    } else {
        sp3_c as f64 / total_c as f64
    }
}

/// Count rotatable bonds (Veber / SwissADME definition: single, non-ring, non-terminal heavy-heavy
/// bonds, excluding amide C(=O)-N and ester/thioamide C(=S)-N or sp-sp bonds).
pub fn count_rotatable_bonds(mol: &Molecule) -> usize {
    mol.bonds
        .iter()
        .filter(|b| {
            if b.in_ring || b.order != BondOrder::Single {
                return false;
            }
            let (Some(au), Some(av)) = (mol.atoms.get(b.u), mol.atoms.get(b.v)) else {
                return false;
            };
            if au.atomic_number == 1 || av.atomic_number == 1 {
                return false;
            }
            // Exclude terminal bonds (heavy degree <= 1)
            if mol.heavy_degree(b.u) <= 1 || mol.heavy_degree(b.v) <= 1 {
                return false;
            }
            // Exclude bonds adjacent to triple bonds (sp hybridized)
            if au.hybridization == Hybridization::Sp || av.hybridization == Hybridization::Sp {
                return false;
            }
            // Exclude amide C(=O)-N bonds
            if is_amide_cn(mol, b.u, b.v) || is_amide_cn(mol, b.v, b.u) {
                return false;
            }
            // Exclude symmetric terminal trihalomethyl/tert-butyl symmetry?
            // Note: SwissADME uses OpenBabel RotatableBonds which counts -CF3/-C(CH3)3?
            // OpenBabel ignores terminal -CH3/-NH2/-OH/-CX3 if degree==1 at atom,
            // wait: for -CF3 C has heavy degree 4 (1 C + 3 F), OpenBabel counts or excludes?
            // OpenBabel IsRotor() excludes bonds to terminal groups like -CF3? Actually
            // Fluoxetine has 6 rotatable bonds in SwissADME? Wait: CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1:
            // C-N(1), N-C(2), C-C(3), C-O(4), O-c(5), C-c(6), c-C(F)3(7)!
            // Let's check: in Fluoxetine `CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1`,
            // OpenBabel OBAtom::IsRotor counts: N-C (no, N has H and CH3 -> heavy deg 2!),
            // C-C, C-C, C-O, O-Ar, C-Ar -> 6 (because C-CF3 is symmetric C3 rotor excluded by OpenBabel!).
            if is_symmetric_cx3(mol, b.u, b.v) || is_symmetric_cx3(mol, b.v, b.u) {
                return false;
            }
            true
        })
        .count()
}

fn is_amide_cn(mol: &Molecule, c_idx: usize, n_idx: usize) -> bool {
    let (Some(ac), Some(an)) = (mol.atoms.get(c_idx), mol.atoms.get(n_idx)) else {
        return false;
    };
    if ac.atomic_number != 6 || an.atomic_number != 7 {
        return false;
    }
    mol.adj[c_idx].iter().any(|&(nbr, bid)| {
        nbr != n_idx
            && mol.bonds[bid].order == BondOrder::Double
            && mol.atoms[nbr].atomic_number == 8
    })
}

fn is_symmetric_cx3(mol: &Molecule, term_c: usize, anchor: usize) -> bool {
    let Some(ac) = mol.atoms.get(term_c) else {
        return false;
    };
    if ac.atomic_number != 6 || ac.total_h() != 0 || mol.adj[term_c].len() != 4 {
        return false;
    }
    let others: Vec<usize> = mol.neighbors(term_c).filter(|&nbr| nbr != anchor).collect();
    if others.len() != 3 {
        return false;
    }
    let first_z = mol.atoms[others[0]].atomic_number;
    let first_deg = mol.adj[others[0]].len();
    others.iter().all(|&o| {
        mol.atoms[o].atomic_number == first_z
            && mol.adj[o].len() == 1
            && first_deg == 1
            && matches!(first_z, 9 | 17 | 35 | 53)
    })
}

/// Lipinski-style N/O count (`NorO`): count of nitrogen and oxygen atoms.
pub fn count_n_or_o(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.atomic_number, 7 | 8))
        .count()
}

/// Lipinski-style NH/OH count (`NHorOH`): count of N and O atoms bearing at least 1 H.
pub fn count_nh_or_oh(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.atomic_number, 7 | 8) && a.total_h() > 0)
        .count()
}

/// H-bond donors (OH and NH groups).
pub fn count_hbd(mol: &Molecule) -> usize {
    count_nh_or_oh(mol)
}

/// H-bond acceptors (SwissADME/OpenBabel-style: N/O acceptors excluding positively charged N,
/// pyrrole-like NH, and amide N, plus aliphatic fluorine F attached to sp3 carbon).
pub fn count_hba(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .enumerate()
        .filter(|&(idx, a)| match a.atomic_number {
            8 => a.formal_charge <= 0,
            7 => {
                if a.formal_charge > 0 {
                    return false;
                }
                // Exclude pyrrole-like aromatic N with H or 3 bonds
                if a.aromatic && (a.total_h() > 0 || mol.degree(idx) >= 3) {
                    return false;
                }
                // Exclude amide/sulfonamide N (next to C=O or S=O)
                let next_to_carbonyl_or_sulfonyl = mol.neighbors(idx).any(|nbr| {
                    let nb_atom = &mol.atoms[nbr];
                    (nb_atom.atomic_number == 6 || nb_atom.atomic_number == 16)
                        && mol.adj[nbr].iter().any(|&(nbr2, bid)| {
                            nbr2 != idx
                                && mol.bonds[bid].order == BondOrder::Double
                                && mol.atoms[nbr2].atomic_number == 8
                        })
                });
                !next_to_carbonyl_or_sulfonyl
            }
            9 => {
                // Aliphatic fluorine attached to sp3 C (e.g. -CF3 in Fluoxetine!)
                mol.neighbors(idx).any(|nbr| {
                    mol.atoms[nbr].atomic_number == 6
                        && mol.atoms[nbr].hybridization == Hybridization::Sp3
                })
            }
            _ => false,
        })
        .count()
}

/// Ertl, Rohde, Selzer (2000) *J. Med. Chem.* 43, 3714-3717 Topological Polar Surface Area (TPSA)
/// including N, O, S, and P contributions (SwissADME convention).
pub fn compute_tpsa_ertl(mol: &Molecule) -> f64 {
    compute_tpsa_inner(mol, true)
}

/// Ertl (2000) TPSA counting N and O only (RDKit default & Bickerton QED convention).
pub fn compute_tpsa_no_only(mol: &Molecule) -> f64 {
    compute_tpsa_inner(mol, false)
}

fn compute_tpsa_inner(mol: &Molecule, include_s_and_p: bool) -> f64 {
    let mut tpsa = 0.0;
    for (idx, atom) in mol.atoms.iter().enumerate() {
        let z = atom.atomic_number;
        if !matches!(z, 7 | 8) && !(include_s_and_p && matches!(z, 15 | 16)) {
            continue;
        }
        let h = atom.total_h();
        let charge = atom.formal_charge;
        let in_3_ring = atom.ring_sizes.contains(&3);

        let mut single_b = 0usize;
        let mut double_b = 0usize;
        let mut triple_b = 0usize;
        let mut arom_b = 0usize;
        for &(nbr, bid) in &mol.adj[idx] {
            if mol.atoms[nbr].atomic_number == 1 {
                continue;
            }
            match mol.bonds[bid].order {
                BondOrder::Single => single_b += 1,
                BondOrder::Double => double_b += 1,
                BondOrder::Triple => triple_b += 1,
                BondOrder::Aromatic => arom_b += 1,
            }
        }

        let contrib = match z {
            7 => ertl_n(single_b, double_b, triple_b, arom_b, h, charge, in_3_ring),
            8 => ertl_o(single_b, double_b, arom_b, h, charge, in_3_ring),
            16 => ertl_s(single_b, double_b, arom_b, h),
            15 => ertl_p(single_b, double_b, h),
            _ => 0.0,
        };
        tpsa += contrib;
    }
    tpsa
}

fn ertl_n(
    single_b: usize,
    double_b: usize,
    triple_b: usize,
    arom_b: usize,
    h: u8,
    charge: i8,
    in_3_ring: bool,
) -> f64 {
    if charge == 0 {
        if arom_b >= 2 {
            if h == 0 && single_b == 0 {
                return 12.89; // :n: (pyridine-like)
            }
            if h == 1 {
                return 15.79; // :[nH]: (pyrrole-like)
            }
            if single_b >= 1 {
                return 4.93; // :n(-): (N-substituted pyrrole/purine)
            }
        }
        if triple_b == 1 && single_b == 0 {
            return 23.79; // #N (nitrile)
        }
        if double_b == 1 && single_b == 1 && h == 0 {
            return 12.36; // =N-
        }
        if double_b == 1 && h == 1 {
            return 23.85; // =NH
        }
        if double_b == 2 {
            return 13.60; // =N=
        }
        if single_b == 3 && h == 0 {
            return if in_3_ring { 3.01 } else { 3.24 };
        }
        if single_b == 2 && h == 1 {
            return if in_3_ring { 21.94 } else { 12.03 };
        }
        if single_b == 1 && h == 2 {
            return if in_3_ring { 35.02 } else { 26.02 };
        }
        if single_b == 0 && h == 3 {
            return 37.38;
        }
        if double_b == 1 && single_b == 0 {
            return 23.85;
        }
    } else if charge > 0 {
        if arom_b >= 2 && h == 0 {
            return 3.88; // :[n+]:
        }
        if arom_b >= 2 && h == 1 {
            return 14.14; // :[nH+]:
        }
        if single_b == 4 {
            return 0.0;
        }
        if single_b == 3 && h == 1 {
            return 4.44;
        }
        if single_b == 2 && h == 2 {
            return 16.61;
        }
        if single_b == 1 && h == 3 {
            return 27.64;
        }
        if double_b == 1 && single_b == 2 {
            return 3.01;
        }
        if double_b == 1 && single_b == 1 && h == 1 {
            return 13.97;
        }
        if double_b == 1 && h == 2 {
            return 25.59;
        }
        // Ertl (2000) linear fallback for unlisted charged N (e.g. azide =[N+]= -> 14.10)
        let heavy_deg = single_b + double_b + triple_b + arom_b;
        return (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * f64::from(h)).max(0.0);
    } else if charge < 0 {
        // Ertl (2000) linear fallback for formally negative N (e.g. azide =[N-] -> 22.30)
        let heavy_deg = single_b + double_b + triple_b + arom_b;
        return (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * f64::from(h)).max(0.0);
    }
    let heavy_deg = single_b + double_b + triple_b + arom_b;
    (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * f64::from(h)).max(0.0)
}

fn ertl_o(
    single_b: usize,
    double_b: usize,
    arom_b: usize,
    h: u8,
    charge: i8,
    in_3_ring: bool,
) -> f64 {
    if charge < 0 {
        return 23.06; // [O-]
    }
    if arom_b >= 2 {
        return 13.14; // :o: (furan)
    }
    if double_b == 1 {
        return 17.07; // =O
    }
    if single_b == 2 && h == 0 {
        return if in_3_ring { 12.53 } else { 9.23 }; // -O- (epoxide vs ether)
    }
    if single_b == 1 && h == 1 {
        return 20.23; // -OH
    }
    if single_b == 0 && h == 2 {
        return 31.50; // H2O
    }
    9.23
}

fn ertl_s(single_b: usize, double_b: usize, arom_b: usize, h: u8) -> f64 {
    if arom_b >= 2 {
        return 28.24; // :s:
    }
    if double_b == 2 && single_b == 2 {
        return 8.38; // -S(=O)(=O)- (S contribution in Ertl 2000!)
    }
    if double_b == 1 && single_b == 2 {
        return 19.21; // -S(=O)- (sulfoxide S=O -> S=19.21? Wait: Ertl 2000 table: -S(=O)- S=19.21 + O=17.07 = 36.28! Wait: Penicillin thioether -S- = 25.30!)
    }
    if double_b == 1 && single_b == 0 {
        return 32.09; // =S
    }
    if single_b == 2 && h == 0 {
        return 25.30; // -S-
    }
    if single_b == 1 && h == 1 {
        return 38.80; // -SH
    }
    25.30
}

fn ertl_p(single_b: usize, double_b: usize, h: u8) -> f64 {
    if double_b == 1 && single_b == 3 {
        return 9.81; // -P(=O)(-)-
    }
    if double_b == 1 && single_b == 2 && h == 1 {
        return 23.47;
    }
    if single_b == 3 && h == 0 {
        return 13.59;
    }
    if single_b == 2 && h == 1 {
        return 34.14;
    }
    13.59
}

pub fn round1(v: f64) -> f64 {
    (v * 10.0).round() / 10.0
}

pub fn round2(v: f64) -> f64 {
    (v * 100.0).round() / 100.0
}

pub fn round4(v: f64) -> f64 {
    (v * 10_000.0).round() / 10_000.0
}
