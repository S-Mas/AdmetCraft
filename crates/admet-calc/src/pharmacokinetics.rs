use chem_core::{BondOrder, Molecule};

use crate::physicochem::round2;
use crate::types::{AbsorptionLevel, BoiledEggPoint, Lipophilicity, Pharmacokinetics, PhysicochemicalProperties};

/// Daina & Zoete (2016) BOILED-Egg White (GI) ellipse center and rotated semi-axes.
/// Exposed publicly so `admet-report` and `admetcraft` GUI render the exact same mathematical ellipses.
pub const BOILED_EGG_WHITE_CENTER: (f64, f64) = (71.051, 2.292);
pub const BOILED_EGG_WHITE_AXES: (f64, f64) = (71.04, 4.65);
pub const BOILED_EGG_WHITE_ANGLE_RAD: f64 = -0.025; // slight negative tilt in (TPSA, WLOGP)

/// Daina & Zoete (2016) BOILED-Egg Yolk (BBB) ellipse center and semi-axes.
pub const BOILED_EGG_YOLK_CENTER: (f64, f64) = (38.117, 3.177);
pub const BOILED_EGG_YOLK_AXES: (f64, f64) = (39.50, 2.70);
pub const BOILED_EGG_YOLK_ANGLE_RAD: f64 = 0.0;

/// Test whether `(tpsa, wlogp)` lies inside a rotated ellipse.
pub fn point_in_ellipse(
    tpsa: f64,
    wlogp: f64,
    center: (f64, f64),
    axes: (f64, f64),
    angle_rad: f64,
) -> bool {
    let dx = tpsa - center.0;
    let dy = wlogp - center.1;
    let cos_a = angle_rad.cos();
    let sin_a = angle_rad.sin();
    let xr = dx * cos_a + dy * sin_a;
    let yr = -dx * sin_a + dy * cos_a;
    (xr / axes.0).powi(2) + (yr / axes.1).powi(2) <= 1.0
}

/// Evaluate BOILED-Egg coordinates and region membership.
pub fn compute_boiled_egg(
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
    pgp_substrate: bool,
) -> BoiledEggPoint {
    let in_white = point_in_ellipse(
        phys.tpsa,
        lipo.wlogp,
        BOILED_EGG_WHITE_CENTER,
        BOILED_EGG_WHITE_AXES,
        BOILED_EGG_WHITE_ANGLE_RAD,
    );
    let in_yolk = point_in_ellipse(
        phys.tpsa,
        lipo.wlogp,
        BOILED_EGG_YOLK_CENTER,
        BOILED_EGG_YOLK_AXES,
        BOILED_EGG_YOLK_ANGLE_RAD,
    );
    let in_plot_bounds =
        (0.0..=200.0).contains(&phys.tpsa) && (-4.0..=8.0).contains(&lipo.wlogp);

    BoiledEggPoint {
        tpsa: phys.tpsa,
        wlogp: lipo.wlogp,
        in_white,
        in_yolk,
        pgp_substrate,
        in_plot_bounds,
    }
}

/// Compute Pharmacokinetics predictions (BOILED-Egg GI/BBB, P-gp, CYP450, Log Kp).
pub fn compute_pharmacokinetics(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> (Pharmacokinetics, BoiledEggPoint) {
    let pgp = predict_pgp_substrate(mol, phys, lipo);
    let egg = compute_boiled_egg(phys, lipo, pgp);

    let gi_absorption = if egg.in_white || egg.in_yolk {
        AbsorptionLevel::High
    } else {
        AbsorptionLevel::Low
    };
    let bbb_permeant = egg.in_yolk;

    let has_basic_amine = has_aliphatic_amine(mol);
    let has_acid_group = has_carboxylic_or_sulfonic_acid(mol);
    let arom_rings = count_aromatic_rings(mol);

    // Descriptor-calibrated CYP450 inhibition predictors (reflecting known SAR rules):
    // - CYP1A2: planar polyaromatics / low MW / moderate lipophilicity
    let cyp1a2 = !has_acid_group
        && phys.aromatic_heavy_atoms >= 10
        && phys.fraction_csp3 < 0.35
        && phys.molecular_weight <= 420.0
        && lipo.consensus_log_p >= 1.5;

    // - CYP2C19: medium MW, aromatic heterocycles / diaryl / moderately lipophilic
    let cyp2c19 = !has_acid_group
        && phys.molecular_weight >= 240.0
        && phys.molecular_weight <= 500.0
        && lipo.consensus_log_p >= 2.0
        && arom_rings >= 2
        && phys.tpsa <= 105.0;

    // - CYP2C9: lipophilic aromatic molecules (often acidic or sulfonyl/heterocyclic)
    let cyp2c9 = phys.molecular_weight >= 250.0
        && lipo.consensus_log_p >= 2.4
        && phys.aromatic_heavy_atoms >= 10
        && phys.tpsa <= 125.0;

    // - CYP2D6: basic amine + aromatic ring + lipophilic (classical CYP2D6 pharmacophore)
    let cyp2d6 = has_basic_amine
        && !has_acid_group
        && phys.aromatic_heavy_atoms >= 6
        && phys.molecular_weight >= 220.0
        && lipo.consensus_log_p >= 2.0;

    // - CYP3A4: larger MW, flexible, lipophilic, H-bond acceptors >= 3
    let cyp3a4 = phys.molecular_weight >= 340.0
        && lipo.consensus_log_p >= 1.8
        && phys.h_bond_acceptors >= 3
        && phys.aromatic_heavy_atoms >= 10;

    // Potts & Guy (1992) skin permeation coefficient log Kp (cm/s): -2.72 - log10(3600) = -6.2763
    let log_kp = round2(0.71 * lipo.xlogp3 - 0.0061 * phys.molecular_weight - 6.276);

    let pk = Pharmacokinetics {
        gi_absorption,
        bbb_permeant,
        pgp_substrate: pgp,
        cyp1a2_inhibitor: cyp1a2,
        cyp2c19_inhibitor: cyp2c19,
        cyp2c9_inhibitor: cyp2c9,
        cyp2d6_inhibitor: cyp2d6,
        cyp3a4_inhibitor: cyp3a4,
        log_kp_cm_s: log_kp,
    };
    (pk, egg)
}

fn predict_pgp_substrate(
    mol: &Molecule,
    phys: &PhysicochemicalProperties,
    lipo: &Lipophilicity,
) -> bool {
    let has_basic_amine = has_aliphatic_amine(mol);
    // Classical P-gp substrate profile: high MW / high HBA+HBD / cationic amphiphile or macro/polycyclic
    let high_polar_size = phys.molecular_weight >= 400.0 && (phys.h_bond_acceptors + phys.h_bond_donors) >= 5;
    let basic_amphiphile = has_basic_amine && phys.molecular_weight >= 275.0 && (phys.h_bond_donors >= 1 || lipo.consensus_log_p >= 2.0);
    let zwitterion_fluoroquinolone = has_basic_amine && has_carboxylic_or_sulfonic_acid(mol) && phys.molecular_weight >= 280.0;
    high_polar_size || basic_amphiphile || zwitterion_fluoroquinolone
}

fn has_aliphatic_amine(mol: &Molecule) -> bool {
    mol.atoms.iter().enumerate().any(|(idx, a)| {
        if a.atomic_number != 7 || a.aromatic {
            return false;
        }
        // Exclude amide, sulfonamide, nitro, nitrile
        let next_to_double_o = mol.neighbors(idx).any(|nbr| {
            matches!(mol.atoms[nbr].atomic_number, 6 | 16)
                && mol.adj[nbr].iter().any(|&(o, bid)| {
                    o != idx
                        && mol.bonds[bid].order == BondOrder::Double
                        && mol.atoms[o].atomic_number == 8
                })
        });
        let has_multi = mol.adj[idx]
            .iter()
            .any(|&(_, bid)| matches!(mol.bonds[bid].order, BondOrder::Double | BondOrder::Triple));
        let on_arom = mol.neighbors(idx).any(|nbr| mol.atoms[nbr].aromatic);
        !next_to_double_o && !has_multi && !on_arom
    })
}

pub fn has_carboxylic_or_sulfonic_acid(mol: &Molecule) -> bool {
    mol.atoms.iter().enumerate().any(|(idx, a)| {
        if !matches!(a.atomic_number, 6 | 16) {
            return false;
        }
        let has_double_o = mol.adj[idx].iter().any(|&(nbr, bid)| {
            mol.bonds[bid].order == BondOrder::Double && mol.atoms[nbr].atomic_number == 8
        });
        let has_oh = mol.adj[idx].iter().any(|&(nbr, bid)| {
            mol.bonds[bid].order == BondOrder::Single
                && mol.atoms[nbr].atomic_number == 8
                && (mol.atoms[nbr].total_h() > 0 || mol.atoms[nbr].formal_charge < 0)
        });
        has_double_o && has_oh
    })
}

fn count_aromatic_rings(mol: &Molecule) -> usize {
    mol.rings
        .iter()
        .filter(|r| r.iter().all(|&aid| mol.atoms[aid].aromatic))
        .count()
}
