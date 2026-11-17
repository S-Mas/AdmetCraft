use std::f64::consts::PI;
use std::fmt::Write as _;

use admet_calc::MoleculeEvaluation;
use admet_calc::pharmacokinetics::{
    BOILED_EGG_WHITE_ANGLE_RAD, BOILED_EGG_WHITE_AXES, BOILED_EGG_WHITE_CENTER,
    BOILED_EGG_YOLK_ANGLE_RAD, BOILED_EGG_YOLK_AXES, BOILED_EGG_YOLK_CENTER,
};
use chem_core::{BondOrder, ensure_2d_coords, parse_smiles};

use crate::ReportError;

const PAGE_W: f64 = 595.0;
const PAGE_H: f64 = 842.0;

/// Generate a publication-ready multi-page PDF 1.4 report in pure Rust.
pub fn export_pdf_bytes(evals: &[MoleculeEvaluation]) -> Result<Vec<u8>, ReportError> {
    if evals.is_empty() {
        return Err(ReportError::EmptyReport);
    }

    let mut page_streams: Vec<String> = Vec::new();
    page_streams.push(build_summary_page(evals));
    for eval in evals {
        page_streams.push(build_molecule_page(eval));
    }

    Ok(assemble_pdf(&page_streams))
}

fn build_summary_page(evals: &[MoleculeEvaluation]) -> String {
    let mut s = String::new();
    // Header banner
    draw_rect_fill(&mut s, 0.0, PAGE_H - 64.0, PAGE_W, 64.0, (0.22, 0.22, 0.24));
    draw_text(
        &mut s,
        "F2",
        18.0,
        40.0,
        PAGE_H - 36.0,
        (1.0, 1.0, 1.0),
        "AdmetCraft  |  SwissADME-Compatible ADMET Report",
    );
    draw_text(
        &mut s,
        "F1",
        9.5,
        40.0,
        PAGE_H - 52.0,
        (0.85, 0.85, 0.85),
        &format!("Local-First Pure-Rust Report  |  {} molecule(s) evaluated", evals.len()),
    );

    // Summary table
    let mut y = PAGE_H - 92.0;
    draw_text(&mut s, "F2", 11.5, 40.0, y, (0.15, 0.15, 0.15), "Summary Table");
    y -= 16.0;
    draw_rect_fill(&mut s, 40.0, y - 4.0, 515.0, 16.0, (0.92, 0.92, 0.94));
    draw_text(&mut s, "F2", 8.0, 44.0, y, (0.1, 0.1, 0.1), "Molecule");
    draw_text(&mut s, "F2", 8.0, 155.0, y, (0.1, 0.1, 0.1), "Formula");
    draw_text(&mut s, "F2", 8.0, 225.0, y, (0.1, 0.1, 0.1), "MW");
    draw_text(&mut s, "F2", 8.0, 270.0, y, (0.1, 0.1, 0.1), "TPSA");
    draw_text(&mut s, "F2", 8.0, 315.0, y, (0.1, 0.1, 0.1), "cLogP");
    draw_text(&mut s, "F2", 8.0, 360.0, y, (0.1, 0.1, 0.1), "ESOL");
    draw_text(&mut s, "F2", 8.0, 405.0, y, (0.1, 0.1, 0.1), "GI");
    draw_text(&mut s, "F2", 8.0, 442.0, y, (0.1, 0.1, 0.1), "BBB");
    draw_text(&mut s, "F2", 8.0, 478.0, y, (0.1, 0.1, 0.1), "Lipinski");
    draw_text(&mut s, "F2", 8.0, 525.0, y, (0.1, 0.1, 0.1), "SA");
    y -= 15.0;

    for (idx, e) in evals.iter().take(18).enumerate() {
        if idx % 2 == 1 {
            draw_rect_fill(&mut s, 40.0, y - 3.0, 515.0, 13.0, (0.97, 0.97, 0.98));
        }
        draw_text(&mut s, "F1", 8.0, 44.0, y, (0.1, 0.1, 0.1), &truncate(&e.name, 20));
        draw_text(&mut s, "F1", 8.0, 155.0, y, (0.1, 0.1, 0.1), &truncate(&e.physicochemical.formula, 12));
        draw_text(&mut s, "F1", 8.0, 225.0, y, (0.1, 0.1, 0.1), &format!("{:.1}", e.physicochemical.molecular_weight));
        draw_text(&mut s, "F1", 8.0, 270.0, y, (0.1, 0.1, 0.1), &format!("{:.1}", e.physicochemical.tpsa));
        draw_text(&mut s, "F1", 8.0, 315.0, y, (0.1, 0.1, 0.1), &format!("{:.2}", e.lipophilicity.consensus_log_p));
        draw_text(&mut s, "F1", 8.0, 360.0, y, (0.1, 0.1, 0.1), &format!("{:.2}", e.water_solubility.esol.log_s));
        draw_text(&mut s, "F1", 8.0, 405.0, y, (0.1, 0.1, 0.1), e.pharmacokinetics.gi_absorption.as_str());
        draw_text(&mut s, "F1", 8.0, 442.0, y, (0.1, 0.1, 0.1), if e.pharmacokinetics.bbb_permeant { "Yes" } else { "No" });
        draw_text(&mut s, "F1", 8.0, 478.0, y, (0.1, 0.1, 0.1), &format!("{} ({})", if e.druglikeness.lipinski.pass { "Yes" } else { "No" }, e.druglikeness.lipinski.violations));
        draw_text(&mut s, "F1", 8.0, 525.0, y, (0.1, 0.1, 0.1), &format!("{:.2}", e.medicinal_chemistry.synthetic_accessibility));
        y -= 13.0;
    }

    // BOILED-Egg vector plot on lower half of summary page
    let egg_x = 65.0;
    let egg_y = 75.0;
    let egg_w = 370.0;
    let egg_h = 250.0;
    draw_text(&mut s, "F2", 11.5, 40.0, egg_y + egg_h + 20.0, (0.15, 0.15, 0.15), "BOILED-Egg (WLOGP vs TPSA)");
    draw_rect_fill_stroke(&mut s, (egg_x, egg_y, egg_w, egg_h), (0.85, 0.85, 0.85), (0.35, 0.35, 0.35), 0.8);

    let map_x = |tpsa: f64| egg_x + (tpsa / 200.0) * egg_w;
    let map_y = |wlogp: f64| egg_y + ((wlogp + 4.0) / 12.0) * egg_h;

    let white_poly = sample_ellipse_pdf(
        BOILED_EGG_WHITE_CENTER,
        BOILED_EGG_WHITE_AXES,
        BOILED_EGG_WHITE_ANGLE_RAD,
        egg_x,
        egg_x + egg_w,
        egg_y,
        egg_y + egg_h,
        &map_x,
        &map_y,
    );
    draw_polygon(&mut s, &white_poly, Some((1.0, 1.0, 1.0)), Some((0.5, 0.5, 0.5)), 0.8);

    let yolk_poly = sample_ellipse_pdf(
        BOILED_EGG_YOLK_CENTER,
        BOILED_EGG_YOLK_AXES,
        BOILED_EGG_YOLK_ANGLE_RAD,
        egg_x,
        egg_x + egg_w,
        egg_y,
        egg_y + egg_h,
        &map_x,
        &map_y,
    );
    draw_polygon(&mut s, &yolk_poly, Some((0.96, 0.77, 0.10)), Some((0.75, 0.58, 0.0)), 0.8);

    for t in [0, 50, 100, 150, 200] {
        let px = map_x(f64::from(t));
        draw_line(&mut s, px, egg_y, px, egg_y - 4.0, (0.3, 0.3, 0.3), 0.7);
        draw_text(&mut s, "F1", 7.5, px - 6.0, egg_y - 13.0, (0.2, 0.2, 0.2), &t.to_string());
    }
    for lp in [-4, -2, 0, 2, 4, 6, 8] {
        let py = map_y(f64::from(lp));
        draw_line(&mut s, egg_x - 4.0, py, egg_x, py, (0.3, 0.3, 0.3), 0.7);
        draw_text(&mut s, "F1", 7.5, egg_x - 18.0, py - 2.5, (0.2, 0.2, 0.2), &lp.to_string());
    }
    draw_text(&mut s, "F2", 8.5, egg_x + egg_w * 0.42, egg_y - 26.0, (0.15, 0.15, 0.15), "TPSA (A^2)");
    draw_text(&mut s, "F2", 8.5, 22.0, egg_y + egg_h * 0.5, (0.15, 0.15, 0.15), "WLOGP");

    for e in evals {
        if !e.boiled_egg.in_plot_bounds {
            continue;
        }
        let px = map_x(e.boiled_egg.tpsa);
        let py = map_y(e.boiled_egg.wlogp);
        let col = if e.boiled_egg.pgp_substrate {
            (0.11, 0.49, 0.84)
        } else {
            (0.88, 0.19, 0.19)
        };
        draw_circle(&mut s, px, py, 3.2, col);
        draw_text(&mut s, "F1", 6.5, px + 4.5, py - 2.0, (0.1, 0.1, 0.1), &truncate(&e.name, 14));
    }

    // Legend box
    let lx = egg_x + egg_w + 15.0;
    draw_rect_fill_stroke(&mut s, (lx, egg_y + 110.0, 105.0, 110.0), (0.98, 0.98, 0.98), (0.78, 0.78, 0.78), 0.6);
    draw_rect_fill_stroke(&mut s, (lx + 8.0, egg_y + 195.0, 10.0, 10.0), (1.0, 1.0, 1.0), (0.5, 0.5, 0.5), 0.6);
    draw_text(&mut s, "F1", 8.0, lx + 24.0, egg_y + 197.0, (0.1, 0.1, 0.1), "HIA (White)");
    draw_rect_fill_stroke(&mut s, (lx + 8.0, egg_y + 175.0, 10.0, 10.0), (0.96, 0.77, 0.10), (0.75, 0.58, 0.0), 0.6);
    draw_text(&mut s, "F1", 8.0, lx + 24.0, egg_y + 177.0, (0.1, 0.1, 0.1), "BBB (Yolk)");
    draw_circle(&mut s, lx + 13.0, egg_y + 158.0, 3.5, (0.11, 0.49, 0.84));
    draw_text(&mut s, "F1", 8.0, lx + 24.0, egg_y + 155.0, (0.1, 0.1, 0.1), "PGP+");
    draw_circle(&mut s, lx + 13.0, egg_y + 138.0, 3.5, (0.88, 0.19, 0.19));
    draw_text(&mut s, "F1", 8.0, lx + 24.0, egg_y + 135.0, (0.1, 0.1, 0.1), "PGP-");

    s
}

fn build_molecule_page(eval: &MoleculeEvaluation) -> String {
    let mut s = String::new();
    draw_rect_fill(&mut s, 0.0, PAGE_H - 54.0, PAGE_W, 54.0, (0.24, 0.24, 0.27));
    draw_text(&mut s, "F2", 15.0, 40.0, PAGE_H - 28.0, (1.0, 1.0, 1.0), &truncate(&eval.name, 52));
    draw_text(
        &mut s,
        "F1",
        8.5,
        40.0,
        PAGE_H - 43.0,
        (0.86, 0.86, 0.86),
        &format!("SMILES: {}", truncate(&eval.canonical_smiles, 86)),
    );

    // Top row: 2D Structure box (left) & Bioavailability Radar (right)
    draw_rect_fill_stroke(&mut s, (40.0, PAGE_H - 250.0, 250.0, 180.0), (0.99, 0.99, 0.99), (0.82, 0.82, 0.82), 0.7);
    draw_text(&mut s, "F2", 9.0, 48.0, PAGE_H - 83.0, (0.3, 0.3, 0.3), "2D Structure");
    draw_2d_mol_pdf(&mut s, &eval.canonical_smiles, 48.0, PAGE_H - 242.0, 234.0, 150.0);

    draw_rect_fill_stroke(&mut s, (305.0, PAGE_H - 250.0, 250.0, 180.0), (0.99, 0.99, 0.99), (0.82, 0.82, 0.82), 0.7);
    draw_text(&mut s, "F2", 9.0, 313.0, PAGE_H - 83.0, (0.3, 0.3, 0.3), "Bioavailability Radar");
    draw_radar_pdf(&mut s, eval, 430.0, PAGE_H - 165.0, 56.0);

    // Two-column property cards
    let p = &eval.physicochemical;
    let l = &eval.lipophilicity;
    let w = &eval.water_solubility;
    let pk = &eval.pharmacokinetics;
    let d = &eval.druglikeness;
    let mc = &eval.medicinal_chemistry;

    let mut y_left = PAGE_H - 275.0;
    y_left = draw_kv_section(
        &mut s,
        40.0,
        y_left,
        250.0,
        "Physicochemical Properties",
        &[
            ("Formula", p.formula.clone()),
            ("Molecular weight", format!("{:.2} g/mol", p.molecular_weight)),
            ("Num. heavy atoms", p.heavy_atoms.to_string()),
            ("Num. arom. heavy atoms", p.aromatic_heavy_atoms.to_string()),
            ("Fraction Csp3", format!("{:.2}", p.fraction_csp3)),
            ("Num. rotatable bonds", p.rotatable_bonds.to_string()),
            ("Num. H-bond acceptors", p.h_bond_acceptors.to_string()),
            ("Num. H-bond donors", p.h_bond_donors.to_string()),
            ("Molar Refractivity", format!("{:.2}", p.molar_refractivity)),
            ("TPSA", format!("{:.2} A^2", p.tpsa)),
        ],
    );
    y_left -= 12.0;
    y_left = draw_kv_section(
        &mut s,
        40.0,
        y_left,
        250.0,
        "Lipophilicity (Log Po/w)",
        &[
            ("iLOGP", format!("{:.2}", l.ilogp)),
            ("XLOGP3", format!("{:.2}", l.xlogp3)),
            ("WLOGP", format!("{:.2}", l.wlogp)),
            ("MLOGP", format!("{:.2}", l.mlogp)),
            ("Silicos-IT Log P", format!("{:.2}", l.silicos_it)),
            ("Consensus Log Po/w", format!("{:.2}", l.consensus_log_p)),
        ],
    );
    y_left -= 12.0;
    let _ = draw_kv_section(
        &mut s,
        40.0,
        y_left,
        250.0,
        "Water Solubility",
        &[
            ("Log S (ESOL)", format!("{:.2} ({})", w.esol.log_s, w.esol.class.as_str())),
            ("Log S (Ali)", format!("{:.2} ({})", w.ali.log_s, w.ali.class.as_str())),
            ("Log S (Silicos-IT)", format!("{:.2} ({})", w.silicos_it.log_s, w.silicos_it.class.as_str())),
        ],
    );

    let mut y_right = PAGE_H - 275.0;
    y_right = draw_kv_section(
        &mut s,
        305.0,
        y_right,
        250.0,
        "Pharmacokinetics",
        &[
            ("GI absorption", pk.gi_absorption.as_str().to_string()),
            ("BBB permeant", yes_no(pk.bbb_permeant).to_string()),
            ("P-gp substrate", yes_no(pk.pgp_substrate).to_string()),
            ("CYP1A2 inhibitor", yes_no(pk.cyp1a2_inhibitor).to_string()),
            ("CYP2C19 inhibitor", yes_no(pk.cyp2c19_inhibitor).to_string()),
            ("CYP2C9 inhibitor", yes_no(pk.cyp2c9_inhibitor).to_string()),
            ("CYP2D6 inhibitor", yes_no(pk.cyp2d6_inhibitor).to_string()),
            ("CYP3A4 inhibitor", yes_no(pk.cyp3a4_inhibitor).to_string()),
            ("Log Kp (skin perm.)", format!("{:.2} cm/s", pk.log_kp_cm_s)),
        ],
    );
    y_right -= 12.0;
    y_right = draw_kv_section(
        &mut s,
        305.0,
        y_right,
        250.0,
        "Druglikeness",
        &[
            ("Lipinski", format!("{} ({} viol.)", yes_no(d.lipinski.pass), d.lipinski.violations)),
            ("Ghose", format!("{} ({} viol.)", yes_no(d.ghose.pass), d.ghose.violations)),
            ("Veber", format!("{} ({} viol.)", yes_no(d.veber.pass), d.veber.violations)),
            ("Egan", format!("{} ({} viol.)", yes_no(d.egan.pass), d.egan.violations)),
            ("Muegge", format!("{} ({} viol.)", yes_no(d.muegge.pass), d.muegge.violations)),
            ("Bioavailability Score", format!("{:.2}", d.bioavailability_score)),
        ],
    );
    y_right -= 12.0;
    let pains_desc = if mc.pains_alerts.is_empty() {
        "0 alerts".to_string()
    } else {
        format!("{} ({})", mc.pains_alerts.len(), truncate(&mc.pains_alerts.join(", "), 24))
    };
    let brenk_desc = if mc.brenk_alerts.is_empty() {
        "0 alerts".to_string()
    } else {
        format!("{} ({})", mc.brenk_alerts.len(), truncate(&mc.brenk_alerts.join(", "), 24))
    };
    let _ = draw_kv_section(
        &mut s,
        305.0,
        y_right,
        250.0,
        "Medicinal Chemistry",
        &[
            ("PAINS", pains_desc),
            ("Brenk", brenk_desc),
            ("Leadlikeness", format!("{} ({} viol.)", yes_no(mc.leadlikeness.pass), mc.leadlikeness.violations)),
            ("Synthetic accessibility", format!("{:.2}", mc.synthetic_accessibility)),
        ],
    );

    s
}

fn draw_kv_section(
    s: &mut String,
    x: f64,
    mut y: f64,
    w: f64,
    title: &str,
    rows: &[(&str, String)],
) -> f64 {
    draw_rect_fill(s, x, y - 4.0, w, 15.0, (0.36, 0.36, 0.38));
    draw_text(s, "F2", 9.0, x + 6.0, y, (1.0, 1.0, 1.0), title);
    y -= 15.0;
    for (idx, (k, v)) in rows.iter().enumerate() {
        if idx % 2 == 0 {
            draw_rect_fill(s, x, y - 3.0, w, 12.5, (0.96, 0.96, 0.97));
        }
        draw_text(s, "F1", 8.0, x + 6.0, y, (0.2, 0.2, 0.2), k);
        draw_text(s, "F2", 8.0, x + w * 0.54, y, (0.08, 0.08, 0.08), &truncate(v, 24));
        y -= 12.5;
    }
    y
}

fn draw_2d_mol_pdf(s: &mut String, smiles: &str, x: f64, y: f64, w: f64, h: f64) {
    let Ok(mut mol) = parse_smiles(smiles) else {
        return;
    };
    ensure_2d_coords(&mut mol);
    let coords: Vec<[f64; 2]> = mol
        .atoms
        .iter()
        .filter_map(|a| a.coords_2d)
        .collect();
    if coords.len() != mol.atoms.len() || coords.is_empty() {
        return;
    }
    let min_x = coords.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
    let max_x = coords.iter().map(|c| c[0]).fold(f64::NEG_INFINITY, f64::max);
    let min_y = coords.iter().map(|c| c[1]).fold(f64::INFINITY, f64::min);
    let max_y = coords.iter().map(|c| c[1]).fold(f64::NEG_INFINITY, f64::max);
    let span_x = (max_x - min_x).max(1.0);
    let span_y = (max_y - min_y).max(1.0);
    let scale = (w * 0.78 / span_x).min(h * 0.78 / span_y).min(24.0);
    let cx = x + w * 0.5;
    let cy = y + h * 0.5;
    let mid_x = 0.5 * (min_x + max_x);
    let mid_y = 0.5 * (min_y + max_y);

    let pt = |i: usize| -> (f64, f64) {
        (
            cx + (coords[i][0] - mid_x) * scale,
            cy + (coords[i][1] - mid_y) * scale,
        )
    };

    for b in &mol.bonds {
        let (x1, y1) = pt(b.u);
        let (x2, y2) = pt(b.v);
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = (dx * dx + dy * dy).sqrt().max(1e-6);
        let nx = -dy / len * 1.8;
        let ny = dx / len * 1.8;
        match b.order {
            BondOrder::Single => draw_line(s, x1, y1, x2, y2, (0.15, 0.15, 0.15), 1.0),
            BondOrder::Double | BondOrder::Aromatic => {
                draw_line(s, x1 + nx, y1 + ny, x2 + nx, y2 + ny, (0.15, 0.15, 0.15), 0.9);
                draw_line(s, x1 - nx, y1 - ny, x2 - nx, y2 - ny, (0.15, 0.15, 0.15), 0.9);
            }
            BondOrder::Triple => {
                draw_line(s, x1, y1, x2, y2, (0.15, 0.15, 0.15), 0.9);
                draw_line(s, x1 + nx * 1.4, y1 + ny * 1.4, x2 + nx * 1.4, y2 + ny * 1.4, (0.15, 0.15, 0.15), 0.8);
                draw_line(s, x1 - nx * 1.4, y1 - ny * 1.4, x2 - nx * 1.4, y2 - ny * 1.4, (0.15, 0.15, 0.15), 0.8);
            }
        }
    }

    for (i, a) in mol.atoms.iter().enumerate() {
        if a.atomic_number == 6 {
            continue;
        }
        let (ax, ay) = pt(i);
        let col = match a.atomic_number {
            7 => (0.10, 0.35, 0.85),
            8 => (0.82, 0.12, 0.12),
            16 => (0.72, 0.55, 0.05),
            9 | 17 => (0.12, 0.62, 0.25),
            35 => (0.55, 0.18, 0.12),
            _ => (0.3, 0.3, 0.3),
        };
        draw_rect_fill(s, ax - 4.5, ay - 3.5, 9.0, 7.5, (0.99, 0.99, 0.99));
        draw_text(s, "F2", 7.0, ax - 3.2, ay - 2.3, col, &a.symbol);
    }
}

fn draw_radar_pdf(s: &mut String, eval: &MoleculeEvaluation, cx: f64, cy: f64, max_r: f64) {
    let hex = |r: f64| -> Vec<(f64, f64)> {
        (0..6)
            .map(|i| {
                let ang = PI / 2.0 - (i as f64) * PI / 3.0;
                (cx + r * ang.cos(), cy + r * ang.sin())
            })
            .collect()
    };
    draw_polygon(s, &hex(max_r * 0.85), Some((0.99, 0.89, 0.90)), Some((0.95, 0.65, 0.70)), 0.7);
    draw_polygon(s, &hex(max_r * 0.25), Some((0.99, 0.99, 0.99)), Some((0.95, 0.65, 0.70)), 0.7);
    for frac in [0.5, 1.0] {
        draw_polygon(s, &hex(max_r * frac), None, Some((0.80, 0.80, 0.80)), 0.5);
    }
    let labels = ["LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"];
    for (i, label) in labels.iter().enumerate() {
        let ang = PI / 2.0 - (i as f64) * PI / 3.0;
        let x2 = cx + max_r * ang.cos();
        let y2 = cy + max_r * ang.sin();
        draw_line(s, cx, cy, x2, y2, (0.80, 0.80, 0.80), 0.5);
        let lx = cx + (max_r + 12.0) * ang.cos() - 10.0;
        let ly = cy + (max_r + 12.0) * ang.sin() - 2.5;
        draw_text(s, "F2", 6.5, lx, ly, (0.25, 0.25, 0.25), label);
    }
    let pts: Vec<(f64, f64)> = eval
        .bioavailability_radar
        .normalized_axes
        .iter()
        .enumerate()
        .map(|(i, &val)| {
            let ang = PI / 2.0 - (i as f64) * PI / 3.0;
            let r = max_r * val.clamp(0.05, 1.0);
            (cx + r * ang.cos(), cy + r * ang.sin())
        })
        .collect();
    draw_polygon(s, &pts, None, Some((0.79, 0.16, 0.16)), 1.5);
    for (px, py) in pts {
        draw_circle(s, px, py, 2.2, (0.79, 0.16, 0.16));
    }
}

#[allow(clippy::too_many_arguments)]
fn sample_ellipse_pdf(
    center: (f64, f64),
    axes: (f64, f64),
    angle: f64,
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
    map_x: &dyn Fn(f64) -> f64,
    map_y: &dyn Fn(f64) -> f64,
) -> Vec<(f64, f64)> {
    let steps = 64;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    (0..steps)
        .map(|i| {
            let t = 2.0 * PI * (i as f64) / (steps as f64);
            let lx = axes.0 * t.cos();
            let ly = axes.1 * t.sin();
            let gx = center.0 + lx * cos_a - ly * sin_a;
            let gy = center.1 + lx * sin_a + ly * cos_a;
            (map_x(gx).clamp(min_x, max_x), map_y(gy).clamp(min_y, max_y))
        })
        .collect()
}

fn draw_rect_fill(s: &mut String, x: f64, y: f64, w: f64, h: f64, rgb: (f64, f64, f64)) {
    let _ = writeln!(s, "{:.3} {:.3} {:.3} rg {:.2} {:.2} {:.2} {:.2} re f", rgb.0, rgb.1, rgb.2, x, y, w, h);
}

fn draw_rect_fill_stroke(
    s: &mut String,
    rect: (f64, f64, f64, f64),
    fill: (f64, f64, f64),
    stroke: (f64, f64, f64),
    lw: f64,
) {
    let (x, y, w, h) = rect;
    let _ = writeln!(
        s,
        "{:.3} {:.3} {:.3} rg {:.3} {:.3} {:.3} RG {:.2} w {:.2} {:.2} {:.2} {:.2} re B",
        fill.0, fill.1, fill.2, stroke.0, stroke.1, stroke.2, lw, x, y, w, h
    );
}

fn draw_line(s: &mut String, x1: f64, y1: f64, x2: f64, y2: f64, rgb: (f64, f64, f64), lw: f64) {
    let _ = writeln!(
        s,
        "{:.3} {:.3} {:.3} RG {:.2} w {:.2} {:.2} m {:.2} {:.2} l S",
        rgb.0, rgb.1, rgb.2, lw, x1, y1, x2, y2
    );
}

fn draw_polygon(
    s: &mut String,
    pts: &[(f64, f64)],
    fill: Option<(f64, f64, f64)>,
    stroke: Option<(f64, f64, f64)>,
    lw: f64,
) {
    let Some(&(x0, y0)) = pts.first() else {
        return;
    };
    if let Some(f) = fill {
        let _ = write!(s, "{:.3} {:.3} {:.3} rg ", f.0, f.1, f.2);
    }
    if let Some(st) = stroke {
        let _ = write!(s, "{:.3} {:.3} {:.3} RG {:.2} w ", st.0, st.1, st.2, lw);
    }
    let _ = write!(s, "{x0:.2} {y0:.2} m ");
    for &(px, py) in &pts[1..] {
        let _ = write!(s, "{px:.2} {py:.2} l ");
    }
    let op = match (fill.is_some(), stroke.is_some()) {
        (true, true) => "b",
        (true, false) => "f",
        (false, true) => "s",
        (false, false) => "n",
    };
    let _ = writeln!(s, "{op}");
}

fn draw_circle(s: &mut String, cx: f64, cy: f64, r: f64, fill: (f64, f64, f64)) {
    let pts: Vec<(f64, f64)> = (0..16)
        .map(|i| {
            let a = 2.0 * PI * (i as f64) / 16.0;
            (cx + r * a.cos(), cy + r * a.sin())
        })
        .collect();
    draw_polygon(s, &pts, Some(fill), Some((0.1, 0.1, 0.1)), 0.5);
}

fn draw_text(
    s: &mut String,
    font: &str,
    size: f64,
    x: f64,
    y: f64,
    rgb: (f64, f64, f64),
    text: &str,
) {
    let escaped = escape_pdf_str(text);
    let _ = writeln!(
        s,
        "BT /{font} {size:.1} Tf {:.3} {:.3} {:.3} rg {x:.2} {y:.2} Td ({escaped}) Tj ET",
        rgb.0, rgb.1, rgb.2
    );
}

fn escape_pdf_str(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            c if c.is_ascii() && !c.is_ascii_control() => out.push(c),
            _ => out.push('?'),
        }
    }
    out
}

fn truncate(s: &str, max_len: usize) -> String {
    if s.chars().count() <= max_len {
        s.to_string()
    } else {
        let mut out: String = s.chars().take(max_len.saturating_sub(3)).collect();
        out.push_str("...");
        out
    }
}

fn yes_no(b: bool) -> &'static str {
    if b { "Yes" } else { "No" }
}

fn assemble_pdf(page_streams: &[String]) -> Vec<u8> {
    // Object numbering:
    // 1: Catalog
    // 2: Pages
    // 3: Font F1 (Helvetica)
    // 4: Font F2 (Helvetica-Bold)
    // For each page i (0..n):
    //   page_obj = 5 + 2*i
    //   content_obj = 6 + 2*i
    let mut out = Vec::<u8>::new();
    out.extend_from_slice(b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");

    let mut offsets: Vec<usize> = Vec::new();
    let push_obj = |buf: &mut Vec<u8>, offsets: &mut Vec<usize>, id: usize, body: &[u8]| {
        offsets.push(buf.len());
        let header = format!("{id} 0 obj\n");
        buf.extend_from_slice(header.as_bytes());
        buf.extend_from_slice(body);
        buf.extend_from_slice(b"\nendobj\n");
    };

    push_obj(&mut out, &mut offsets, 1, b"<< /Type /Catalog /Pages 2 0 R >>");

    let kids: Vec<String> = (0..page_streams.len())
        .map(|i| format!("{} 0 R", 5 + 2 * i))
        .collect();
    let pages_body = format!(
        "<< /Type /Pages /Kids [{}] /Count {} >>",
        kids.join(" "),
        page_streams.len()
    );
    push_obj(&mut out, &mut offsets, 2, pages_body.as_bytes());

    push_obj(
        &mut out,
        &mut offsets,
        3,
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    );
    push_obj(
        &mut out,
        &mut offsets,
        4,
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
    );

    for (i, stream) in page_streams.iter().enumerate() {
        let page_id = 5 + 2 * i;
        let content_id = 6 + 2 * i;
        let page_body = format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {PAGE_W:.0} {PAGE_H:.0}] \
             /Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> /Contents {content_id} 0 R >>"
        );
        push_obj(&mut out, &mut offsets, page_id, page_body.as_bytes());

        let content_bytes = stream.as_bytes();
        let mut stream_body = format!("<< /Length {} >>\nstream\n", content_bytes.len()).into_bytes();
        stream_body.extend_from_slice(content_bytes);
        stream_body.extend_from_slice(b"\nendstream");
        push_obj(&mut out, &mut offsets, content_id, &stream_body);
    }

    let xref_offset = out.len();
    let total_objs = offsets.len() + 1;
    let _ = write_xref(&mut out, &offsets, total_objs, xref_offset);
    out
}

fn write_xref(
    out: &mut Vec<u8>,
    offsets: &[usize],
    total_objs: usize,
    xref_offset: usize,
) -> std::fmt::Result {
    let mut s = format!("xref\n0 {total_objs}\n0000000000 65535 f \n");
    for &off in offsets {
        writeln!(s, "{off:010} 00000 n ")?;
    }
    writeln!(s, "trailer\n<< /Size {total_objs} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF")?;
    out.extend_from_slice(s.as_bytes());
    Ok(())
}
