use std::f64::consts::PI;
use std::fmt::Write as _;

use admet_calc::MoleculeEvaluation;
use admet_calc::pharmacokinetics::{
    BOILED_EGG_WHITE_ANGLE_RAD, BOILED_EGG_WHITE_AXES, BOILED_EGG_WHITE_CENTER,
    BOILED_EGG_YOLK_ANGLE_RAD, BOILED_EGG_YOLK_AXES, BOILED_EGG_YOLK_CENTER,
};

/// Render a 6-axis Bioavailability Radar chart as standalone SVG.
pub fn render_radar_svg(eval: &MoleculeEvaluation) -> String {
    let w = 280.0;
    let h = 260.0;
    let cx = 140.0;
    let cy = 132.0;
    let max_r = 84.0;
    let labels = ["LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"];

    let mut out = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
    );

    // Optimal pink zone between r=0.25 and r=0.85
    let outer_pts = hex_points(cx, cy, max_r * 0.85);
    let inner_pts = hex_points(cx, cy, max_r * 0.25);
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"#fde2e4\" stroke=\"#f4a6b0\" stroke-width=\"1\"/>",
        fmt_pts(&outer_pts)
    );
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"#ffffff\" stroke=\"#f4a6b0\" stroke-width=\"1\"/>",
        fmt_pts(&inner_pts)
    );

    // Concentric hexagon grid
    for frac in [0.25, 0.5, 0.75, 1.0] {
        let pts = hex_points(cx, cy, max_r * frac);
        let _ = writeln!(
            out,
            "<polygon points=\"{}\" fill=\"none\" stroke=\"#d0d0d0\" stroke-width=\"0.8\"/>",
            fmt_pts(&pts)
        );
    }

    // Spokes & labels
    for (i, label) in labels.iter().enumerate() {
        let ang = -PI / 2.0 + (i as f64) * PI / 3.0;
        let x2 = cx + max_r * ang.cos();
        let y2 = cy + max_r * ang.sin();
        let lx = cx + (max_r + 18.0) * ang.cos();
        let ly = cy + (max_r + 18.0) * ang.sin() + 4.0;
        let _ = writeln!(
            out,
            "<line x1=\"{cx:.1}\" y1=\"{cy:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#cccccc\" stroke-width=\"0.8\"/>"
        );
        let _ = writeln!(
            out,
            "<text x=\"{lx:.1}\" y=\"{ly:.1}\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" font-weight=\"700\" text-anchor=\"middle\" fill=\"#444444\">{label}</text>"
        );
    }

    // Molecule polygon
    let mut mol_pts = Vec::with_capacity(6);
    for (i, &val) in eval.bioavailability_radar.normalized_axes.iter().enumerate() {
        let ang = -PI / 2.0 + (i as f64) * PI / 3.0;
        let r = max_r * val.clamp(0.05, 1.0);
        mol_pts.push((cx + r * ang.cos(), cy + r * ang.sin()));
    }
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"rgba(217, 72, 72, 0.25)\" stroke=\"#c92a2a\" stroke-width=\"2\"/>",
        fmt_pts(&mol_pts)
    );
    for (px, py) in &mol_pts {
        let _ = writeln!(
            out,
            "<circle cx=\"{px:.1}\" cy=\"{py:.1}\" r=\"3\" fill=\"#c92a2a\"/>"
        );
    }

    let _ = writeln!(out, "</svg>");
    out
}

/// Render an interactive-style BOILED-Egg SVG scatter plot for a batch of evaluations.
pub fn render_boiled_egg_svg(evals: &[MoleculeEvaluation]) -> String {
    let w = 640.0;
    let h = 440.0;
    let pad_l = 60.0;
    let pad_r = 160.0;
    let pad_t = 30.0;
    let pad_b = 50.0;
    let plot_w = w - pad_l - pad_r;
    let plot_h = h - pad_t - pad_b;

    let tpsa_min = 0.0;
    let tpsa_max = 200.0;
    let wlogp_min = -4.0;
    let wlogp_max = 8.0;

    let map_x = |tpsa: f64| pad_l + ((tpsa - tpsa_min) / (tpsa_max - tpsa_min)) * plot_w;
    let map_y = |wlogp: f64| pad_t + plot_h - ((wlogp - wlogp_min) / (wlogp_max - wlogp_min)) * plot_h;

    let mut out = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n\
         <rect x=\"{pad_l}\" y=\"{pad_t}\" width=\"{plot_w}\" height=\"{plot_h}\" fill=\"#d9d9d9\" stroke=\"#555555\" stroke-width=\"1\"/>\n"
    );

    // White ellipse polygon (sampled)
    let white_pts = sample_ellipse_svg(
        BOILED_EGG_WHITE_CENTER,
        BOILED_EGG_WHITE_AXES,
        BOILED_EGG_WHITE_ANGLE_RAD,
        &map_x,
        &map_y,
    );
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"#ffffff\" stroke=\"#888888\" stroke-width=\"1\"/>",
        fmt_pts(&white_pts)
    );

    // Yolk ellipse polygon (sampled)
    let yolk_pts = sample_ellipse_svg(
        BOILED_EGG_YOLK_CENTER,
        BOILED_EGG_YOLK_AXES,
        BOILED_EGG_YOLK_ANGLE_RAD,
        &map_x,
        &map_y,
    );
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"#f5c518\" stroke=\"#c99a00\" stroke-width=\"1\"/>",
        fmt_pts(&yolk_pts)
    );

    // Axes ticks
    for t in [0, 40, 80, 120, 160, 200] {
        let x = map_x(f64::from(t));
        let y0 = pad_t + plot_h;
        let _ = writeln!(
            out,
            "<line x1=\"{x:.1}\" y1=\"{y0:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\" stroke=\"#444444\"/>",
            y0 + 5.0
        );
        let _ = writeln!(
            out,
            "<text x=\"{x:.1}\" y=\"{:.1}\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" text-anchor=\"middle\" fill=\"#333333\">{t}</text>",
            y0 + 18.0
        );
    }
    for lp in [-4, -2, 0, 2, 4, 6, 8] {
        let y = map_y(f64::from(lp));
        let _ = writeln!(
            out,
            "<line x1=\"{:.1}\" y1=\"{y:.1}\" x2=\"{pad_l:.1}\" y2=\"{y:.1}\" stroke=\"#444444\"/>",
            pad_l - 5.0
        );
        let _ = writeln!(
            out,
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" text-anchor=\"end\" fill=\"#333333\">{lp}</text>",
            pad_l - 8.0,
            y + 3.0
        );
    }

    let _ = writeln!(
        out,
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Helvetica,sans-serif\" font-size=\"11\" font-weight=\"700\" text-anchor=\"middle\" fill=\"#222222\">TPSA (A^2)</text>",
        pad_l + plot_w * 0.5,
        h - 12.0
    );
    let _ = writeln!(
        out,
        "<text x=\"18\" y=\"{:.1}\" transform=\"rotate(-90 18 {:.1})\" font-family=\"Helvetica,sans-serif\" font-size=\"11\" font-weight=\"700\" text-anchor=\"middle\" fill=\"#222222\">WLOGP</text>",
        pad_t + plot_h * 0.5,
        pad_t + plot_h * 0.5
    );

    // Plot molecules
    let mut out_of_range = 0usize;
    for e in evals {
        if !e.boiled_egg.in_plot_bounds {
            out_of_range += 1;
            continue;
        }
        let px = map_x(e.boiled_egg.tpsa);
        let py = map_y(e.boiled_egg.wlogp);
        let color = if e.boiled_egg.pgp_substrate {
            "#1c7ed6"
        } else {
            "#e03131"
        };
        let _ = writeln!(
            out,
            "<circle cx=\"{px:.1}\" cy=\"{py:.1}\" r=\"5\" fill=\"{color}\" stroke=\"#111111\" stroke-width=\"0.9\"/>"
        );
        let safe_name = escape_xml(&e.name);
        let _ = writeln!(
            out,
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Helvetica,sans-serif\" font-size=\"9\" fill=\"#111111\">{safe_name}</text>",
            px + 7.0,
            py + 3.0
        );
    }

    // Legend
    let lx = pad_l + plot_w + 15.0;
    let _ = writeln!(
        out,
        "<rect x=\"{lx}\" y=\"40\" width=\"135\" height=\"150\" fill=\"#fbfbfb\" stroke=\"#cccccc\"/>\n\
         <rect x=\"{:.1}\" y=\"54\" width=\"12\" height=\"12\" fill=\"#ffffff\" stroke=\"#888888\"/>\n\
         <text x=\"{:.1}\" y=\"64\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" fill=\"#222\">HIA (White)</text>\n\
         <rect x=\"{:.1}\" y=\"74\" width=\"12\" height=\"12\" fill=\"#f5c518\" stroke=\"#c99a00\"/>\n\
         <text x=\"{:.1}\" y=\"84\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" fill=\"#222\">BBB (Yolk)</text>\n\
         <circle cx=\"{:.1}\" cy=\"104\" r=\"4.5\" fill=\"#1c7ed6\" stroke=\"#111\"/>\n\
         <text x=\"{:.1}\" y=\"108\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" fill=\"#222\">PGP+</text>\n\
         <circle cx=\"{:.1}\" cy=\"124\" r=\"4.5\" fill=\"#e03131\" stroke=\"#111\"/>\n\
         <text x=\"{:.1}\" y=\"128\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" fill=\"#222\">PGP-</text>\n\
         <text x=\"{:.1}\" y=\"158\" font-family=\"Helvetica,sans-serif\" font-size=\"9\" fill=\"#555\">Out of range: {out_of_range}</text>",
        lx + 10.0,
        lx + 28.0,
        lx + 10.0,
        lx + 28.0,
        lx + 16.0,
        lx + 28.0,
        lx + 16.0,
        lx + 28.0,
        lx + 10.0
    );

    let _ = writeln!(out, "</svg>");
    out
}

fn hex_points(cx: f64, cy: f64, r: f64) -> Vec<(f64, f64)> {
    (0..6)
        .map(|i| {
            let ang = -PI / 2.0 + (i as f64) * PI / 3.0;
            (cx + r * ang.cos(), cy + r * ang.sin())
        })
        .collect()
}

fn sample_ellipse_svg(
    center: (f64, f64),
    axes: (f64, f64),
    angle: f64,
    map_x: &dyn Fn(f64) -> f64,
    map_y: &dyn Fn(f64) -> f64,
) -> Vec<(f64, f64)> {
    let steps = 72;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    (0..steps)
        .map(|i| {
            let t = 2.0 * PI * (i as f64) / (steps as f64);
            let lx = axes.0 * t.cos();
            let ly = axes.1 * t.sin();
            let gx = center.0 + lx * cos_a - ly * sin_a;
            let gy = center.1 + lx * sin_a + ly * cos_a;
            (map_x(gx).clamp(60.0, 480.0), map_y(gy).clamp(30.0, 390.0))
        })
        .collect()
}

fn fmt_pts(pts: &[(f64, f64)]) -> String {
    pts.iter()
        .map(|(x, y)| format!("{x:.1},{y:.1}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
