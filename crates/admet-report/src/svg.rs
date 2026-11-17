use std::fmt::Write as _;

use admet_calc::pharmacokinetics::{
    BOILED_EGG_WHITE_ANGLE_RAD, BOILED_EGG_WHITE_AXES, BOILED_EGG_WHITE_CENTER,
    BOILED_EGG_YOLK_ANGLE_RAD, BOILED_EGG_YOLK_AXES, BOILED_EGG_YOLK_CENTER,
};
use admet_calc::{BioavailabilityRadar, BoiledEggGeometry, MoleculeEvaluation};

/// Render a 6-axis Bioavailability Radar chart as standalone SVG using shared radar geometry.
pub fn render_radar_svg(eval: &MoleculeEvaluation) -> String {
    let w = 280.0;
    let h = 260.0;
    let cx = 140.0;
    let cy = 132.0;
    let max_r = 84.0;

    let mut out = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
    );

    let to_svg = |(ux, uy): (f64, f64)| (cx + ux * max_r, cy - uy * max_r);

    let outer_pts: Vec<(f64, f64)> = BioavailabilityRadar::unit_hexagon(0.85)
        .into_iter()
        .map(to_svg)
        .collect();
    let inner_pts: Vec<(f64, f64)> = BioavailabilityRadar::unit_hexagon(0.25)
        .into_iter()
        .map(to_svg)
        .collect();
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

    for frac in [0.25, 0.5, 0.75, 1.0] {
        let pts: Vec<(f64, f64)> = BioavailabilityRadar::unit_hexagon(frac)
            .into_iter()
            .map(to_svg)
            .collect();
        let _ = writeln!(
            out,
            "<polygon points=\"{}\" fill=\"none\" stroke=\"#d0d0d0\" stroke-width=\"0.8\"/>",
            fmt_pts(&pts)
        );
    }

    let full_hex = BioavailabilityRadar::unit_hexagon(1.0);
    let label_hex = BioavailabilityRadar::unit_hexagon((max_r + 18.0) / max_r);
    for (i, label) in BioavailabilityRadar::AXIS_LABELS.iter().enumerate() {
        let (x2, y2) = to_svg(full_hex[i]);
        let (lx, ly) = to_svg(label_hex[i]);
        let _ = writeln!(
            out,
            "<line x1=\"{cx:.1}\" y1=\"{cy:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#cccccc\" stroke-width=\"0.8\"/>"
        );
        let _ = writeln!(
            out,
            "<text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"Helvetica,sans-serif\" font-size=\"10\" font-weight=\"700\" text-anchor=\"middle\" fill=\"#444444\">{label}</text>",
            ly + 4.0
        );
    }

    let mol_pts: Vec<(f64, f64)> = eval
        .bioavailability_radar
        .unit_polygon()
        .into_iter()
        .map(to_svg)
        .collect();
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

/// Render a BOILED-Egg SVG scatter plot using shared `BoiledEggGeometry`.
pub fn render_boiled_egg_svg(evals: &[MoleculeEvaluation]) -> String {
    let w = 640.0;
    let h = 440.0;
    let pad_l = 60.0;
    let pad_r = 160.0;
    let pad_t = 30.0;
    let pad_b = 50.0;
    let plot_w = w - pad_l - pad_r;
    let plot_h = h - pad_t - pad_b;

    let map_x = |tpsa: f64| (pad_l + (tpsa / 200.0) * plot_w).clamp(pad_l, pad_l + plot_w);
    let map_y = |wlogp: f64| {
        (pad_t + plot_h - ((wlogp + 4.0) / 12.0) * plot_h).clamp(pad_t, pad_t + plot_h)
    };

    let mut out = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n\
         <rect x=\"{pad_l}\" y=\"{pad_t}\" width=\"{plot_w}\" height=\"{plot_h}\" fill=\"#d9d9d9\" stroke=\"#555555\" stroke-width=\"1\"/>\n"
    );

    let white_pts: Vec<(f64, f64)> = BoiledEggGeometry::sample_ellipse(
        BOILED_EGG_WHITE_CENTER,
        BOILED_EGG_WHITE_AXES,
        BOILED_EGG_WHITE_ANGLE_RAD,
        72,
    )
    .into_iter()
    .map(|(x, y)| (map_x(x), map_y(y)))
    .collect();
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"#ffffff\" stroke=\"#888888\" stroke-width=\"1\"/>",
        fmt_pts(&white_pts)
    );

    let yolk_pts: Vec<(f64, f64)> = BoiledEggGeometry::sample_ellipse(
        BOILED_EGG_YOLK_CENTER,
        BOILED_EGG_YOLK_AXES,
        BOILED_EGG_YOLK_ANGLE_RAD,
        72,
    )
    .into_iter()
    .map(|(x, y)| (map_x(x), map_y(y)))
    .collect();
    let _ = writeln!(
        out,
        "<polygon points=\"{}\" fill=\"#f5c518\" stroke=\"#c99a00\" stroke-width=\"1\"/>",
        fmt_pts(&yolk_pts)
    );

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

fn fmt_pts(pts: &[(f64, f64)]) -> String {
    pts.iter()
        .map(|(x, y)| format!("{x:.1},{y:.1}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn escape_xml(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        if c.is_control() {
            out.push(' ');
            continue;
        }
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}
