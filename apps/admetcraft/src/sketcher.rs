use std::f64::consts::PI;

use admet_calc::{MetabolismRisk, MetabolismSite};
use chem_core::{Atom, BondOrder, BondStereo, Molecule, ensure_2d_coords, parse_smiles};
use egui::{Color32, Pos2, Rect, Sense, Stroke, Ui, Vec2};

/// Active drawing tool in the 2D molecular sketcher.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SketchTool {
    Atom,
    SingleBond,
    DoubleBond,
    TripleBond,
    BenzeneRing,
    CyclohexaneRing,
    CyclopentaneRing,
    Erase,
}

/// Interactive 2D chemical structure editor synchronized with SMILES.
#[derive(Debug, Clone)]
pub struct MoleculeSketcher {
    pub mol: Molecule,
    pub tool: SketchTool,
    pub element: &'static str,
    pub selected_atom: Option<usize>,
    pub drag_from_atom: Option<usize>,
    pub cached_smiles: String,
}

impl Default for MoleculeSketcher {
    fn default() -> Self {
        let mut s = Self {
            mol: Molecule::new(),
            tool: SketchTool::SingleBond,
            element: "C",
            selected_atom: None,
            drag_from_atom: None,
            cached_smiles: String::new(),
        };
        let _ = s.load_smiles("c1ccccc1");
        s
    }
}

impl MoleculeSketcher {
    /// Load a SMILES string into the 2D sketcher canvas.
    pub fn load_smiles(&mut self, smiles: &str) -> Result<(), chem_core::ChemError> {
        let mut parsed = parse_smiles(smiles)?;
        ensure_2d_coords(&mut parsed);
        self.mol = parsed;
        self.selected_atom = None;
        self.drag_from_atom = None;
        self.refresh_smiles();
        Ok(())
    }

    /// Recompute perception & canonical SMILES after an interactive edit.
    pub fn refresh_smiles(&mut self) {
        let _ = self.mol.perceive();
        self.cached_smiles = self.mol.to_canonical_smiles();
    }

    /// Clear all atoms and bonds from the sketcher.
    pub fn clear(&mut self) {
        self.mol = Molecule::new();
        self.selected_atom = None;
        self.drag_from_atom = None;
        self.cached_smiles.clear();
    }

    /// Render the sketcher toolbar + interactive canvas, returning `true` if the graph changed.
    pub fn ui(&mut self, ui: &mut Ui, dark_mode: bool) -> bool {
        let mut changed = false;

        ui.horizontal_wrapped(|ui| {
            let tools = [
                (SketchTool::SingleBond, "-"),
                (SketchTool::DoubleBond, "="),
                (SketchTool::TripleBond, "==="),
                (SketchTool::BenzeneRing, "Benz"),
                (SketchTool::CyclohexaneRing, "Hex"),
                (SketchTool::CyclopentaneRing, "Pent"),
                (SketchTool::Atom, "Atom"),
                (SketchTool::Erase, "Del"),
            ];
            for (t, label) in tools {
                if ui.selectable_label(self.tool == t, label).clicked() {
                    self.tool = t;
                }
            }
        });

        ui.horizontal_wrapped(|ui| {
            for sym in ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"] {
                if ui.selectable_label(self.element == sym, sym).clicked() {
                    self.element = sym;
                    if self.tool == SketchTool::Erase {
                        self.tool = SketchTool::Atom;
                    }
                }
            }
        });

        let avail = Vec2::new(ui.available_width(), 220.0);
        let (rect, response) = ui.allocate_exact_size(avail, Sense::click_and_drag());
        let bg = if dark_mode {
            Color32::from_rgb(30, 32, 36)
        } else {
            Color32::WHITE
        };
        let border = if dark_mode {
            Color32::from_rgb(70, 74, 82)
        } else {
            Color32::from_rgb(200, 200, 200)
        };
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 4.0, bg);
        painter.rect_stroke(
            rect,
            4.0,
            Stroke::new(1.0_f32, border),
            egui::StrokeKind::Outside,
        );

        let view = CanvasTransform::fit_molecule(&self.mol, rect);

        // Handle pointer interactions
        if response.clicked()
            && let Some(pos) = response.interact_pointer_pos()
        {
            changed |= self.handle_click(pos, &view);
        } else if response.drag_started()
            && let Some(pos) = response.interact_pointer_pos()
        {
            self.drag_from_atom = self.hit_test_atom(pos, &view, 14.0);
        } else if response.drag_stopped()
            && let Some(pos) = response.interact_pointer_pos()
        {
            changed |= self.handle_drag_release(pos, &view);
        }

        // Draw molecule
        draw_molecule_in_rect(&painter, &self.mol, rect, dark_mode, self.selected_atom);

        // Drag preview line
        if response.dragged()
            && let (Some(src), Some(cur_pos)) =
                (self.drag_from_atom, response.interact_pointer_pos())
            && let Some(c) = self.mol.atoms.get(src).and_then(|a| a.coords_2d)
        {
            let p0 = view.to_screen(c);
            painter.line_segment(
                [p0, cur_pos],
                Stroke::new(1.5_f32, Color32::from_rgb(220, 90, 90)),
            );
        }

        if changed {
            self.refresh_smiles();
        }
        changed
    }

    fn hit_test_atom(&self, pos: Pos2, view: &CanvasTransform, radius: f32) -> Option<usize> {
        let mut best: Option<(usize, f32)> = None;
        for (idx, a) in self.mol.atoms.iter().enumerate() {
            let Some(c) = a.coords_2d else {
                continue;
            };
            let sp = view.to_screen(c);
            let d = sp.distance(pos);
            if d <= radius && best.is_none_or(|(_, bd)| d < bd) {
                best = Some((idx, d));
            }
        }
        best.map(|(i, _)| i)
    }

    fn hit_test_bond(&self, pos: Pos2, view: &CanvasTransform, tol: f32) -> Option<usize> {
        for b in &self.mol.bonds {
            let (Some(ca), Some(cb)) = (
                self.mol.atoms.get(b.u).and_then(|a| a.coords_2d),
                self.mol.atoms.get(b.v).and_then(|a| a.coords_2d),
            ) else {
                continue;
            };
            let p1 = view.to_screen(ca);
            let p2 = view.to_screen(cb);
            if dist_point_to_segment(pos, p1, p2) <= tol {
                return Some(b.id);
            }
        }
        None
    }

    fn handle_click(&mut self, pos: Pos2, view: &CanvasTransform) -> bool {
        let hit_atom = self.hit_test_atom(pos, view, 14.0);
        let hit_bond = if hit_atom.is_none() {
            self.hit_test_bond(pos, view, 8.0)
        } else {
            None
        };

        match self.tool {
            SketchTool::Erase => {
                if let Some(aid) = hit_atom {
                    self.remove_atom(aid);
                    return true;
                }
                if let Some(bid) = hit_bond {
                    self.remove_bond(bid);
                    return true;
                }
                false
            }
            SketchTool::Atom => {
                if let Some(aid) = hit_atom {
                    self.set_atom_element(aid, self.element);
                    return true;
                }
                let world = view.to_world(pos);
                self.add_free_atom(self.element, world);
                true
            }
            SketchTool::SingleBond | SketchTool::DoubleBond | SketchTool::TripleBond => {
                let order = match self.tool {
                    SketchTool::DoubleBond => BondOrder::Double,
                    SketchTool::TripleBond => BondOrder::Triple,
                    _ => BondOrder::Single,
                };
                if let Some(bid) = hit_bond {
                    if let Some(b) = self.mol.bonds.get_mut(bid) {
                        b.order = match b.order {
                            BondOrder::Single => BondOrder::Double,
                            BondOrder::Double => BondOrder::Triple,
                            BondOrder::Triple => BondOrder::Single,
                            BondOrder::Aromatic => BondOrder::Single,
                        };
                    }
                    return true;
                }
                if let Some(aid) = hit_atom {
                    if let Some(prev) = self.selected_atom
                        && prev != aid
                    {
                        let _ = self.mol.add_bond(prev, aid, order, BondStereo::None);
                        self.selected_atom = None;
                        return true;
                    }
                    // Extend a new substituent from clicked atom at ~1.5 distance
                    let new_pos = self.suggest_extension_pos(aid);
                    let new_id = self.add_free_atom(self.element, new_pos);
                    let _ = self.mol.add_bond(aid, new_id, order, BondStereo::None);
                    self.selected_atom = None;
                    return true;
                }
                // Empty canvas click: place atom (or bond from selected)
                let world = view.to_world(pos);
                let new_id = self.add_free_atom(self.element, world);
                if let Some(prev) = self.selected_atom {
                    let _ = self.mol.add_bond(prev, new_id, order, BondStereo::None);
                }
                self.selected_atom = Some(new_id);
                true
            }
            SketchTool::BenzeneRing
            | SketchTool::CyclohexaneRing
            | SketchTool::CyclopentaneRing => {
                let sides = match self.tool {
                    SketchTool::CyclopentaneRing => 5,
                    _ => 6,
                };
                let aromatic = self.tool == SketchTool::BenzeneRing;
                if let Some(aid) = hit_atom {
                    self.attach_ring_to_atom(aid, sides, aromatic);
                } else {
                    let world = view.to_world(pos);
                    self.place_ring_at(world, sides, aromatic);
                }
                true
            }
        }
    }

    fn handle_drag_release(&mut self, pos: Pos2, view: &CanvasTransform) -> bool {
        let Some(src) = self.drag_from_atom.take() else {
            return false;
        };
        let order = match self.tool {
            SketchTool::DoubleBond => BondOrder::Double,
            SketchTool::TripleBond => BondOrder::Triple,
            _ => BondOrder::Single,
        };
        if let Some(dst) = self.hit_test_atom(pos, view, 14.0) {
            if dst != src {
                if let Some(existing) = self
                    .mol
                    .bonds
                    .iter_mut()
                    .find(|b| (b.u == src && b.v == dst) || (b.u == dst && b.v == src))
                {
                    existing.order = order;
                } else {
                    let _ = self.mol.add_bond(src, dst, order, BondStereo::None);
                }
                return true;
            }
            return false;
        }
        let world = view.to_world(pos);
        let new_id = self.add_free_atom(self.element, world);
        let _ = self.mol.add_bond(src, new_id, order, BondStereo::None);
        true
    }

    fn add_free_atom(&mut self, sym: &str, world: [f64; 2]) -> usize {
        if let Ok(mut atom) = Atom::new(0, sym) {
            atom.coords_2d = Some(world);
            self.mol.add_atom(atom)
        } else {
            0
        }
    }

    fn set_atom_element(&mut self, aid: usize, sym: &str) {
        if let Ok(new_atom) = Atom::new(aid, sym)
            && let Some(a) = self.mol.atoms.get_mut(aid)
        {
            a.atomic_number = new_atom.atomic_number;
            a.symbol = new_atom.symbol;
            a.aromatic = false;
        }
    }

    fn remove_atom(&mut self, target: usize) {
        if target >= self.mol.atoms.len() {
            return;
        }
        let mut rebuilt = Molecule::new();
        let mut map = vec![None; self.mol.atoms.len()];
        for (i, a) in self.mol.atoms.iter().enumerate() {
            if i != target {
                let mut cloned = a.clone();
                cloned.aromatic = false;
                let nid = rebuilt.add_atom(cloned);
                map[i] = Some(nid);
            }
        }
        for b in &self.mol.bonds {
            if let (Some(nu), Some(nv)) = (map[b.u], map[b.v]) {
                let _ = rebuilt.add_bond(nu, nv, b.order, b.stereo);
            }
        }
        self.mol = rebuilt;
        self.selected_atom = None;
    }

    fn remove_bond(&mut self, target_bid: usize) {
        let mut rebuilt = Molecule::new();
        for a in &self.mol.atoms {
            let mut cloned = a.clone();
            cloned.aromatic = false;
            rebuilt.add_atom(cloned);
        }
        for b in &self.mol.bonds {
            if b.id != target_bid {
                let _ = rebuilt.add_bond(b.u, b.v, b.order, b.stereo);
            }
        }
        self.mol = rebuilt;
    }

    fn suggest_extension_pos(&self, aid: usize) -> [f64; 2] {
        let base = self.mol.atoms[aid].coords_2d.unwrap_or([0.0, 0.0]);
        let nbrs: Vec<[f64; 2]> = self
            .mol
            .neighbors(aid)
            .filter_map(|n| self.mol.atoms[n].coords_2d)
            .collect();
        let ang = if nbrs.is_empty() {
            0.0
        } else if nbrs.len() == 1 {
            let dx = base[0] - nbrs[0][0];
            let dy = base[1] - nbrs[0][1];
            dy.atan2(dx) + PI / 3.0
        } else {
            let sx: f64 = nbrs.iter().map(|c| c[0] - base[0]).sum();
            let sy: f64 = nbrs.iter().map(|c| c[1] - base[1]).sum();
            (-sy).atan2(-sx)
        };
        [base[0] + 1.5 * ang.cos(), base[1] + 1.5 * ang.sin()]
    }

    fn place_ring_at(&mut self, center: [f64; 2], sides: usize, aromatic: bool) {
        let r = 1.5 / (2.0 * (PI / sides as f64).sin());
        let mut ids = Vec::with_capacity(sides);
        for i in 0..sides {
            let ang = (i as f64) * 2.0 * PI / (sides as f64);
            let id =
                self.add_free_atom("C", [center[0] + r * ang.cos(), center[1] + r * ang.sin()]);
            ids.push(id);
        }
        for i in 0..sides {
            let u = ids[i];
            let v = ids[(i + 1) % sides];
            let order = if aromatic && i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            let _ = self.mol.add_bond(u, v, order, BondStereo::None);
        }
    }

    fn attach_ring_to_atom(&mut self, aid: usize, sides: usize, aromatic: bool) {
        let ext = self.suggest_extension_pos(aid);
        let base = self.mol.atoms[aid].coords_2d.unwrap_or([0.0, 0.0]);
        let dir = (ext[1] - base[1]).atan2(ext[0] - base[0]);
        let r = 1.5 / (2.0 * (PI / sides as f64).sin());
        let center = [ext[0] + r * dir.cos(), ext[1] + r * dir.sin()];
        let first_id = self.add_free_atom("C", ext);
        let _ = self
            .mol
            .add_bond(aid, first_id, BondOrder::Single, BondStereo::None);
        let mut ids = vec![first_id];
        let base_ang = (ext[1] - center[1]).atan2(ext[0] - center[0]);
        for i in 1..sides {
            let ang = base_ang + (i as f64) * 2.0 * PI / (sides as f64);
            ids.push(
                self.add_free_atom("C", [center[0] + r * ang.cos(), center[1] + r * ang.sin()]),
            );
        }
        for i in 0..sides {
            let order = if aromatic && i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            let _ = self
                .mol
                .add_bond(ids[i], ids[(i + 1) % sides], order, BondStereo::None);
        }
    }
}

struct CanvasTransform {
    cx: f32,
    cy: f32,
    mid_x: f64,
    mid_y: f64,
    scale: f64,
}

impl CanvasTransform {
    fn fit_molecule(mol: &Molecule, rect: Rect) -> Self {
        let coords: Vec<[f64; 2]> = mol.atoms.iter().filter_map(|a| a.coords_2d).collect();
        let (mid_x, mid_y, scale) = if coords.is_empty() {
            (0.0, 0.0, 26.0)
        } else {
            let min_x = coords.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
            let max_x = coords
                .iter()
                .map(|c| c[0])
                .fold(f64::NEG_INFINITY, f64::max);
            let min_y = coords.iter().map(|c| c[1]).fold(f64::INFINITY, f64::min);
            let max_y = coords
                .iter()
                .map(|c| c[1])
                .fold(f64::NEG_INFINITY, f64::max);
            let span_x = (max_x - min_x).max(2.0);
            let span_y = (max_y - min_y).max(2.0);
            let sc = ((f64::from(rect.width()) * 0.72) / span_x)
                .min((f64::from(rect.height()) * 0.72) / span_y)
                .clamp(14.0, 34.0);
            (0.5 * (min_x + max_x), 0.5 * (min_y + max_y), sc)
        };
        Self {
            cx: rect.center().x,
            cy: rect.center().y,
            mid_x,
            mid_y,
            scale,
        }
    }

    fn to_screen(&self, world: [f64; 2]) -> Pos2 {
        Pos2::new(
            self.cx + ((world[0] - self.mid_x) * self.scale) as f32,
            self.cy - ((world[1] - self.mid_y) * self.scale) as f32,
        )
    }

    fn to_world(&self, screen: Pos2) -> [f64; 2] {
        [
            self.mid_x + f64::from(screen.x - self.cx) / self.scale,
            self.mid_y - f64::from(screen.y - self.cy) / self.scale,
        ]
    }
}

/// Render a pre-parsed, pre-laid-out `Molecule` inside `rect` (zero per-frame parsing).
pub fn draw_cached_molecule_2d(
    painter: &egui::Painter,
    mol: &Molecule,
    rect: Rect,
    dark_mode: bool,
) {
    draw_molecule_in_rect(painter, mol, rect, dark_mode, None);
}

/// Render a pre-laid-out `Molecule` with translucent atomic Sites of Metabolism (SoM) halos.
pub fn draw_cached_molecule_2d_with_som(
    painter: &egui::Painter,
    mol: &Molecule,
    rect: Rect,
    dark_mode: bool,
    som_sites: &[MetabolismSite],
) {
    let view = CanvasTransform::fit_molecule(mol, rect);
    for site in som_sites {
        let Some(c) = mol.atoms.get(site.atom_index).and_then(|a| a.coords_2d) else {
            continue;
        };
        let p = view.to_screen(c);
        let (fill, stroke, r) = match site.risk {
            MetabolismRisk::High => (
                Color32::from_rgba_unmultiplied(220, 45, 45, 78),
                Color32::from_rgb(205, 40, 40),
                11.5_f32,
            ),
            MetabolismRisk::Moderate => (
                Color32::from_rgba_unmultiplied(230, 145, 25, 68),
                Color32::from_rgb(205, 120, 15),
                9.5_f32,
            ),
            MetabolismRisk::Low => (
                Color32::from_rgba_unmultiplied(180, 180, 50, 45),
                Color32::from_rgb(160, 150, 40),
                8.0_f32,
            ),
        };
        painter.circle_filled(p, r, fill);
        painter.circle_stroke(p, r, Stroke::new(1.2_f32, stroke));
    }
    draw_molecule_in_rect(painter, mol, rect, dark_mode, None);
}

/// Render a 2D chemical structure inside `rect` (fallback convenience wrapper).
#[allow(dead_code)]
pub fn draw_molecule_2d(painter: &egui::Painter, smiles: &str, rect: Rect, dark_mode: bool) {
    let Ok(mut mol) = parse_smiles(smiles) else {
        return;
    };
    ensure_2d_coords(&mut mol);
    draw_molecule_in_rect(painter, &mol, rect, dark_mode, None);
}

fn draw_molecule_in_rect(
    painter: &egui::Painter,
    mol: &Molecule,
    rect: Rect,
    dark_mode: bool,
    selected_atom: Option<usize>,
) {
    let view = CanvasTransform::fit_molecule(mol, rect);
    let bond_color = if dark_mode {
        Color32::from_rgb(220, 224, 230)
    } else {
        Color32::from_rgb(35, 35, 40)
    };
    let bg_color = if dark_mode {
        Color32::from_rgb(30, 32, 36)
    } else {
        Color32::WHITE
    };

    for b in &mol.bonds {
        let (Some(ca), Some(cb)) = (
            mol.atoms.get(b.u).and_then(|a| a.coords_2d),
            mol.atoms.get(b.v).and_then(|a| a.coords_2d),
        ) else {
            continue;
        };
        let p1 = view.to_screen(ca);
        let p2 = view.to_screen(cb);
        let dir = p2 - p1;
        let len = dir.length().max(1e-3);
        let n = Vec2::new(-dir.y / len, dir.x / len) * 2.5;
        match b.order {
            BondOrder::Single | BondOrder::Aromatic => {
                painter.line_segment([p1, p2], Stroke::new(1.6_f32, bond_color));
            }
            BondOrder::Double => {
                painter.line_segment([p1 + n, p2 + n], Stroke::new(1.4_f32, bond_color));
                painter.line_segment([p1 - n, p2 - n], Stroke::new(1.4_f32, bond_color));
            }
            BondOrder::Triple => {
                painter.line_segment([p1, p2], Stroke::new(1.3_f32, bond_color));
                painter.line_segment(
                    [p1 + n * 1.4, p2 + n * 1.4],
                    Stroke::new(1.2_f32, bond_color),
                );
                painter.line_segment(
                    [p1 - n * 1.4, p2 - n * 1.4],
                    Stroke::new(1.2_f32, bond_color),
                );
            }
        }
    }

    // Draw clean inscribed aromatic circle for fully aromatic rings
    let arom_circle_col = if dark_mode {
        Color32::from_rgb(175, 180, 190)
    } else {
        Color32::from_rgb(65, 65, 75)
    };
    for ring in &mol.rings {
        if ring.len() >= 5 && ring.iter().all(|&aid| mol.atoms[aid].aromatic) {
            let pts: Vec<Pos2> = ring
                .iter()
                .filter_map(|&aid| mol.atoms[aid].coords_2d.map(|c| view.to_screen(c)))
                .collect();
            if pts.len() == ring.len() {
                let rcx = pts.iter().map(|p| p.x).sum::<f32>() / pts.len() as f32;
                let rcy = pts.iter().map(|p| p.y).sum::<f32>() / pts.len() as f32;
                let center = Pos2::new(rcx, rcy);
                let avg_r = pts.iter().map(|p| p.distance(center)).sum::<f32>() / pts.len() as f32;
                painter.circle_stroke(center, avg_r * 0.56, Stroke::new(1.1_f32, arom_circle_col));
            }
        }
    }

    for (idx, a) in mol.atoms.iter().enumerate() {
        let Some(c) = a.coords_2d else {
            continue;
        };
        let p = view.to_screen(c);
        if selected_atom == Some(idx) {
            painter.circle_stroke(p, 9.0, Stroke::new(1.6_f32, Color32::from_rgb(220, 70, 70)));
        }
        if a.atomic_number == 6 && mol.degree(idx) > 0 {
            continue;
        }
        let col = match a.atomic_number {
            7 => Color32::from_rgb(45, 110, 230),
            8 => Color32::from_rgb(215, 45, 45),
            16 => Color32::from_rgb(190, 145, 15),
            9 | 17 => Color32::from_rgb(35, 165, 65),
            35 => Color32::from_rgb(160, 50, 35),
            53 => Color32::from_rgb(130, 40, 165),
            _ => bond_color,
        };
        painter.circle_filled(p, 7.5, bg_color);
        painter.text(
            p,
            egui::Align2::CENTER_CENTER,
            &a.symbol,
            egui::FontId::proportional(11.5),
            col,
        );
    }
}

fn dist_point_to_segment(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let len2 = ab.length_sq();
    if len2 < 1e-4 {
        return p.distance(a);
    }
    let t = ((p - a).dot(ab) / len2).clamp(0.0, 1.0);
    let proj = a + ab * t;
    p.distance(proj)
}
