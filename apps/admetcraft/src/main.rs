#![forbid(unsafe_code)]

mod charts;
mod sketcher;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use admet_calc::MoleculeEvaluation;
use admet_engine::{
    AdmetSession, Command, CommandRegistry, FilterCriteria, MAX_READ_BYTES, PathSandbox,
};
use admet_report::{ExportFormat, write_report_to_writer};
use chem_core::{
    Molecule, ensure_2d_coords, parse_csv_batch, parse_sdf, parse_smiles, parse_smiles_batch,
};
use eframe::egui;
use egui::{Color32, RichText, Sense, Vec2};

use charts::{draw_bioavailability_radar, draw_boiled_egg_interactive};
use sketcher::{MoleculeSketcher, draw_cached_molecule_2d};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MainTab {
    Cards,
    BoiledEgg,
    Table,
}

struct AdmetCraftApp {
    session: AdmetSession,
    mol_cache: HashMap<String, Molecule>,
    sketcher: MoleculeSketcher,
    active_tab: MainTab,
    dark_mode: bool,
    show_egg_names: bool,
    palette_open: bool,
    palette_query: String,
    import_path: String,
    export_dir: String,
    status_msg: String,
    async_status: Arc<Mutex<Option<String>>>,
}

impl Default for AdmetCraftApp {
    fn default() -> Self {
        let mut session = AdmetSession::new();
        let status_msg = match CommandRegistry::execute(&mut session, Command::LoadExampleBatch) {
            Ok(out) => out.message,
            Err(e) => format!("Init warning: {e}"),
        };
        let export_dir = std::env::current_dir()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|_| ".".to_string());
        let mut sketcher = MoleculeSketcher::default();
        if let Some(first) = session.evaluations.first() {
            let _ = sketcher.load_smiles(&first.canonical_smiles);
        }
        let mut app = Self {
            session,
            mol_cache: HashMap::new(),
            sketcher,
            active_tab: MainTab::Cards,
            dark_mode: false,
            show_egg_names: true,
            palette_open: false,
            palette_query: String::new(),
            import_path: String::new(),
            export_dir,
            status_msg,
            async_status: Arc::new(Mutex::new(None)),
        };
        app.rebuild_mol_cache();
        app
    }
}

impl AdmetCraftApp {
    fn rebuild_mol_cache(&mut self) {
        self.mol_cache.clear();
        for ev in &self.session.evaluations {
            if self.mol_cache.contains_key(&ev.canonical_smiles) {
                continue;
            }
            if let Ok(mut mol) = parse_smiles(&ev.canonical_smiles) {
                ensure_2d_coords(&mut mol);
                self.mol_cache.insert(ev.canonical_smiles.clone(), mol);
            }
        }
    }

    fn run_command(&mut self, cmd: Command) {
        let is_history_nav = matches!(cmd, Command::Undo | Command::Redo);
        match CommandRegistry::execute(&mut self.session, cmd) {
            Ok(out) => {
                self.status_msg = out.message;
                self.rebuild_mol_cache();
                if is_history_nav {
                    let selected_eval = self
                        .session
                        .selected_index
                        .and_then(|idx| self.session.evaluations.get(idx))
                        .or_else(|| self.session.evaluations.first());
                    if let Some(ev) = selected_eval {
                        let _ = self.sketcher.load_smiles(&ev.canonical_smiles);
                    } else {
                        self.sketcher.clear();
                    }
                }
            }
            Err(e) => self.status_msg = format!("Error: {e}"),
        }
    }

    /// Export asynchronously on a background thread using `PathSandbox::open_write` so the UI
    /// thread never freezes on large PDFs and symlink/traversal escapes are blocked.
    fn export(&mut self, format: ExportFormat, filename: &str) {
        if self.session.evaluations.is_empty() {
            self.status_msg = "No molecules to export".to_string();
            return;
        }
        let export_root = PathBuf::from(&self.export_dir);
        let sandbox = match PathSandbox::new(None, Some(&export_root)) {
            Ok(sb) => sb,
            Err(e) => {
                self.status_msg = format!("Export directory error: {e}");
                return;
            }
        };
        let target_name = PathBuf::from(filename);
        let evals = self.session.evaluations.clone();
        let status_slot = Arc::clone(&self.async_status);
        self.status_msg = format!("Exporting {filename} in background...");

        std::thread::spawn(move || {
            let msg = match sandbox.open_write(&target_name) {
                Ok((resolved, mut file)) => match write_report_to_writer(&evals, &mut file, format)
                {
                    Ok(()) => format!("Exported {}", resolved.display()),
                    Err(e) => format!("Export error: {e}"),
                },
                Err(e) => format!("Sandbox export error: {e}"),
            };
            if let Ok(mut guard) = status_slot.lock() {
                *guard = Some(msg);
            }
        });
    }

    /// Import a `.smi`, `.csv`, or `.sdf` file inside `export_dir` (or the file's parent for
    /// explicit drag-and-drop) using `PathSandbox::read_to_string_limited` (rejects symlinks/FIFOs).
    fn import_file_path(&mut self, path: &Path, from_drop: bool) {
        let root = if from_drop {
            path.parent().unwrap_or_else(|| Path::new("."))
        } else {
            Path::new(&self.export_dir)
        };
        let sandbox = match PathSandbox::new(Some(root), None) {
            Ok(sb) => sb,
            Err(e) => {
                self.status_msg = format!("Import root error: {e}");
                return;
            }
        };
        let rel_or_abs = if from_drop {
            path.file_name().map(Path::new).unwrap_or(path)
        } else {
            path
        };
        let (resolved, content) = match sandbox.read_to_string_limited(rel_or_abs, MAX_READ_BYTES) {
            Ok(res) => res,
            Err(e) => {
                self.status_msg = format!("Sandbox import error: {e}");
                return;
            }
        };
        let ext = resolved
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        let cmd = if ext == "sdf" || ext == "mol" {
            Command::ImportSdf { content }
        } else if ext == "csv" {
            Command::ImportCsv { content }
        } else {
            Command::EvaluateBatchText { input: content }
        };
        self.run_command(cmd);
    }

    fn import_dropped_files(&mut self, paths: &[PathBuf]) {
        if paths.is_empty() {
            return;
        }
        if paths.len() == 1 {
            self.import_file_path(&paths[0], true);
            return;
        }

        let mut combined_lines = Vec::new();
        let mut loaded_files = 0usize;
        for path in paths {
            let root = path.parent().unwrap_or_else(|| Path::new("."));
            let sandbox = match PathSandbox::new(Some(root), None) {
                Ok(sb) => sb,
                Err(e) => {
                    self.status_msg = format!("Import root error: {e}");
                    return;
                }
            };
            let rel = path.file_name().map(Path::new).unwrap_or(path.as_path());
            let (resolved, content) = match sandbox.read_to_string_limited(rel, MAX_READ_BYTES) {
                Ok(res) => res,
                Err(e) => {
                    self.status_msg = format!("Sandbox import error ({}): {e}", path.display());
                    return;
                }
            };
            let ext = resolved
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            let parsed_mols: Vec<Molecule> = if ext == "sdf" || ext == "mol" {
                parse_sdf(&content)
                    .into_iter()
                    .filter_map(Result::ok)
                    .collect()
            } else if ext == "csv" {
                match parse_csv_batch(&content) {
                    Ok(items) => items.into_iter().filter_map(Result::ok).collect(),
                    Err(e) => {
                        self.status_msg = format!("CSV parse error ({}): {e}", resolved.display());
                        return;
                    }
                }
            } else {
                parse_smiles_batch(&content)
                    .into_iter()
                    .filter_map(|(_, _, r)| r.ok())
                    .collect()
            };

            for (idx, mol) in parsed_mols.into_iter().enumerate() {
                let Ok(smi) = mol.try_to_canonical_smiles() else {
                    continue;
                };
                let name = mol
                    .name
                    .unwrap_or_else(|| format!("Mol_{}_{}", loaded_files + 1, idx + 1));
                combined_lines.push(format!("{smi} {name}"));
            }
            loaded_files += 1;
        }

        let input = combined_lines.join("\n");
        self.run_command(Command::EvaluateBatchText { input });
        self.status_msg = format!(
            "Imported {} molecule(s) from {} dropped files",
            self.session.evaluations.len(),
            loaded_files
        );
    }

    fn handle_shortcuts(&mut self, ctx: &egui::Context) {
        let wants_kb = ctx.wants_keyboard_input();
        let mut open_palette = false;
        let mut run_batch = false;
        let mut load_examples = false;
        let mut do_undo = false;
        let mut do_redo = false;
        ctx.input(|i| {
            if i.modifiers.ctrl && i.key_pressed(egui::Key::K) {
                open_palette = true;
            }
            if i.modifiers.ctrl && i.key_pressed(egui::Key::Enter) {
                run_batch = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::E) {
                load_examples = true;
            }
            // Only intercept Ctrl+Z / Ctrl+Y for session undo when no TextEdit has focus
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::Z) {
                do_undo = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::Y) {
                do_redo = true;
            }
        });
        if open_palette {
            self.palette_open = !self.palette_open;
            self.palette_query.clear();
        }
        if run_batch {
            let input = self.session.smiles_input.clone();
            self.run_command(Command::EvaluateBatchText { input });
        }
        if load_examples {
            self.run_command(Command::LoadExampleBatch);
        }
        if do_undo {
            self.run_command(Command::Undo);
        }
        if do_redo {
            self.run_command(Command::Redo);
        }

        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .filter_map(|f| f.path.clone())
                .collect()
        });
        self.import_dropped_files(&dropped);
    }

    fn render_palette(&mut self, ctx: &egui::Context) {
        if !self.palette_open {
            return;
        }
        let mut close = false;
        let mut exec_id: Option<&'static str> = None;
        egui::Window::new("Command Palette (Ctrl+K)")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_TOP, [0.0, 64.0])
            .fixed_size([460.0, 300.0])
            .show(ctx, |ui| {
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut self.palette_query)
                        .hint_text("Type a command...")
                        .desired_width(f32::INFINITY),
                );
                resp.request_focus();
                ui.separator();

                let q = self.palette_query.to_ascii_lowercase();
                egui::ScrollArea::vertical().show(ui, |ui| {
                    for desc in CommandRegistry::catalog() {
                        if !q.is_empty()
                            && !desc.title.to_ascii_lowercase().contains(&q)
                            && !desc.category.to_ascii_lowercase().contains(&q)
                        {
                            continue;
                        }
                        ui.horizontal(|ui| {
                            if ui
                                .selectable_label(
                                    false,
                                    format!("[{}] {}", desc.category, desc.title),
                                )
                                .clicked()
                            {
                                exec_id = Some(desc.id);
                                close = true;
                            }
                            if !desc.shortcut.is_empty() {
                                ui.with_layout(
                                    egui::Layout::right_to_left(egui::Align::Center),
                                    |ui| {
                                        ui.weak(desc.shortcut);
                                    },
                                );
                            }
                        });
                    }
                });
                ui.separator();
                if ui.button("Close (Esc)").clicked()
                    || ui.input(|i| i.key_pressed(egui::Key::Escape))
                {
                    close = true;
                }
            });

        if let Some(id) = exec_id {
            match id {
                "load_examples" => self.run_command(Command::LoadExampleBatch),
                "run_batch" => {
                    let input = self.session.smiles_input.clone();
                    self.run_command(Command::EvaluateBatchText { input });
                }
                "filter_lipinski_0" => self.run_command(Command::FilterBatch {
                    criteria: FilterCriteria {
                        max_lipinski_violations: Some(0),
                        require_bbb: None,
                        require_gi_high: None,
                        no_pains: true,
                        no_brenk: false,
                    },
                }),
                "filter_bbb_yes" => self.run_command(Command::FilterBatch {
                    criteria: FilterCriteria {
                        max_lipinski_violations: None,
                        require_bbb: Some(true),
                        require_gi_high: None,
                        no_pains: false,
                        no_brenk: false,
                    },
                }),
                "clear_session" => self.run_command(Command::ClearSession),
                "undo" => self.run_command(Command::Undo),
                "redo" => self.run_command(Command::Redo),
                _ => {}
            }
        }
        if close {
            self.palette_open = false;
        }
    }
}

impl eframe::App for AdmetCraftApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if let Ok(mut guard) = self.async_status.lock()
            && let Some(msg) = guard.take()
        {
            self.status_msg = msg;
        }

        if self.dark_mode {
            ctx.set_visuals(egui::Visuals::dark());
        } else {
            ctx.set_visuals(egui::Visuals::light());
        }

        self.handle_shortcuts(ctx);
        self.render_palette(ctx);

        egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading(RichText::new("AdmetCraft").strong());
                ui.label(RichText::new("Local-First ADMET Workbench").weak());
                ui.separator();

                if ui
                    .button(
                        RichText::new("▶ Run!")
                            .strong()
                            .color(Color32::from_rgb(180, 40, 40)),
                    )
                    .clicked()
                {
                    let input = self.session.smiles_input.clone();
                    self.run_command(Command::EvaluateBatchText { input });
                }
                if ui.button("Examples").clicked() {
                    self.run_command(Command::LoadExampleBatch);
                }
                if ui.button("Clear").clicked() {
                    self.run_command(Command::ClearSession);
                }
                ui.separator();
                if ui
                    .add_enabled(self.session.can_undo(), egui::Button::new("Undo"))
                    .clicked()
                {
                    self.run_command(Command::Undo);
                }
                if ui
                    .add_enabled(self.session.can_redo(), egui::Button::new("Redo"))
                    .clicked()
                {
                    self.run_command(Command::Redo);
                }
                ui.separator();
                if ui.button("CSV").clicked() {
                    self.export(ExportFormat::Csv, "admetcraft_report.csv");
                }
                if ui.button("JSON").clicked() {
                    self.export(ExportFormat::Json, "admetcraft_report.json");
                }
                if ui.button("PDF").clicked() {
                    self.export(ExportFormat::Pdf, "admetcraft_report.pdf");
                }
                if ui.button("Egg SVG").clicked() {
                    self.export(ExportFormat::SvgBoiledEgg, "admetcraft_boiled_egg.svg");
                }
                ui.separator();
                if ui
                    .button(if self.dark_mode { "Light" } else { "Dark" })
                    .clicked()
                {
                    self.dark_mode = !self.dark_mode;
                }
                if ui.button("Palette (Ctrl+K)").clicked() {
                    self.palette_open = true;
                    self.palette_query.clear();
                }
            });
        });

        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(format!(
                    "Molecules: {}  |  Errors: {}  |  {}",
                    self.session.evaluations.len(),
                    self.session.failures.len(),
                    self.status_msg
                ));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.add(
                        egui::TextEdit::singleline(&mut self.export_dir)
                            .desired_width(180.0)
                            .hint_text("Sandbox dir"),
                    );
                    ui.label("Sandbox dir:");
                });
            });
        });

        egui::SidePanel::left("input_panel")
            .resizable(true)
            .default_width(360.0)
            .min_width(280.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.strong("2D Molecular Sketcher");
                    self.sketcher.ui(ui, self.dark_mode);
                    ui.horizontal(|ui| {
                        ui.label("Sketcher SMILES:");
                        ui.monospace(&self.sketcher.cached_smiles);
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Add to SMILES List >>").clicked()
                            && !self.sketcher.cached_smiles.is_empty()
                        {
                            let line = format!(
                                "{} Sketched_{}",
                                self.sketcher.cached_smiles,
                                self.session.evaluations.len() + 1
                            );
                            if !self.session.smiles_input.is_empty()
                                && !self.session.smiles_input.ends_with('\n')
                            {
                                self.session.smiles_input.push('\n');
                            }
                            self.session.smiles_input.push_str(&line);
                            let input = self.session.smiles_input.clone();
                            self.run_command(Command::EvaluateBatchText { input });
                        }
                        if ui.button("Clear Sketcher").clicked() {
                            self.sketcher.clear();
                        }
                    });

                    ui.separator();
                    ui.strong("SMILES List (SMILES [Optional Name])");
                    ui.add(
                        egui::TextEdit::multiline(&mut self.session.smiles_input)
                            .desired_rows(9)
                            .desired_width(f32::INFINITY)
                            .font(egui::TextStyle::Monospace),
                    );

                    ui.separator();
                    ui.strong("Import .smi / .csv / .sdf (sandboxed or drag & drop)");
                    ui.horizontal(|ui| {
                        ui.add(
                            egui::TextEdit::singleline(&mut self.import_path)
                                .desired_width(210.0)
                                .hint_text("molecules.smi"),
                        );
                        if ui.button("Load").clicked() && !self.import_path.trim().is_empty() {
                            let p = PathBuf::from(self.import_path.trim());
                            self.import_file_path(&p, false);
                        }
                    });

                    if !self.session.failures.is_empty() {
                        ui.separator();
                        ui.colored_label(
                            Color32::from_rgb(200, 60, 60),
                            format!("Parse warnings ({}):", self.session.failures.len()),
                        );
                        for f in &self.session.failures {
                            ui.small(format!("Line {}: {}", f.line, f.error));
                        }
                    }
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.selectable_value(&mut self.active_tab, MainTab::Cards, "Molecule Cards");
                ui.selectable_value(&mut self.active_tab, MainTab::BoiledEgg, "BOILED-Egg Plot");
                ui.selectable_value(&mut self.active_tab, MainTab::Table, "Batch Summary Table");
                ui.separator();
                if ui.small_button("Filter: Lipinski 0").clicked() {
                    self.run_command(Command::FilterBatch {
                        criteria: FilterCriteria {
                            max_lipinski_violations: Some(0),
                            require_bbb: None,
                            require_gi_high: None,
                            no_pains: false,
                            no_brenk: false,
                        },
                    });
                }
                if ui.small_button("Filter: BBB+").clicked() {
                    self.run_command(Command::FilterBatch {
                        criteria: FilterCriteria {
                            max_lipinski_violations: None,
                            require_bbb: Some(true),
                            require_gi_high: None,
                            no_pains: false,
                            no_brenk: false,
                        },
                    });
                }
                if ui.small_button("Filter: No PAINS").clicked() {
                    self.run_command(Command::FilterBatch {
                        criteria: FilterCriteria {
                            max_lipinski_violations: None,
                            require_bbb: None,
                            require_gi_high: None,
                            no_pains: true,
                            no_brenk: false,
                        },
                    });
                }
            });
            ui.separator();

            if self.session.evaluations.is_empty() {
                ui.centered_and_justified(|ui| {
                    ui.label("No molecules evaluated yet. Click 'Examples' or enter SMILES and press 'Run!'.");
                });
                return;
            }

            match self.active_tab {
                MainTab::Cards => self.render_cards_tab(ui),
                MainTab::BoiledEgg => {
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut self.show_egg_names, "Show Molecule Names");
                        ui.weak("(Hover a point for 2D structure; click to inspect its card)");
                    });
                    if let Some(clicked) = draw_boiled_egg_interactive(
                        ui,
                        &self.session.evaluations,
                        &self.mol_cache,
                        self.show_egg_names,
                        self.session.selected_index,
                        self.dark_mode,
                    ) {
                        self.session.selected_index = Some(clicked);
                        if let Some(ev) = self.session.evaluations.get(clicked) {
                            let _ = self.sketcher.load_smiles(&ev.canonical_smiles);
                        }
                        self.active_tab = MainTab::Cards;
                    }
                }
                MainTab::Table => self.render_table_tab(ui),
            }
        });
    }
}

impl AdmetCraftApp {
    fn render_cards_tab(&mut self, ui: &mut egui::Ui) {
        let mut next_selected: Option<usize> = None;
        egui::ScrollArea::horizontal()
            .id_salt("mol_chips")
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    for (idx, ev) in self.session.evaluations.iter().enumerate() {
                        let is_sel = self.session.selected_index == Some(idx);
                        if ui.selectable_label(is_sel, &ev.name).clicked() {
                            next_selected = Some(idx);
                        }
                    }
                });
            });
        if let Some(idx) = next_selected {
            self.session.selected_index = Some(idx);
            if let Some(ev) = self.session.evaluations.get(idx) {
                let _ = self.sketcher.load_smiles(&ev.canonical_smiles);
            }
        }
        ui.separator();

        let sel_idx = self.session.selected_index.unwrap_or(0);
        let Some(eval) = self.session.evaluations.get(sel_idx).cloned() else {
            return;
        };

        egui::ScrollArea::vertical()
            .id_salt("card_body")
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.heading(&eval.name);
                    ui.label(format!("({})", eval.physicochemical.formula));
                    if ui.small_button("Load into Sketcher").clicked() {
                        let _ = self.sketcher.load_smiles(&eval.canonical_smiles);
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Canonical SMILES:");
                    ui.monospace(&eval.canonical_smiles);
                });
                if !eval.warnings.is_empty() {
                    for w in &eval.warnings {
                        ui.colored_label(Color32::from_rgb(180, 95, 10), format!("Note: {w}"));
                    }
                }
                ui.add_space(6.0);

                ui.horizontal(|ui| {
                    ui.vertical(|ui| {
                        ui.strong("2D Structure");
                        let (rect, _) =
                            ui.allocate_exact_size(Vec2::new(260.0, 220.0), Sense::hover());
                        let bg = if self.dark_mode {
                            Color32::from_rgb(30, 32, 36)
                        } else {
                            Color32::WHITE
                        };
                        ui.painter_at(rect).rect_filled(rect, 4.0, bg);
                        if let Some(mol) = self.mol_cache.get(&eval.canonical_smiles) {
                            draw_cached_molecule_2d(
                                &ui.painter_at(rect),
                                mol,
                                rect,
                                self.dark_mode,
                            );
                        }
                    });
                    ui.add_space(16.0);
                    ui.vertical(|ui| {
                        ui.strong("Bioavailability Radar");
                        draw_bioavailability_radar(ui, &eval, 220.0, self.dark_mode);
                    });
                });

                ui.add_space(8.0);
                render_six_sections(ui, &eval);
            });
    }

    fn render_table_tab(&mut self, ui: &mut egui::Ui) {
        let mut picked: Option<usize> = None;
        egui::ScrollArea::both().show(ui, |ui| {
            egui::Grid::new("batch_summary_grid")
                .striped(true)
                .spacing([14.0, 6.0])
                .show(ui, |ui| {
                    for h in [
                        "Molecule",
                        "Formula",
                        "MW",
                        "TPSA",
                        "cLogP",
                        "ESOL LogS",
                        "GI",
                        "BBB",
                        "P-gp (exp.)",
                        "Lipinski",
                        "PAINS",
                        "SA (proxy)",
                    ] {
                        ui.strong(h);
                    }
                    ui.end_row();

                    for (idx, e) in self.session.evaluations.iter().enumerate() {
                        if ui
                            .selectable_label(self.session.selected_index == Some(idx), &e.name)
                            .clicked()
                        {
                            picked = Some(idx);
                        }
                        ui.label(&e.physicochemical.formula);
                        ui.label(format!("{:.2}", e.physicochemical.molecular_weight));
                        ui.label(format!("{:.2}", e.physicochemical.tpsa));
                        ui.label(format!("{:.2}", e.lipophilicity.consensus_log_p));
                        ui.label(format!("{:.2}", e.water_solubility.esol.log_s));
                        ui.label(e.pharmacokinetics.gi_absorption.as_str());
                        ui.label(yes_no(e.pharmacokinetics.bbb_permeant));
                        ui.label(yes_no(e.pharmacokinetics.pgp_substrate));
                        ui.label(format!(
                            "{} ({})",
                            yes_no(e.druglikeness.lipinski.pass),
                            e.druglikeness.lipinski.violations
                        ));
                        ui.label(format!(
                            "{}/{}",
                            e.medicinal_chemistry.pains_alerts.len(),
                            e.medicinal_chemistry.pains_catalog_size
                        ));
                        ui.label(format!(
                            "{:.2}",
                            e.medicinal_chemistry.synthetic_accessibility
                        ));
                        ui.end_row();
                    }
                });
        });
        if let Some(idx) = picked {
            self.session.selected_index = Some(idx);
            if let Some(ev) = self.session.evaluations.get(idx) {
                let _ = self.sketcher.load_smiles(&ev.canonical_smiles);
            }
            self.active_tab = MainTab::Cards;
        }
    }
}

fn render_six_sections(ui: &mut egui::Ui, e: &MoleculeEvaluation) {
    let p = &e.physicochemical;
    let l = &e.lipophilicity;
    let w = &e.water_solubility;
    let pk = &e.pharmacokinetics;
    let d = &e.druglikeness;
    let mc = &e.medicinal_chemistry;

    egui::CollapsingHeader::new("Physicochemical Properties")
        .default_open(true)
        .show(ui, |ui| {
            kv_grid(
                ui,
                "phys_grid",
                &[
                    ("Formula", p.formula.clone()),
                    (
                        "Molecular weight",
                        format!("{:.2} g/mol", p.molecular_weight),
                    ),
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
        });

    egui::CollapsingHeader::new("Lipophilicity")
        .default_open(true)
        .show(ui, |ui| {
            kv_grid(
                ui,
                "lipo_grid",
                &[
                    ("Log Po/w (WLOGP)", format!("{:.2}", l.wlogp)),
                    ("Log Po/w (XLOGP3 heur.)", format!("{:.2}", l.xlogp3)),
                    ("Log Po/w (MLOGP)", format!("{:.2}", l.mlogp)),
                    (
                        "Solvation est. (exp.)",
                        format!("{:.2}", l.solvation_logp_est),
                    ),
                    (
                        "Topological est. (exp.)",
                        format!("{:.2}", l.topological_logp_est),
                    ),
                    (
                        "Consensus Log Po/w (3-model)",
                        format!("{:.2}", l.consensus_log_p),
                    ),
                ],
            );
        });

    egui::CollapsingHeader::new("Water Solubility")
        .default_open(true)
        .show(ui, |ui| {
            kv_grid(
                ui,
                "sol_grid",
                &[
                    (
                        "Log S (ESOL)",
                        format!(
                            "{:.2} | {:.2e} mg/ml; {:.2e} mol/l | {}",
                            w.esol.log_s,
                            w.esol.solubility_mg_ml,
                            w.esol.solubility_mol_l,
                            w.esol.class.as_str()
                        ),
                    ),
                    (
                        "Log S (Ali)",
                        format!(
                            "{:.2} | {:.2e} mg/ml; {:.2e} mol/l | {}",
                            w.ali.log_s,
                            w.ali.solubility_mg_ml,
                            w.ali.solubility_mol_l,
                            w.ali.class.as_str()
                        ),
                    ),
                    (
                        "Log S (Topo exp.)",
                        format!(
                            "{:.2} | {:.2e} mg/ml; {:.2e} mol/l | {}",
                            w.topological_est.log_s,
                            w.topological_est.solubility_mg_ml,
                            w.topological_est.solubility_mol_l,
                            w.topological_est.class.as_str()
                        ),
                    ),
                ],
            );
        });

    egui::CollapsingHeader::new("Pharmacokinetics")
        .default_open(true)
        .show(ui, |ui| {
            kv_grid(
                ui,
                "pk_grid",
                &[
                    ("GI absorption", pk.gi_absorption.as_str().to_string()),
                    ("BBB permeant", yes_no(pk.bbb_permeant).to_string()),
                    (
                        "P-gp substrate (exp. heur.)",
                        yes_no(pk.pgp_substrate).to_string(),
                    ),
                    (
                        "CYP1A2 inhibitor (exp. heur.)",
                        yes_no(pk.cyp1a2_inhibitor).to_string(),
                    ),
                    (
                        "CYP2C19 inhibitor (exp. heur.)",
                        yes_no(pk.cyp2c19_inhibitor).to_string(),
                    ),
                    (
                        "CYP2C9 inhibitor (exp. heur.)",
                        yes_no(pk.cyp2c9_inhibitor).to_string(),
                    ),
                    (
                        "CYP2D6 inhibitor (exp. heur.)",
                        yes_no(pk.cyp2d6_inhibitor).to_string(),
                    ),
                    (
                        "CYP3A4 inhibitor (exp. heur.)",
                        yes_no(pk.cyp3a4_inhibitor).to_string(),
                    ),
                    (
                        "Log Kp (skin permeation)",
                        format!("{:.2} cm/s", pk.log_kp_cm_s),
                    ),
                ],
            );
        });

    egui::CollapsingHeader::new("Druglikeness")
        .default_open(true)
        .show(ui, |ui| {
            kv_grid(
                ui,
                "drug_grid",
                &[
                    ("Lipinski", fmt_rule(&d.lipinski)),
                    ("Ghose", fmt_rule(&d.ghose)),
                    ("Veber", fmt_rule(&d.veber)),
                    ("Egan", fmt_rule(&d.egan)),
                    ("Muegge", fmt_rule(&d.muegge)),
                    (
                        "Bioavailability Score",
                        format!("{:.2}", d.bioavailability_score),
                    ),
                ],
            );
        });

    egui::CollapsingHeader::new("Medicinal Chemistry")
        .default_open(true)
        .show(ui, |ui| {
            let pains_txt = if mc.pains_alerts.is_empty() {
                format!(
                    "0 alerts ({} curated patterns checked)",
                    mc.pains_catalog_size
                )
            } else {
                format!(
                    "{}/{} alert(s): {}",
                    mc.pains_alerts.len(),
                    mc.pains_catalog_size,
                    mc.pains_alerts.join(", ")
                )
            };
            let brenk_txt = if mc.brenk_alerts.is_empty() {
                format!(
                    "0 alerts ({} curated patterns checked)",
                    mc.brenk_catalog_size
                )
            } else {
                format!(
                    "{}/{} alert(s): {}",
                    mc.brenk_alerts.len(),
                    mc.brenk_catalog_size,
                    mc.brenk_alerts.join(", ")
                )
            };
            kv_grid(
                ui,
                "mc_grid",
                &[
                    ("PAINS", pains_txt),
                    ("Brenk", brenk_txt),
                    ("Leadlikeness", fmt_rule(&mc.leadlikeness)),
                    (
                        "SA score (heuristic proxy)",
                        format!("{:.2}", mc.synthetic_accessibility),
                    ),
                ],
            );
        });
}

fn kv_grid(ui: &mut egui::Ui, id: &str, rows: &[(&str, String)]) {
    egui::Grid::new(id)
        .striped(true)
        .spacing([20.0, 4.0])
        .show(ui, |ui| {
            for (k, v) in rows {
                ui.label(*k);
                ui.strong(v);
                ui.end_row();
            }
        });
}

fn fmt_rule(r: &admet_calc::RuleResult) -> String {
    if r.details.is_empty() {
        format!("{}; {} violations", yes_no(r.pass), r.violations)
    } else {
        format!(
            "{}; {} violation(s): {}",
            yes_no(r.pass),
            r.violations,
            r.details.join(", ")
        )
    }
}

fn yes_no(b: bool) -> &'static str {
    if b { "Yes" } else { "No" }
}

fn main() -> Result<(), eframe::Error> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("AdmetCraft — Local-First ADMET Workbench")
            .with_inner_size([1280.0, 820.0]),
        ..Default::default()
    };
    eframe::run_native(
        "AdmetCraft",
        options,
        Box::new(|_cc| Ok(Box::new(AdmetCraftApp::default()))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sketcher_loads_and_syncs_smiles() -> Result<(), chem_core::ChemError> {
        let mut sk = MoleculeSketcher::default();
        sk.load_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        assert!(!sk.cached_smiles.is_empty());
        assert_eq!(sk.mol.heavy_atom_count(), 13);
        Ok(())
    }
}
