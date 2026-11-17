# AdmetCraft

**Local-First & Agent-Native ADMET, Physicochemical, Druglikeness, and Medicinal Chemistry Suite in 100% Pure Rust.**

Inspired by the SwissADME workflow (`http://www.swissadme.ch`), **AdmetCraft** computes molecular descriptors, lipophilicity, water solubility, pharmacokinetics (*BOILED-Egg*, P-gp, CYP450, $\log K_p$), druglikeness filters, and medicinal chemistry alerts completely offline on your own CPU cores—with zero data leaving your machine and zero rate limits.

## Key Features

- **100% Pure Rust & Safe:** `#![forbid(unsafe_code)]`, `#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)]` across the entire workspace.
- **Layered Architecture:**
  - `crates/chem-core` (L0): SMILES & V2000 MOL/SDF parsers, molecular graph, GF(2) SSSR ring perception, Hückel $4n+2$ aromaticity perception, canonical SMILES, 2D coordinate layout, and recursive SMARTS matcher.
  - `crates/admet-calc` (L1): Parallel (`rayon`) calculation of physicochemical descriptors (including Ertl 2000 TPSA with S & P and Wildman-Crippen MR), 5 $\log P_{o/w}$ models (`iLOGP`, `XLOGP3`, `WLOGP`, `MLOGP`, `SILICOS-IT`) + Consensus $\log P$, 3 $\log S$ models (`ESOL`, `Ali`, `SILICOS-IT`), Pharmacokinetics (*BOILED-Egg* GI/BBB, P-gp, CYP1A2/2C19/2C9/2D6/3A4, Potts & Guy $\log K_p$), Druglikeness (`Lipinski`, `Ghose`, `Veber`, `Egan`, `Muegge`, Abbott Bioavailability Score), Medicinal Chemistry (`PAINS` A/B/C, `Brenk`, `Leadlikeness`, `Synthetic Accessibility`), and 6-axis Bioavailability Radar.
  - `crates/admet-report` (L2): CSV, JSON, SVG (Bioavailability Radar & BOILED-Egg), and pure-Rust PDF 1.4 vector report export.
  - `crates/admet-engine` (L3): Session state, undo/redo history, unified `CommandRegistry`, and strict `PathSandbox` canonicalization.
  - `apps/admetcraft-cli` (L4): Headless CLI (`eval`, `batch`) and sandboxed **Model Context Protocol (MCP)** server (`stdio` and `--tcp` bridge).
  - `apps/admetcraft` (L5): Native `egui`/`eframe` desktop app with interactive 2D sketcher, batch/file import (`.smi`, `.csv`, `.sdf` + drag-and-drop), 6-axis Bioavailability Radar, interactive *BOILED-Egg* chart, Light/Dark themes, and `Ctrl+K` Command Palette.

## Building & Testing

```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo build --workspace --release
```

## CLI & MCP Usage

```bash
# Single SMILES evaluation (JSON output)
admetcraft-cli eval "CC(=O)Oc1ccccc1C(=O)O Aspirin"

# Parallel batch evaluation (.smi, .csv, .sdf -> .csv, .json, .pdf, .svg)
admetcraft-cli batch molecules.smi --out report.pdf --format pdf

# Sandboxed MCP server over stdio
admetcraft-cli mcp --read-root ./inputs --write-root ./outputs
```

### Exposed MCP Tools
- `evaluate_smiles` — Evaluate SMILES strings or authorized files inside `--read-root`.
- `filter_druglikeness` — Filter compounds by Lipinski violations, BBB permeation, GI absorption, and PAINS/Brenk alerts.
- `inspect_molecule` — Return structural graph breakdown, rings, 2D coordinates, and PAINS/Brenk alerts.
- `export_report` — Export CSV, JSON, PDF, or BOILED-Egg SVG reports strictly inside `--write-root`.

## Packaging & Installation

```bash
./packaging/install.sh
```

This builds release binaries and installs `admetcraft` and `admetcraft-cli` into `~/.local/bin/` along with a FreeDesktop `.desktop` entry and SVG icon. An Arch Linux `packaging/PKGBUILD` is also provided.

## License

Dual-licensed under MIT OR Apache-2.0.
