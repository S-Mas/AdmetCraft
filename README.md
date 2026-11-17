# AdmetCraft

**Local-First & Agent-Native Cheminformatics, ADMET, and Druglikeness Workbench in 100% Pure Rust.**

**AdmetCraft** computes molecular descriptors, lipophilicity, aqueous solubility, pharmacokinetics, druglikeness rules, and medicinal chemistry substructure alerts completely offline across all CPU cores—with zero data leaving your machine and zero batch-size server limits.

## Architecture (Layered Cargo Workspace)

- `crates/chem-core` (L0): SMILES & V2000 MOL/SDF parsers, RFC-4180 CSV importer, molecular graph, GF(2) SSSR ring perception, Hückel $4n+2$ aromaticity perception, iterative stack-safe canonical SMILES (with reusable ring closure numbers `1..=99`), automatic desalting (with warnings), stereochemistry warnings (`@/@@`, `/\`), 2D coordinate layout, and recursive SMARTS matcher.
- `crates/admet-calc` (L1): Parallel (`rayon`) descriptor and ADMET calculation engine + shared Bioavailability Radar and BOILED-Egg geometry helpers.
- `crates/admet-report` (L2): CSV (roundtrip-compatible), JSON, SVG (Radar & BOILED-Egg), and pure-Rust multi-page PDF 1.4 vector report writer (WinAnsiEncoding for European accents and paginated summary tables).
- `crates/admet-engine` (L3): Session state, undo/redo history, unified `CommandRegistry`, and `PathSandbox` (`O_NOFOLLOW`, component-wise symlink rejection, FIFO/non-regular-file rejection, atomic `open_read` / `open_write` handles).
- `apps/admetcraft-cli` (L4): Headless CLI (`eval`, `batch`) and sandboxed **Model Context Protocol (MCP)** server (`stdio` and authenticated loopback-default `--tcp`).
- `apps/admetcraft` (L5): Native `egui`/`eframe` desktop GUI with 2D molecular sketcher, cached 2D depictions, background async exports, interactive *BOILED-Egg* scatter plot, Light/Dark themes, and `Ctrl+K` Command Palette.

## Scientific Scope & Model Transparency

AdmetCraft distinguishes clearly between published equations/parameters and experimental heuristic estimators:

| Category | Implemented Models | Notes & Limitations |
|---|---|---|
| Physicochemical | Formula, IUPAC MW, heavy/aromatic atoms, $Fsp^3$, rotatable bonds, HBA/HBD, Wildman-Crippen MR, Ertl (2000) TPSA (incl. S & P) | Desalts multi-component inputs (retaining largest fragment) and emits warnings when `@/@@` or `/\` stereo markers are present. |
| Lipophilicity | `WLOGP` (Wildman & Crippen 1999), `MLOGP` (Moriguchi et al. 1992), `XLOGP3 (heur.)` (atom-additive heuristic inspired by Cheng et al. 2007), `Consensus Log P` (mean of the 3 primary models) | `solvation_logp_est` and `topological_logp_est` are labeled experimental linear heuristics—**not** SwissADME iLOGP or Silicos-IT. |
| Solubility | `ESOL` (Delaney 2004), `Ali` (Ali et al. 2012), `topological_est` (experimental heuristic) | Quantitative $\log S$, mg/mL, mol/L, and solubility class. |
| Pharmacokinetics | `BOILED-Egg` (Daina & Zoete 2016) GI & BBB ellipses, Potts & Guy (1992) $\log K_p$ (cm/s) | **P-gp substrate** and **CYP1A2/2C19/2C9/2D6/3A4** flags are **experimental rule-based heuristics**, not trained SVMs. |
| Druglikeness | `Lipinski`, `Ghose`, `Veber`, `Egan`, `Muegge`, Abbott `Bioavailability Score` (Martin 2005) | Reports pass/fail, violation count, and exact violated criteria. |
| Medicinal Chemistry | `PAINS` (26 curated SMARTS from Baell & Holloway), `Brenk` (25 curated SMARTS), `Leadlikeness` (Teague et al.), `SA score (proxy)` | Catalog sizes (`26` PAINS / `25` Brenk, subsets of ~480 / ~105 full catalogs) are reported explicitly in every output. SA score is a heuristic complexity proxy without PubChem's 13M fragment table. |

## Security Guarantees

- `#![forbid(unsafe_code)]`, `#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)]` workspace-wide.
- **Directory Sandbox Isolation:** `PathSandbox` pins `--read-root` / `--write-root` directory file descriptors and opens files atomically via kernel `openat2(..., RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS)` (with `O_NOFOLLOW | O_NONBLOCK` and regular-file verification before truncation), preventing intermediate-directory symlink races and FIFOs.
- **DoS & Network Guards:** SMILES parser enforces `DEFAULT_MAX_ATOMS = 2000` with iterative stack-safe canonicalization; MCP caps requests at `1 MiB`, batches at `1,000` molecules, requires `initialize` lifecycle, binds `--tcp` to loopback only unless `--allow-remote` is passed, and requires `--mcp-token` for TCP clients.

## Building & Testing

```bash
cargo fmt --check
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

# Authenticated loopback TCP MCP server
admetcraft-cli mcp --read-root ./inputs --write-root ./outputs --tcp 127.0.0.1:7878 --mcp-token "<secret>"
```

## License

Dual-licensed under [MIT](LICENSE-MIT) OR [Apache-2.0](LICENSE-APACHE).
