# AdmetCraft

**Local-First & Agent-Native Cheminformatics, ADMET, Druglikeness & Safety Workbench in 100% Pure Rust.**

**AdmetCraft** computes molecular descriptors, ionization / $\log D_{7.4}$, lipophilicity, aqueous solubility, pharmacokinetics, quantitative druglikeness (**QED**, **Pfizer CNS MPO**, **Pfizer 3/75**, **GSK 4/400**, **Golden Triangle**, **Ro3**, **bRo5**), medicinal chemistry filters (**PAINS** / **Brenk**), and early safety alerts (**hERG** / **AMES** / **DILI**) completely offline across all CPU cores—with zero data leaving your machine and zero batch-size server limits.

## Architecture (Layered Cargo Workspace)

- `crates/chem-core` (L0): SMILES & V2000 MOL/SDF parsers, RFC-4180 CSV importer, molecular graph, GF(2) SSSR ring perception, Hückel $4n+2$ aromaticity perception, iterative stack-safe canonical SMILES (reusable ring closure numbers `1..=99`), automatic desalting (with warnings), stereochemistry warnings (`@/@@`, `/\`), 2D coordinate layout, and Daylight recursive SMARTS matcher (with full 2-letter bracket element support `[Br]`, `[Cr]`, etc.).
- `crates/admet-calc` (L1): Parallel (`rayon`) descriptor, ionization ($\text{p}K_a$/$\log D_{7.4}$), ADMET, QED/CNS MPO, industrial medchem rules, PAINS/Brenk/hERG/AMES/DILI, matched-pair comparison (`compare_smiles`), and shared Bioavailability Radar / BOILED-Egg geometry engine.
- `crates/admet-report` (L2): CSV (roundtrip-compatible, 60+ columns), JSON, SVG (Radar & BOILED-Egg with 36-char label truncation and XML noncharacter sanitization), and pure-Rust multi-page PDF 1.4 vector report writer with **zlib `/Filter /FlateDecode`** stream compression (~4.9× smaller files), WinAnsiEncoding, and paginated summary tables.
- `crates/admet-engine` (L3): Session state, undo/redo history, unified `CommandRegistry`, multi-criteria filtering (`Lipinski`, `QED`, `CNS MPO`, `BBB`, `GI`, `PAINS`, `Brenk`, `hERG`, `AMES`, `DILI`), and `PathSandbox` (`openat2(RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS)`, `O_NOFOLLOW`, FIFO rejection, explicit `--create-write-root` guard).
- `apps/admetcraft-cli` (L4): Headless CLI (`eval`, `profile`, `compare`, `batch`) and strict **JSON-RPC 2.0 / MCP (2025-06-18)** server (`stdio` and authenticated loopback-default `--tcp`) with `outputSchema` on all 5 tools.
- `apps/admetcraft` (L5): Native `egui`/`eframe` desktop GUI with 2D molecular sketcher (synced across session/filter/history operations), KPI badge strip (`QED`, `CNS MPO`, `LogD7.4`, `hERG`, alerts), native `zenity`/`kdialog` file picker, background async exports with `ctx.request_repaint()`, interactive *BOILED-Egg* scatter plot, expanded Batch Summary Table, Light/Dark themes, and `Ctrl+K` Command Palette.

## Scientific Scope, Benchmark Parity & Model Transparency

Differential verification against **RDKit 2026.9.2**, **EPAM Indigo 1.49.0**, and the **15-drug SwissADME reference corpus** across 45 diverse molecules (clinical drugs, fragments, bRo5 macrocycles, acids/bases/zwitterions, and toxicophores):

| Category | Implemented Models | Verified Parity / Notes & Limitations |
|---|---|---|
| Physicochemical | Formula, IUPAC MW, monoisotopic Exact Mass, heavy/aromatic atoms, rings/aromatic rings, DBE, $Fsp^3$, rotatable bonds, HBA/HBD, Wildman-Crippen MR, Ertl (2000) TPSA (both **S & P** SwissADME and **N & O only** RDKit) | **Exact Mass** max err `0.00005 Da`; **MW** MAE `0.0038 Da` vs RDKit (`max 0.005 Da` vs Indigo); **TPSA (S&P & N/O)** **`MAE = 0.0000` (45/45 exact)**; **MR** $r = 0.9997$; **HBD / Rings / Arom. Rings** 45/45 exact. |
| Lipophilicity & Ionization | `WLOGP` (Wildman & Crippen 1999), `MLOGP` (Moriguchi 1992), `XLOGP3 (heur.)`, `Consensus Log P`, plus physiological pH 7.4 `IonizationClass`, heuristic acidic/basic $\text{p}K_a$, charge@7.4, and Henderson-Hasselbalch `LogD 7.4` | **WLOGP** MAE `0.238`, $r = 0.9844$ vs RDKit. `solvation_logp_est`, `topological_logp_est`, and $\text{p}K_a$ are clearly labeled heuristics (not quantum/QM solvers). |
| Solubility | `ESOL` (Delaney 2004), `Ali` (Ali et al. 2012), `topological_est` (experimental heuristic) | Quantitative $\log S$, mg/mL, mol/L, and solubility class. |
| Pharmacokinetics | `BOILED-Egg` (Daina & Zoete 2016) GI & BBB ellipses, Potts & Guy (1992) $\log K_p$ (cm/s) | **P-gp substrate** and **CYP1A2/2C19/2C9/2D6/3A4** flags are **experimental rule-based heuristics**, not trained SVMs. |
| Druglikeness & MedChem Rules | **QED** (exact Bickerton et al. 2012 ADS curves), **Pfizer CNS MPO** (Wager et al. 2010, `0..6`), `Lipinski`, `Ghose`, `Veber`, `Egan`, `Muegge`, **Pfizer 3/75**, **GSK 4/400**, **Golden Triangle**, **Rule of Three (Ro3)**, **Beyond Ro5 (`bRo5`)**, Abbott `Bioavailability Score` (Martin 2005) | **QED** MAE `0.0267`, $r = 0.9723$ vs RDKit `Chem.QED.qed`. |
| Medicinal Chemistry & Safety | `PAINS` (40 curated A/B/C SMARTS), `Brenk` (56 curated filters incl. `>_2_ester_groups` count rule), `Leadlikeness` (Teague et al.), `SA score (proxy)`, plus heuristic **hERG**, **AMES** (15 Benigni-Bossa rules), and **DILI** (14 bioactivation rules) | **PAINS** `45/45` (100%) and **Brenk** `45/45` (100%) exact binary agreement with RDKit `FilterCatalog` on the 45-molecule benchmark. Catalog sizes are reported explicitly; hERG/AMES/DILI and SA are structural/descriptor heuristics. |

## Security Guarantees

- `#![forbid(unsafe_code)]`, `#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)]` workspace-wide.
- **Directory Sandbox Isolation:** `PathSandbox` pins `--read-root` / `--write-root` directory file descriptors and opens files atomically via kernel `openat2(..., RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS)` (with `O_NOFOLLOW | O_NONBLOCK` and regular-file verification before truncation), preventing intermediate-directory symlink races and FIFOs. Non-existent `--write-root` paths are refused unless `--create-write-root` is passed.
- **Strict JSON-RPC 2.0 / MCP & DoS Guards:** Enforces `"jsonrpc": "2.0"` (`-32600`), strict parameter validation (`-32602`), silent notifications without `id`, preserved `"id": null`, `1 MiB` message cap, `1,000`-molecule MCP batch cap, `initialize` lifecycle, loopback-only `--tcp` by default, and mandatory `--mcp-token` for TCP clients.

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

# Concise human-readable ADMET & Safety card
admetcraft-cli profile "CC(C)Cc1ccc(C(C)C(=O)O)cc1" --name Ibuprofen

# Matched Molecular Pair / Lead optimization comparison (JSON deltas)
admetcraft-cli compare "CC(C)Cc1ccc(C(C)C(=O)O)cc1" "CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1" --name-a Ibuprofen --name-b Ibu_Tetrazole

# Parallel batch evaluation (.smi, .csv, .sdf -> .csv, .json, .pdf, .svg)
admetcraft-cli batch molecules.smi --read-root ./inputs --write-root ./outputs --out report.pdf --format pdf

# Sandboxed MCP server over stdio (5 tools with outputSchema)
admetcraft-cli mcp --read-root ./inputs --write-root ./outputs

# Authenticated loopback TCP MCP server
admetcraft-cli mcp --read-root ./inputs --write-root ./outputs --tcp 127.0.0.1:7878 --mcp-token "<secret>"
```

## License

Dual-licensed under [MIT](LICENSE-MIT) OR [Apache-2.0](LICENSE-APACHE).
