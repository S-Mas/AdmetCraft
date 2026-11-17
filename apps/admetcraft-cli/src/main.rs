#![forbid(unsafe_code)]

mod mcp;

use std::path::{Path, PathBuf};

use admet_calc::{evaluate_batch, evaluate_smiles};
use admet_engine::PathSandbox;
use admet_report::{ExportFormat, write_report_file};
use chem_core::{parse_csv_batch, parse_sdf, parse_smiles_batch};
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Debug, Parser)]
#[command(
    name = "admetcraft-cli",
    version,
    about = "Local-first pure-Rust ADMET calculator and MCP server"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Evaluate a single SMILES string and print complete ADMET JSON.
    Eval {
        /// SMILES string (or "SMILES Name").
        smiles: String,
        /// Optional molecule name.
        #[arg(long)]
        name: Option<String>,
    },
    /// Process a batch file (.smi, .csv, .sdf) in parallel across all CPU cores.
    Batch {
        /// Input file (.smi, .csv, .sdf).
        input: PathBuf,
        /// Output file path.
        #[arg(long)]
        out: PathBuf,
        /// Output format.
        #[arg(long, value_enum, default_value = "csv")]
        format: CliFormat,
        /// Optional read sandbox root.
        #[arg(long)]
        read_root: Option<PathBuf>,
        /// Optional write sandbox root.
        #[arg(long)]
        write_root: Option<PathBuf>,
    },
    /// Start the Model Context Protocol (MCP) server on stdio or TCP.
    Mcp {
        /// Authorized directory for reading input files.
        #[arg(long)]
        read_root: PathBuf,
        /// Authorized directory for writing output reports.
        #[arg(long)]
        write_root: PathBuf,
        /// Optional TCP bridge address (e.g. 127.0.0.1:7878); defaults to stdio.
        #[arg(long)]
        tcp: Option<String>,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum CliFormat {
    Csv,
    Json,
    Pdf,
    Svg,
}

impl From<CliFormat> for ExportFormat {
    fn from(value: CliFormat) -> Self {
        match value {
            CliFormat::Csv => ExportFormat::Csv,
            CliFormat::Json => ExportFormat::Json,
            CliFormat::Pdf => ExportFormat::Pdf,
            CliFormat::Svg => ExportFormat::SvgBoiledEgg,
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Eval { smiles, name } => {
            let mut parts = smiles.trim().splitn(2, char::is_whitespace);
            let smi_token = parts.next().unwrap_or("");
            let inline_name = parts.next().map(str::trim).filter(|s| !s.is_empty());
            let final_name = name.as_deref().or(inline_name);
            let eval = evaluate_smiles(smi_token, final_name)?;
            let out = serde_json::to_string_pretty(&eval)?;
            println!("{out}");
        }
        Commands::Batch {
            input,
            out,
            format,
            read_root,
            write_root,
        } => {
            let resolved_in = if let Some(rr) = read_root.as_deref() {
                let sb = PathSandbox::new(Some(rr), None)?;
                sb.resolve_read(&input)?
            } else {
                input
            };
            let resolved_out = if let Some(wr) = write_root.as_deref() {
                let sb = PathSandbox::new(None, Some(wr))?;
                sb.resolve_write(&out)?
            } else {
                out
            };
            let evals = load_molecules_from_path(&resolved_in)?;
            write_report_file(&evals, &resolved_out, format.into())?;
            eprintln!(
                "Exported {} molecule(s) to {}",
                evals.len(),
                resolved_out.display()
            );
        }
        Commands::Mcp {
            read_root,
            write_root,
            tcp,
        } => {
            let sandbox = PathSandbox::new(Some(&read_root), Some(&write_root))?;
            if let Some(addr) = tcp {
                mcp::run_mcp_tcp(&addr, &sandbox)?;
            } else {
                mcp::run_mcp_stdio(&sandbox)?;
            }
        }
    }
    Ok(())
}

fn load_molecules_from_path(
    path: &Path,
) -> Result<Vec<admet_calc::MoleculeEvaluation>, Box<dyn std::error::Error>> {
    let content = std::fs::read_to_string(path)?;
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let mut mols = Vec::new();
    if ext == "sdf" || ext == "mol" {
        for (i, res) in parse_sdf(&content).into_iter().enumerate() {
            match res {
                Ok(m) => mols.push(m),
                Err(e) => eprintln!("Warning: SDF record {} skipped: {e}", i + 1),
            }
        }
    } else if ext == "csv" {
        for (i, res) in parse_csv_batch(&content)?.into_iter().enumerate() {
            match res {
                Ok(m) => mols.push(m),
                Err(e) => eprintln!("Warning: CSV row {} skipped: {e}", i + 2),
            }
        }
    } else {
        for (line, _, res) in parse_smiles_batch(&content) {
            match res {
                Ok(m) => mols.push(m),
                Err(e) => eprintln!("Warning: line {line} skipped: {e}"),
            }
        }
    }
    Ok(evaluate_batch(&mols))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn mcp_server_handles_evaluate_and_rejects_outside_sandbox() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("admetcraft_mcp_cli_test");
        let rdir = base.join("read");
        let wdir = base.join("write");
        let forbidden = base.join("forbidden");
        std::fs::create_dir_all(&rdir)?;
        std::fs::create_dir_all(&wdir)?;
        std::fs::create_dir_all(&forbidden)?;
        std::fs::write(forbidden.join("secret.smi"), "CCO Secret\n")?;

        let sb = PathSandbox::new(Some(&rdir), Some(&wdir))?;
        let requests = [
            serde_json::json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"evaluate_smiles","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O Aspirin"}}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"evaluate_smiles","arguments":{"file":"../forbidden/secret.smi"}}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"export_report","arguments":{"smiles":["CCO Ethanol"],"output_path":"../forbidden/evil.csv","format":"csv"}}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"export_report","arguments":{"smiles":["CCO Ethanol"],"output_path":"report.csv","format":"csv"}}}).to_string(),
        ].join("\n") + "\n";

        let reader = BufReader::new(requests.as_bytes());
        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(reader, &mut out_buf, &sb)?;
        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<serde_json::Value> = out_str
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 5);
        assert_eq!(lines[1]["result"]["isError"], false);
        assert_eq!(lines[1]["result"]["structuredContent"]["count"], 1);
        assert_eq!(lines[2]["result"]["isError"], true);
        assert_eq!(lines[3]["result"]["isError"], true);
        assert_eq!(lines[4]["result"]["isError"], false);
        assert!(wdir.join("report.csv").exists());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
