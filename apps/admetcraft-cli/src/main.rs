#![forbid(unsafe_code)]

mod mcp;

use std::path::PathBuf;

use admet_calc::{evaluate_batch, evaluate_smiles};
use admet_engine::{MAX_READ_BYTES, PathSandbox};
use admet_report::{ExportFormat, write_report_file, write_report_to_writer};
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
        /// Optional TCP bridge address (loopback by default; requires --mcp-token).
        #[arg(long)]
        tcp: Option<String>,
        /// Allow binding non-loopback IP addresses for --tcp (must also provide --mcp-token).
        #[arg(long)]
        allow_remote: bool,
        /// Bearer token required by --tcp clients during `initialize` (or env ADMETCRAFT_MCP_TOKEN).
        #[arg(long)]
        mcp_token: Option<String>,
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
            let (ext, content) = if let Some(rr) = read_root.as_deref() {
                let sb = PathSandbox::new(Some(rr), None)?;
                let (resolved_in, text) = sb.read_to_string_limited(&input, MAX_READ_BYTES)?;
                let ext = resolved_in
                    .extension()
                    .and_then(|e| e.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                (ext, text)
            } else {
                let ext = input
                    .extension()
                    .and_then(|e| e.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                (ext, std::fs::read_to_string(&input)?)
            };
            let evals = load_molecules_from_text(&ext, &content)?;
            if let Some(wr) = write_root.as_deref() {
                let sb = PathSandbox::new(None, Some(wr))?;
                let (resolved_out, mut file) = sb.open_write(&out)?;
                write_report_to_writer(&evals, &mut file, format.into())?;
                eprintln!(
                    "Exported {} molecule(s) to {}",
                    evals.len(),
                    resolved_out.display()
                );
            } else {
                write_report_file(&evals, &out, format.into())?;
                eprintln!("Exported {} molecule(s) to {}", evals.len(), out.display());
            }
        }
        Commands::Mcp {
            read_root,
            write_root,
            tcp,
            allow_remote,
            mcp_token,
        } => {
            let sandbox = PathSandbox::new(Some(&read_root), Some(&write_root))?;
            if let Some(addr) = tcp {
                let token = mcp_token.or_else(|| std::env::var("ADMETCRAFT_MCP_TOKEN").ok());
                mcp::run_mcp_tcp(&addr, &sandbox, allow_remote, token)?;
            } else {
                mcp::run_mcp_stdio(&sandbox)?;
            }
        }
    }
    Ok(())
}

fn load_molecules_from_text(
    ext: &str,
    content: &str,
) -> Result<Vec<admet_calc::MoleculeEvaluation>, Box<dyn std::error::Error>> {
    let mut mols = Vec::new();
    if ext == "sdf" || ext == "mol" {
        for (i, res) in parse_sdf(content).into_iter().enumerate() {
            match res {
                Ok(m) => mols.push(m),
                Err(e) => eprintln!("Warning: SDF record {} skipped: {e}", i + 1),
            }
        }
    } else if ext == "csv" {
        for (i, res) in parse_csv_batch(content)?.into_iter().enumerate() {
            match res {
                Ok(m) => mols.push(m),
                Err(e) => eprintln!("Warning: CSV row {} skipped: {e}", i + 2),
            }
        }
    } else {
        for (line, _, res) in parse_smiles_batch(content) {
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
    fn mcp_server_handles_evaluate_and_rejects_outside_sandbox()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("admetcraft_mcp_cli_test_p0");
        let _ = std::fs::remove_dir_all(&base);
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
        mcp::serve_stream_for_test(reader, &mut out_buf, &sb, None)?;
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

    #[test]
    fn mcp_enforces_auth_initialize_size_and_batch_limits() -> Result<(), Box<dyn std::error::Error>>
    {
        let base = std::env::temp_dir().join("admetcraft_mcp_limits_test");
        let _ = std::fs::remove_dir_all(&base);
        let rdir = base.join("read");
        let wdir = base.join("write");
        std::fs::create_dir_all(&rdir)?;
        std::fs::create_dir_all(&wdir)?;
        let sb = PathSandbox::new(Some(&rdir), Some(&wdir))?;

        // 1. Uninitialized call + wrong token must both be rejected
        let uninit_reqs = [
            serde_json::json!({"jsonrpc":"2.0","id":1,"method":"tools/list","params":{}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":2,"method":"initialize","params":{"token":"wrong"}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":3,"method":"initialize","params":{"token":"s3cr3t"}}).to_string(),
            serde_json::json!({"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"evaluate_smiles","arguments":{"smiles": vec!["CCO"; 1005]}}}).to_string(),
        ].join("\n") + "\n";

        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(
            BufReader::new(uninit_reqs.as_bytes()),
            &mut out_buf,
            &sb,
            Some("s3cr3t".to_string()),
        )?;
        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<serde_json::Value> = out_str
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0]["error"]["code"], -32002);
        assert_eq!(lines[1]["error"]["code"], -32001);
        assert!(lines[2]["result"]["serverInfo"]["name"] == "admetcraft-mcp");
        assert_eq!(lines[3]["result"]["isError"], true);

        // 2. Oversized >1 MiB line must be rejected with -32600
        let huge_line = "x".repeat(mcp::MAX_MCP_MESSAGE_BYTES + 128) + "\n";
        let mut huge_out = Vec::new();
        mcp::serve_stream_for_test(
            BufReader::new(huge_line.as_bytes()),
            &mut huge_out,
            &sb,
            None,
        )?;
        let resp: serde_json::Value = serde_json::from_str(String::from_utf8(huge_out)?.trim())?;
        assert_eq!(resp["error"]["code"], -32600);

        // 3. Non-loopback TCP without --allow-remote must be refused
        assert!(mcp::run_mcp_tcp("0.0.0.0:9999", &sb, false, Some("tok".to_string())).is_err());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
