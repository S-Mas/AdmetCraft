use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::path::Path;

use admet_calc::{MoleculeEvaluation, evaluate_batch, evaluate_smiles};
use admet_engine::{FilterCriteria, PathSandbox};
use admet_report::{ExportFormat, write_report_file};
use chem_core::{parse_csv_batch, parse_sdf, parse_smiles, parse_smiles_batch};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

#[derive(Debug, Deserialize)]
struct JsonRpcRequest {
    #[allow(dead_code)]
    jsonrpc: Option<String>,
    id: Option<Value>,
    method: String,
    #[serde(default)]
    params: Value,
}

#[derive(Debug, Serialize)]
struct JsonRpcResponse {
    jsonrpc: &'static str,
    id: Value,
    #[serde(skip_serializing_if = "Option::is_none")]
    result: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<JsonRpcError>,
}

#[derive(Debug, Serialize)]
struct JsonRpcError {
    code: i64,
    message: String,
}

/// Run the MCP server over `stdin`/`stdout`.
pub fn run_mcp_stdio(sandbox: &PathSandbox) -> Result<(), Box<dyn std::error::Error>> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let reader = BufReader::new(stdin.lock());
    let mut writer = stdout.lock();
    serve_stream(reader, &mut writer, sandbox)
}

/// Run the MCP server over a TCP bridge address (e.g. `127.0.0.1:7878`).
pub fn run_mcp_tcp(addr: &str, sandbox: &PathSandbox) -> Result<(), Box<dyn std::error::Error>> {
    let listener = TcpListener::bind(addr)?;
    for stream_res in listener.incoming() {
        let Ok(stream) = stream_res else {
            continue;
        };
        let Ok(reader_stream) = stream.try_clone() else {
            continue;
        };
        let reader = BufReader::new(reader_stream);
        let mut writer = stream;
        let _ = serve_stream(reader, &mut writer, sandbox);
    }
    Ok(())
}

#[cfg(test)]
pub fn serve_stream_for_test<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    sandbox: &PathSandbox,
) -> Result<(), Box<dyn std::error::Error>> {
    serve_stream(reader, writer, sandbox)
}

fn serve_stream<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    sandbox: &PathSandbox,
) -> Result<(), Box<dyn std::error::Error>> {
    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let req: JsonRpcRequest = match serde_json::from_str(trimmed) {
            Ok(r) => r,
            Err(e) => {
                let err_resp = JsonRpcResponse {
                    jsonrpc: "2.0",
                    id: Value::Null,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32700,
                        message: format!("parse error: {e}"),
                    }),
                };
                write_json_line(writer, &err_resp)?;
                continue;
            }
        };

        // Notifications have no `id`
        let Some(id) = req.id.clone() else {
            continue;
        };

        let resp = handle_request(id, &req.method, &req.params, sandbox);
        write_json_line(writer, &resp)?;
    }
    Ok(())
}

fn write_json_line<W: Write, T: Serialize>(
    writer: &mut W,
    value: &T,
) -> Result<(), Box<dyn std::error::Error>> {
    let s = serde_json::to_string(value)?;
    writer.write_all(s.as_bytes())?;
    writer.write_all(b"\n")?;
    writer.flush()?;
    Ok(())
}

fn handle_request(
    id: Value,
    method: &str,
    params: &Value,
    sandbox: &PathSandbox,
) -> JsonRpcResponse {
    match method {
        "initialize" => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({
                "protocolVersion": "2025-06-18",
                "capabilities": {
                    "tools": {}
                },
                "serverInfo": {
                    "name": "admetcraft-mcp",
                    "version": env!("CARGO_PKG_VERSION")
                }
            })),
            error: None,
        },
        "ping" => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({})),
            error: None,
        },
        "tools/list" => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({ "tools": tool_definitions() })),
            error: None,
        },
        "tools/call" => {
            let tool_name = params
                .get("name")
                .and_then(Value::as_str)
                .unwrap_or("");
            let args = params.get("arguments").cloned().unwrap_or_else(|| json!({}));
            let result = call_tool(tool_name, &args, sandbox);
            JsonRpcResponse {
                jsonrpc: "2.0",
                id,
                result: Some(result),
                error: None,
            }
        }
        other => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: None,
            error: Some(JsonRpcError {
                code: -32601,
                message: format!("method not found: {other}"),
            }),
        },
    }
}

fn tool_definitions() -> Value {
    json!([
        {
            "name": "evaluate_smiles",
            "description": "Evaluate one or more SMILES strings (or a batch file inside --read-root) for physicochemical, ADME, druglikeness, and medicinal chemistry metrics.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": {
                        "description": "Single SMILES string or array of SMILES strings (each 'SMILES [Name]').",
                        "oneOf": [
                            { "type": "string" },
                            { "type": "array", "items": { "type": "string" } }
                        ]
                    },
                    "file": {
                        "type": "string",
                        "description": "Optional path to .smi/.csv/.sdf file inside --read-root."
                    }
                }
            }
        },
        {
            "name": "filter_druglikeness",
            "description": "Filter a list of SMILES compounds by Lipinski violations, BBB permeability, GI absorption, and PAINS/Brenk alerts.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "array",
                        "items": { "type": "string" },
                        "description": "List of SMILES strings (each 'SMILES [Name]')."
                    },
                    "max_lipinski_violations": { "type": "integer" },
                    "bbb_permeant": { "type": "boolean" },
                    "gi_high": { "type": "boolean" },
                    "no_pains": { "type": "boolean" },
                    "no_brenk": { "type": "boolean" }
                },
                "required": ["smiles"]
            }
        },
        {
            "name": "inspect_molecule",
            "description": "Return detailed structural breakdown (formula, heavy atoms, rings, aromaticity, 2D coordinates) and PAINS/Brenk alerts for a SMILES.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string" },
                    "name": { "type": "string" }
                },
                "required": ["smiles"]
            }
        },
        {
            "name": "export_report",
            "description": "Generate a CSV, JSON, PDF, or BOILED-Egg SVG report strictly inside --write-root.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "array",
                        "items": { "type": "string" }
                    },
                    "output_path": { "type": "string" },
                    "format": {
                        "type": "string",
                        "enum": ["csv", "json", "pdf", "svg"]
                    }
                },
                "required": ["smiles", "output_path", "format"]
            }
        }
    ])
}

fn tool_ok(payload: Value) -> Value {
    let text = serde_json::to_string_pretty(&payload).unwrap_or_else(|_| "{}".to_string());
    json!({
        "content": [{ "type": "text", "text": text }],
        "structuredContent": payload,
        "isError": false
    })
}

fn tool_err(msg: impl Into<String>) -> Value {
    let text = msg.into();
    json!({
        "content": [{ "type": "text", "text": text }],
        "isError": true
    })
}

fn call_tool(name: &str, args: &Value, sandbox: &PathSandbox) -> Value {
    match name {
        "evaluate_smiles" => tool_evaluate_smiles(args, sandbox),
        "filter_druglikeness" => tool_filter_druglikeness(args),
        "inspect_molecule" => tool_inspect_molecule(args),
        "export_report" => tool_export_report(args, sandbox),
        other => tool_err(format!("unknown tool: {other}")),
    }
}

fn collect_smiles_inputs(args: &Value) -> Vec<String> {
    let mut lines = Vec::new();
    if let Some(s) = args.get("smiles").and_then(Value::as_str) {
        for l in s.lines() {
            let t = l.trim();
            if !t.is_empty() {
                lines.push(t.to_string());
            }
        }
    } else if let Some(arr) = args.get("smiles").and_then(Value::as_array) {
        for item in arr {
            if let Some(s) = item.as_str() {
                let t = s.trim();
                if !t.is_empty() {
                    lines.push(t.to_string());
                }
            }
        }
    }
    lines
}

fn eval_lines(lines: &[String]) -> (Vec<MoleculeEvaluation>, Vec<Value>) {
    let joined = lines.join("\n");
    let parsed = parse_smiles_batch(&joined);
    let mut mols = Vec::new();
    let mut errors = Vec::new();
    for (line_no, raw, res) in parsed {
        match res {
            Ok(m) => mols.push(m),
            Err(e) => errors.push(json!({ "line": line_no, "input": raw, "error": e.to_string() })),
        }
    }
    (evaluate_batch(&mols), errors)
}

fn tool_evaluate_smiles(args: &Value, sandbox: &PathSandbox) -> Value {
    if let Some(file_str) = args.get("file").and_then(Value::as_str) {
        let resolved = match sandbox.resolve_read(Path::new(file_str)) {
            Ok(p) => p,
            Err(e) => return tool_err(format!("read sandbox error: {e}")),
        };
        let content = match std::fs::read_to_string(&resolved) {
            Ok(c) => c,
            Err(e) => return tool_err(format!("failed to read file: {e}")),
        };
        let ext = resolved
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        let mut mols = Vec::new();
        let mut errors = Vec::new();
        if ext == "sdf" || ext == "mol" {
            for (i, res) in parse_sdf(&content).into_iter().enumerate() {
                match res {
                    Ok(m) => mols.push(m),
                    Err(e) => errors.push(json!({ "record": i + 1, "error": e.to_string() })),
                }
            }
        } else if ext == "csv" {
            match parse_csv_batch(&content) {
                Ok(items) => {
                    for (i, res) in items.into_iter().enumerate() {
                        match res {
                            Ok(m) => mols.push(m),
                            Err(e) => errors.push(json!({ "row": i + 2, "error": e.to_string() })),
                        }
                    }
                }
                Err(e) => return tool_err(format!("CSV parse error: {e}")),
            }
        } else {
            for (line_no, raw, res) in parse_smiles_batch(&content) {
                match res {
                    Ok(m) => mols.push(m),
                    Err(e) => errors.push(json!({ "line": line_no, "input": raw, "error": e.to_string() })),
                }
            }
        }
        let evals = evaluate_batch(&mols);
        return tool_ok(json!({
            "count": evals.len(),
            "results": evals,
            "errors": errors,
        }));
    }

    let lines = collect_smiles_inputs(args);
    if lines.is_empty() {
        return tool_err("provide 'smiles' (string or array) or 'file' inside --read-root");
    }
    let (evals, errors) = eval_lines(&lines);
    tool_ok(json!({
        "count": evals.len(),
        "results": evals,
        "errors": errors,
    }))
}

fn tool_filter_druglikeness(args: &Value) -> Value {
    let lines = collect_smiles_inputs(args);
    if lines.is_empty() {
        return tool_err("'smiles' array must not be empty");
    }
    let criteria = FilterCriteria {
        max_lipinski_violations: args
            .get("max_lipinski_violations")
            .and_then(Value::as_u64)
            .map(|v| v as usize),
        require_bbb: args.get("bbb_permeant").and_then(Value::as_bool),
        require_gi_high: args.get("gi_high").and_then(Value::as_bool),
        no_pains: args.get("no_pains").and_then(Value::as_bool).unwrap_or(false),
        no_brenk: args.get("no_brenk").and_then(Value::as_bool).unwrap_or(false),
    };
    let (evals, errors) = eval_lines(&lines);
    let total = evals.len();
    let passed: Vec<MoleculeEvaluation> = evals
        .into_iter()
        .filter(|e| criteria.matches(e))
        .collect();
    tool_ok(json!({
        "total_evaluated": total,
        "matched_count": passed.len(),
        "matched": passed,
        "errors": errors,
    }))
}

fn tool_inspect_molecule(args: &Value) -> Value {
    let Some(smi) = args.get("smiles").and_then(Value::as_str) else {
        return tool_err("missing required 'smiles' string");
    };
    let name = args.get("name").and_then(Value::as_str);
    let mut mol = match parse_smiles(smi) {
        Ok(m) => m,
        Err(e) => return tool_err(format!("invalid SMILES: {e}")),
    };
    chem_core::generate_2d_coords(&mut mol);
    let eval = match evaluate_smiles(smi, name) {
        Ok(e) => e,
        Err(e) => return tool_err(format!("evaluation failed: {e}")),
    };
    tool_ok(json!({
        "name": eval.name,
        "canonical_smiles": eval.canonical_smiles,
        "formula": eval.physicochemical.formula,
        "heavy_atoms": mol.heavy_atom_count(),
        "aromatic_heavy_atoms": mol.aromatic_heavy_atom_count(),
        "rings": mol.rings,
        "atoms": mol.atoms,
        "bonds": mol.bonds,
        "pains_alerts": eval.medicinal_chemistry.pains_alerts,
        "brenk_alerts": eval.medicinal_chemistry.brenk_alerts,
        "synthetic_accessibility": eval.medicinal_chemistry.synthetic_accessibility,
    }))
}

fn tool_export_report(args: &Value, sandbox: &PathSandbox) -> Value {
    let lines = collect_smiles_inputs(args);
    if lines.is_empty() {
        return tool_err("'smiles' list must not be empty");
    }
    let Some(out_str) = args.get("output_path").and_then(Value::as_str) else {
        return tool_err("missing required 'output_path'");
    };
    let fmt_str = args
        .get("format")
        .and_then(Value::as_str)
        .unwrap_or("csv")
        .to_ascii_lowercase();
    let format = match fmt_str.as_str() {
        "csv" => ExportFormat::Csv,
        "json" => ExportFormat::Json,
        "pdf" => ExportFormat::Pdf,
        "svg" => ExportFormat::SvgBoiledEgg,
        other => return tool_err(format!("unsupported format '{other}'")),
    };

    let resolved = match sandbox.resolve_write(Path::new(out_str)) {
        Ok(p) => p,
        Err(e) => return tool_err(format!("write sandbox error: {e}")),
    };

    let (evals, errors) = eval_lines(&lines);
    if evals.is_empty() {
        return tool_err("no valid molecules to export");
    }
    if let Err(e) = write_report_file(&evals, &resolved, format) {
        return tool_err(format!("export error: {e}"));
    }

    tool_ok(json!({
        "written_path": resolved.display().to_string(),
        "format": fmt_str,
        "exported_count": evals.len(),
        "errors": errors,
    }))
}
