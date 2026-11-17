use std::io::{BufRead, BufReader, Write};
use std::net::{SocketAddr, TcpListener};
use std::path::Path;

use admet_calc::{MoleculeEvaluation, evaluate_batch, evaluate_smiles};
use admet_engine::{FilterCriteria, MAX_READ_BYTES, PathSandbox};
use admet_report::{ExportFormat, write_report_to_writer};
use chem_core::{parse_csv_batch, parse_sdf, parse_smiles, parse_smiles_batch};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

/// Maximum incoming JSON-RPC request line size (1 MiB) to prevent memory amplification DoS.
pub const MAX_MCP_MESSAGE_BYTES: usize = 1024 * 1024;
/// Maximum number of molecules accepted in a single MCP tool call.
pub const MAX_MCP_MOLECULES: usize = 1_000;

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

struct SessionAuth {
    initialized: bool,
    required_token: Option<String>,
}

/// Run the MCP server over `stdin`/`stdout` (no token required for local stdio pipe,
/// though `initialize` lifecycle is still enforced).
pub fn run_mcp_stdio(sandbox: &PathSandbox) -> Result<(), Box<dyn std::error::Error>> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let reader = BufReader::new(stdin.lock());
    let mut writer = stdout.lock();
    serve_stream(reader, &mut writer, sandbox, None)
}

/// Run the MCP server over a TCP bridge address.
/// Enforces loopback-only binding unless `allow_remote` is true, and requires `auth_token`.
pub fn run_mcp_tcp(
    addr_str: &str,
    sandbox: &PathSandbox,
    allow_remote: bool,
    auth_token: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let addr: SocketAddr = addr_str.parse()?;
    if !addr.ip().is_loopback() && !allow_remote {
        return Err(format!(
            "refusing to bind non-loopback TCP address {addr} without --allow-remote"
        )
        .into());
    }
    let token = auth_token
        .map(|t| t.trim().to_string())
        .filter(|t| !t.is_empty())
        .ok_or("TCP MCP mode requires a non-empty --mcp-token (or ADMETCRAFT_MCP_TOKEN)")?;

    let listener = TcpListener::bind(addr)?;
    eprintln!("AdmetCraft MCP TCP listening on {addr} (auth required)");
    for stream_res in listener.incoming() {
        let Ok(stream) = stream_res else {
            continue;
        };
        let Ok(reader_stream) = stream.try_clone() else {
            continue;
        };
        let reader = BufReader::new(reader_stream);
        let mut writer = stream;
        let _ = serve_stream(reader, &mut writer, sandbox, Some(token.clone()));
    }
    Ok(())
}

#[cfg(test)]
pub fn serve_stream_for_test<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    sandbox: &PathSandbox,
    required_token: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    serve_stream(reader, writer, sandbox, required_token)
}

enum LineRead {
    Eof,
    Line(String),
    TooLarge,
}

fn read_limited_line<R: BufRead>(reader: &mut R, limit: usize) -> std::io::Result<LineRead> {
    let mut buf = Vec::new();
    let mut exceeded = false;
    loop {
        let avail = reader.fill_buf()?;
        if avail.is_empty() {
            if exceeded {
                return Ok(LineRead::TooLarge);
            }
            if buf.is_empty() {
                return Ok(LineRead::Eof);
            }
            return Ok(LineRead::Line(String::from_utf8_lossy(&buf).into_owned()));
        }
        if let Some(pos) = avail.iter().position(|&b| b == b'\n') {
            if !exceeded {
                if buf.len() + pos > limit {
                    exceeded = true;
                } else {
                    buf.extend_from_slice(&avail[..pos]);
                }
            }
            reader.consume(pos + 1);
            if exceeded {
                return Ok(LineRead::TooLarge);
            }
            return Ok(LineRead::Line(String::from_utf8_lossy(&buf).into_owned()));
        }
        let len = avail.len();
        if !exceeded {
            if buf.len() + len > limit {
                exceeded = true;
                buf.clear();
            } else {
                buf.extend_from_slice(avail);
            }
        }
        reader.consume(len);
    }
}

fn serve_stream<R: BufRead, W: Write>(
    mut reader: R,
    writer: &mut W,
    sandbox: &PathSandbox,
    required_token: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut auth = SessionAuth {
        initialized: false,
        required_token,
    };

    loop {
        match read_limited_line(&mut reader, MAX_MCP_MESSAGE_BYTES)? {
            LineRead::Eof => break,
            LineRead::TooLarge => {
                let err_resp = JsonRpcResponse {
                    jsonrpc: "2.0",
                    id: Value::Null,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32600,
                        message: format!(
                            "request exceeds maximum message size ({MAX_MCP_MESSAGE_BYTES} bytes)"
                        ),
                    }),
                };
                write_json_line(writer, &err_resp)?;
                continue;
            }
            LineRead::Line(line) => {
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

                let Some(id) = req.id.clone() else {
                    continue;
                };

                let resp = handle_request(id, &req.method, &req.params, sandbox, &mut auth);
                write_json_line(writer, &resp)?;
            }
        }
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
    auth: &mut SessionAuth,
) -> JsonRpcResponse {
    if method == "initialize" {
        if let Some(expected) = auth.required_token.as_deref() {
            let provided = params.get("token").and_then(Value::as_str).unwrap_or("");
            if !constant_time_eq(provided.as_bytes(), expected.as_bytes()) {
                return JsonRpcResponse {
                    jsonrpc: "2.0",
                    id,
                    result: None,
                    error: Some(JsonRpcError {
                        code: -32001,
                        message: "unauthorized: invalid or missing MCP auth token".to_string(),
                    }),
                };
            }
        }
        auth.initialized = true;
        return JsonRpcResponse {
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
        };
    }

    if method == "ping" {
        return JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({})),
            error: None,
        };
    }

    if !auth.initialized {
        return JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: None,
            error: Some(JsonRpcError {
                code: -32002,
                message: "session not initialized; call 'initialize' first".to_string(),
            }),
        };
    }

    match method {
        "tools/list" => JsonRpcResponse {
            jsonrpc: "2.0",
            id,
            result: Some(json!({ "tools": tool_definitions() })),
            error: None,
        },
        "tools/call" => {
            let tool_name = params.get("name").and_then(Value::as_str).unwrap_or("");
            let args = params
                .get("arguments")
                .cloned()
                .unwrap_or_else(|| json!({}));
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

fn constant_time_eq(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut diff = 0u8;
    for (&x, &y) in a.iter().zip(b.iter()) {
        diff |= x ^ y;
    }
    diff == 0
}

fn tool_definitions() -> Value {
    json!([
        {
            "name": "evaluate_smiles",
            "description": "Evaluate up to 1000 SMILES strings (or a regular file <=10 MiB inside --read-root) for physicochemical, ADME, druglikeness, and medicinal chemistry metrics.",
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
                        "description": "Optional path to .smi/.csv/.sdf regular file inside --read-root."
                    }
                }
            }
        },
        {
            "name": "filter_druglikeness",
            "description": "Filter up to 1000 SMILES compounds by Lipinski violations, BBB permeability, GI absorption, and PAINS/Brenk alerts.",
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
            "description": "Atomically generate a CSV, JSON, PDF, or BOILED-Egg SVG report strictly inside --write-root using O_NOFOLLOW.",
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

/// Build MCP success response without duplicating large JSON trees twice on the wire.
fn tool_ok(summary: String, payload: Value) -> Value {
    json!({
        "content": [{ "type": "text", "text": summary }],
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

fn collect_smiles_inputs(args: &Value) -> Result<Vec<String>, String> {
    let mut lines = Vec::new();
    if let Some(s) = args.get("smiles").and_then(Value::as_str) {
        for l in s.lines() {
            let t = l.trim();
            if !t.is_empty() {
                lines.push(t.to_string());
                if lines.len() > MAX_MCP_MOLECULES {
                    return Err(format!(
                        "batch exceeds maximum allowed {MAX_MCP_MOLECULES} molecules per MCP call"
                    ));
                }
            }
        }
    } else if let Some(arr) = args.get("smiles").and_then(Value::as_array) {
        for item in arr {
            if let Some(s) = item.as_str() {
                let t = s.trim();
                if !t.is_empty() {
                    lines.push(t.to_string());
                    if lines.len() > MAX_MCP_MOLECULES {
                        return Err(format!(
                            "batch exceeds maximum allowed {MAX_MCP_MOLECULES} molecules per MCP call"
                        ));
                    }
                }
            }
        }
    }
    Ok(lines)
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
        let (resolved, content) =
            match sandbox.read_to_string_limited(Path::new(file_str), MAX_READ_BYTES) {
                Ok(res) => res,
                Err(e) => return tool_err(format!("read sandbox error: {e}")),
            };
        let ext = resolved
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        let mut mols = Vec::new();
        let mut errors = Vec::new();
        if ext == "sdf" || ext == "mol" {
            for (i, res) in parse_sdf(&content)
                .into_iter()
                .take(MAX_MCP_MOLECULES)
                .enumerate()
            {
                match res {
                    Ok(m) => mols.push(m),
                    Err(e) => errors.push(json!({ "record": i + 1, "error": e.to_string() })),
                }
            }
        } else if ext == "csv" {
            match parse_csv_batch(&content) {
                Ok(items) => {
                    for (i, res) in items.into_iter().take(MAX_MCP_MOLECULES).enumerate() {
                        match res {
                            Ok(m) => mols.push(m),
                            Err(e) => errors.push(json!({ "row": i + 2, "error": e.to_string() })),
                        }
                    }
                }
                Err(e) => return tool_err(format!("CSV parse error: {e}")),
            }
        } else {
            for (line_no, raw, res) in parse_smiles_batch(&content)
                .into_iter()
                .take(MAX_MCP_MOLECULES)
            {
                match res {
                    Ok(m) => mols.push(m),
                    Err(e) => errors
                        .push(json!({ "line": line_no, "input": raw, "error": e.to_string() })),
                }
            }
        }
        let evals = evaluate_batch(&mols);
        let summary = format!(
            "Evaluated {} molecule(s) from {}",
            evals.len(),
            resolved.display()
        );
        return tool_ok(
            summary,
            json!({
                "count": evals.len(),
                "results": evals,
                "errors": errors,
            }),
        );
    }

    let lines = match collect_smiles_inputs(args) {
        Ok(l) => l,
        Err(msg) => return tool_err(msg),
    };
    if lines.is_empty() {
        return tool_err("provide 'smiles' (string or array) or 'file' inside --read-root");
    }
    let (evals, errors) = eval_lines(&lines);
    let summary = format!(
        "Evaluated {} molecule(s) ({} error(s))",
        evals.len(),
        errors.len()
    );
    tool_ok(
        summary,
        json!({
            "count": evals.len(),
            "results": evals,
            "errors": errors,
        }),
    )
}

fn tool_filter_druglikeness(args: &Value) -> Value {
    let lines = match collect_smiles_inputs(args) {
        Ok(l) => l,
        Err(msg) => return tool_err(msg),
    };
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
        no_pains: args
            .get("no_pains")
            .and_then(Value::as_bool)
            .unwrap_or(false),
        no_brenk: args
            .get("no_brenk")
            .and_then(Value::as_bool)
            .unwrap_or(false),
    };
    let (evals, errors) = eval_lines(&lines);
    let total = evals.len();
    let passed: Vec<MoleculeEvaluation> =
        evals.into_iter().filter(|e| criteria.matches(e)).collect();
    let summary = format!("Matched {}/{} molecule(s)", passed.len(), total);
    tool_ok(
        summary,
        json!({
            "total_evaluated": total,
            "matched_count": passed.len(),
            "matched": passed,
            "errors": errors,
        }),
    )
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
    let summary = format!("Inspected {} ({})", eval.name, eval.physicochemical.formula);
    tool_ok(
        summary,
        json!({
            "name": eval.name,
            "canonical_smiles": eval.canonical_smiles,
            "formula": eval.physicochemical.formula,
            "warnings": eval.warnings,
            "heavy_atoms": mol.heavy_atom_count(),
            "aromatic_heavy_atoms": mol.aromatic_heavy_atom_count(),
            "rings": mol.rings,
            "atoms": mol.atoms,
            "bonds": mol.bonds,
            "pains_alerts": eval.medicinal_chemistry.pains_alerts,
            "pains_catalog_size": eval.medicinal_chemistry.pains_catalog_size,
            "brenk_alerts": eval.medicinal_chemistry.brenk_alerts,
            "brenk_catalog_size": eval.medicinal_chemistry.brenk_catalog_size,
            "synthetic_accessibility": eval.medicinal_chemistry.synthetic_accessibility,
        }),
    )
}

fn tool_export_report(args: &Value, sandbox: &PathSandbox) -> Value {
    let lines = match collect_smiles_inputs(args) {
        Ok(l) => l,
        Err(msg) => return tool_err(msg),
    };
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

    let (evals, errors) = eval_lines(&lines);
    if evals.is_empty() {
        return tool_err("no valid molecules to export");
    }

    // Atomically validate and open with O_NOFOLLOW, writing directly to the fd (no TOCTOU)
    let (resolved, mut file) = match sandbox.open_write(Path::new(out_str)) {
        Ok(pair) => pair,
        Err(e) => return tool_err(format!("write sandbox error: {e}")),
    };
    if let Err(e) = write_report_to_writer(&evals, &mut file, format) {
        return tool_err(format!("export error: {e}"));
    }

    let summary = format!(
        "Exported {} molecule(s) to {}",
        evals.len(),
        resolved.display()
    );
    tool_ok(
        summary,
        json!({
            "written_path": resolved.display().to_string(),
            "format": fmt_str,
            "exported_count": evals.len(),
            "errors": errors,
        }),
    )
}
