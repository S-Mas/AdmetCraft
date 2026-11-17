//! Layer 2 (`admet-report`): Vector charts (Bioavailability Radar, BOILED-Egg) and
//! CSV, JSON, and publication-ready PDF report generation.

#![forbid(unsafe_code)]

pub mod csv_json;
pub mod pdf;
pub mod svg;

use std::io::Write;
use std::path::Path;

use admet_calc::MoleculeEvaluation;
use thiserror::Error;

pub use csv_json::{export_csv_string, export_json_string};
pub use pdf::export_pdf_bytes;
pub use svg::{render_boiled_egg_svg, render_radar_svg};

/// Errors produced by `admet-report`.
#[derive(Debug, Error)]
pub enum ReportError {
    #[error("cannot generate report for empty molecule list")]
    EmptyReport,
    #[error("CSV serialization error: {0}")]
    Csv(#[from] csv::Error),
    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("formatting error: {0}")]
    Format(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Supported export formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    Csv,
    Json,
    Pdf,
    SvgBoiledEgg,
}

/// Write an evaluation report directly into an open `Write` stream (eliminating path TOCTOU).
pub fn write_report_to_writer<W: Write>(
    evals: &[MoleculeEvaluation],
    writer: &mut W,
    format: ExportFormat,
) -> Result<(), ReportError> {
    if evals.is_empty() {
        return Err(ReportError::EmptyReport);
    }
    match format {
        ExportFormat::Csv => {
            let text = export_csv_string(evals)?;
            writer.write_all(text.as_bytes())?;
        }
        ExportFormat::Json => {
            let text = export_json_string(evals)?;
            writer.write_all(text.as_bytes())?;
        }
        ExportFormat::Pdf => {
            let bytes = export_pdf_bytes(evals)?;
            writer.write_all(&bytes)?;
        }
        ExportFormat::SvgBoiledEgg => {
            let text = render_boiled_egg_svg(evals);
            writer.write_all(text.as_bytes())?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Write an evaluation report to `path` in the requested `format`.
pub fn write_report_file(
    evals: &[MoleculeEvaluation],
    path: &Path,
    format: ExportFormat,
) -> Result<(), ReportError> {
    if evals.is_empty() {
        return Err(ReportError::EmptyReport);
    }
    let mut f = std::fs::File::create(path)?;
    write_report_to_writer(evals, &mut f, format)
}

#[cfg(test)]
mod tests {
    use super::*;
    use admet_calc::evaluate_smiles;
    use chem_core::parse_csv_batch;

    #[test]
    fn exports_csv_json_svg_and_valid_pdf() -> Result<(), Box<dyn std::error::Error>> {
        let evals = vec![
            evaluate_smiles("CC(=O)Oc1ccccc1C(=O)O", Some("Aspirina, ácido"))?,
            evaluate_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C", Some("Cafeína"))?,
        ];
        let csv_str = export_csv_string(&evals)?;
        assert!(csv_str.contains("Aspirina, ácido"));
        assert!(csv_str.contains("Consensus Log P"));

        // Verify exported CSV roundtrips cleanly back through chem_core::parse_csv_batch!
        let reimported = parse_csv_batch(&csv_str)?;
        assert_eq!(reimported.len(), 2);
        let first = reimported
            .into_iter()
            .next()
            .ok_or(ReportError::EmptyReport)??;
        assert_eq!(first.name.as_deref(), Some("Aspirina, ácido"));

        let json_str = export_json_string(&evals)?;
        assert!(json_str.contains("bioavailability_radar"));

        let radar = render_radar_svg(&evals[0]);
        assert!(radar.contains("<svg") && radar.contains("LIPO"));

        let egg = render_boiled_egg_svg(&evals);
        assert!(egg.contains("<svg") && egg.contains("BBB (Yolk)"));

        let pdf_bytes = export_pdf_bytes(&evals)?;
        assert!(pdf_bytes.starts_with(b"%PDF-1.4"));
        assert!(pdf_bytes.ends_with(b"%%EOF\n"));
        let pdf_text = String::from_utf8_lossy(&pdf_bytes);
        assert!(pdf_text.contains("/Encoding /WinAnsiEncoding"));
        // Cafeína uses octal \355 for í instead of '?'
        assert!(pdf_text.contains("Cafe\\355na"));
        Ok(())
    }

    #[test]
    fn pdf_paginates_more_than_18_molecules() -> Result<(), Box<dyn std::error::Error>> {
        let mut evals = Vec::new();
        for i in 1..=30 {
            evals.push(evaluate_smiles("CCO", Some(&format!("Mol_{i}")))?);
        }
        let pdf_bytes = export_pdf_bytes(&evals)?;
        let pdf_text = String::from_utf8_lossy(&pdf_bytes);
        // Continuation summary page must contain Mol_25..Mol_30 (PDF escapes parens)!
        assert!(pdf_text.contains("Summary Table \\(continued: 25..30 of 30\\)"));
        assert!(pdf_text.contains("Mol_30"));
        Ok(())
    }
}
