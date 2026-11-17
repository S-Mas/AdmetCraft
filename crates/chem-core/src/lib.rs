//! Layer 0 (`chem-core`): SMILES & SDF/MOL parsers, molecular graph, ring/aromaticity perception,
//! 2D coordinate layout, and SMARTS substructure matching.

#![forbid(unsafe_code)]

pub mod element;
pub mod error;
pub mod layout;
pub mod mol;
pub mod sdf;
pub mod smarts;
pub mod smiles;

pub use element::{ElementInfo, element_by_number, element_by_symbol};
pub use error::ChemError;
pub use layout::{ensure_2d_coords, generate_2d_coords};
pub use mol::{Atom, Bond, BondOrder, BondStereo, ChiralTag, Hybridization, Molecule};
pub use sdf::{parse_mol_block, parse_sdf};
pub use smarts::SmartsPattern;
pub use smiles::{parse_smiles, parse_smiles_line};

/// Parse batch lines of `SMILES [Name]` (skipping blank lines and `#` comments).
pub fn parse_smiles_batch(input: &str) -> Vec<(usize, String, Result<Molecule, ChemError>)> {
    input
        .lines()
        .enumerate()
        .filter_map(|(line_idx, raw)| {
            let trimmed = raw.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                return None;
            }
            let default_name = format!("Molecule {}", line_idx + 1);
            let parsed = parse_smiles_line(trimmed, &default_name);
            Some((line_idx + 1, trimmed.to_string(), parsed))
        })
        .collect()
}

/// Parse CSV batch text (`smiles` column + optional `name`/`id` column).
pub fn parse_csv_batch(input: &str) -> Result<Vec<Result<Molecule, ChemError>>, ChemError> {
    let mut lines = input.lines().filter(|l| !l.trim().is_empty());
    let Some(header) = lines.next() else {
        return Err(ChemError::EmptyInput);
    };
    let cols: Vec<&str> = header.split(',').map(|s| s.trim().trim_matches('"')).collect();
    let smi_col = cols
        .iter()
        .position(|c| c.eq_ignore_ascii_case("smiles") || c.eq_ignore_ascii_case("canonical_smiles") || c.eq_ignore_ascii_case("smi"))
        .unwrap_or(0);
    let name_col = cols
        .iter()
        .position(|c| c.eq_ignore_ascii_case("name") || c.eq_ignore_ascii_case("molecule") || c.eq_ignore_ascii_case("id"));

    let mut out = Vec::new();
    for (idx, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split(',').map(|s| s.trim().trim_matches('"')).collect();
        let Some(smi) = fields.get(smi_col) else {
            continue;
        };
        let default_name = format!("Molecule {}", idx + 1);
        let name = name_col
            .and_then(|c| fields.get(c))
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .unwrap_or(&default_name);
        let res = parse_smiles(smi).map(|mut m| {
            m.name = Some(name.to_string());
            m
        });
        out.push(res);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_aspirin_and_caffeine() -> Result<(), ChemError> {
        let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        assert_eq!(aspirin.formula(), "C9H8O4");
        assert!((aspirin.molecular_weight() - 180.16).abs() < 0.05);
        assert_eq!(aspirin.heavy_atom_count(), 13);
        assert_eq!(aspirin.aromatic_heavy_atom_count(), 6);

        let caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")?;
        assert_eq!(caffeine.formula(), "C8H10N4O2");
        assert!((caffeine.molecular_weight() - 194.19).abs() < 0.05);
        assert_eq!(caffeine.heavy_atom_count(), 14);
        assert_eq!(caffeine.aromatic_heavy_atom_count(), 9);
        Ok(())
    }

    #[test]
    fn kekule_and_aromatic_benzene_match() -> Result<(), ChemError> {
        let kek = parse_smiles("C1=CC=CC=C1")?;
        let aro = parse_smiles("c1ccccc1")?;
        assert_eq!(kek.formula(), aro.formula());
        assert_eq!(kek.aromatic_heavy_atom_count(), 6);
        assert_eq!(aro.aromatic_heavy_atom_count(), 6);
        assert_eq!(kek.to_canonical_smiles(), aro.to_canonical_smiles());
        Ok(())
    }

    #[test]
    fn smarts_matches_functional_groups() -> Result<(), ChemError> {
        let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        let acid = SmartsPattern::parse("[CX3](=O)[OX2H1]")?;
        assert!(acid.matches(&aspirin));
        let nitro = SmartsPattern::parse("[NX3](=O)=O")?;
        assert!(!nitro.matches(&aspirin));
        Ok(())
    }

    #[test]
    fn sdf_v2000_parses() -> Result<(), ChemError> {
        let mol_block = "\
Ethanol
  AdmetCraft

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0
    2.2500    1.2990    0.0000 O   0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
";
        let res = parse_sdf(mol_block);
        assert_eq!(res.len(), 1);
        let mol = res.into_iter().next().ok_or(ChemError::EmptyInput)??;
        assert_eq!(mol.name.as_deref(), Some("Ethanol"));
        assert_eq!(mol.formula(), "C2H6O");
        Ok(())
    }

    #[test]
    fn layout_generates_finite_2d_coords() -> Result<(), ChemError> {
        let mut sildenafil = parse_smiles("CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12")?;
        generate_2d_coords(&mut sildenafil);
        for a in &sildenafil.atoms {
            let xy = a.coords_2d.ok_or_else(|| ChemError::Graph("missing 2d coords".to_string()))?;
            assert!(xy[0].is_finite() && xy[1].is_finite());
        }

        let mut aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        generate_2d_coords(&mut aspirin);
        for i in 0..aspirin.atoms.len() {
            for j in (i + 1)..aspirin.atoms.len() {
                if aspirin.bond_between(i, j).is_some() {
                    continue;
                }
                let ci = aspirin.atoms[i].coords_2d.ok_or_else(|| ChemError::Graph("missing".to_string()))?;
                let cj = aspirin.atoms[j].coords_2d.ok_or_else(|| ChemError::Graph("missing".to_string()))?;
                let d = ((ci[0] - cj[0]).powi(2) + (ci[1] - cj[1]).powi(2)).sqrt();
                assert!(d > 1.0, "aspirin atoms {i} and {j} overlap (d={d:.2})");
            }
        }
        Ok(())
    }
}
