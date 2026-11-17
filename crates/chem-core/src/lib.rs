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
pub use smiles::{DEFAULT_MAX_ATOMS, parse_smiles, parse_smiles_line, parse_smiles_with_limit};

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

/// Parse CSV batch text (`smiles`/`Canonical SMILES` column + optional `name`/`Molecule`/`id` column).
/// Properly handles RFC-4180 quoted fields containing commas.
pub fn parse_csv_batch(input: &str) -> Result<Vec<Result<Molecule, ChemError>>, ChemError> {
    if input.trim().is_empty() {
        return Err(ChemError::EmptyInput);
    }
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_reader(input.as_bytes());
    let headers = rdr
        .headers()
        .map_err(|e| ChemError::Graph(format!("CSV header error: {e}")))?
        .clone();

    let smi_col = headers
        .iter()
        .position(|c| {
            let norm = c.trim().to_ascii_lowercase().replace(' ', "_");
            matches!(norm.as_str(), "smiles" | "canonical_smiles" | "smi")
        })
        .unwrap_or(0);
    let name_col = headers.iter().position(|c| {
        let norm = c.trim().to_ascii_lowercase().replace(' ', "_");
        matches!(norm.as_str(), "name" | "molecule" | "molecule_name" | "id")
    });

    let mut out = Vec::new();
    for (idx, rec_res) in rdr.records().enumerate() {
        let rec = match rec_res {
            Ok(r) => r,
            Err(e) => {
                out.push(Err(ChemError::Graph(format!("CSV row {}: {e}", idx + 2))));
                continue;
            }
        };
        let Some(smi) = rec.get(smi_col).map(str::trim).filter(|s| !s.is_empty()) else {
            continue;
        };
        let default_name = format!("Molecule {}", idx + 1);
        let name = name_col
            .and_then(|c| rec.get(c))
            .map(str::trim)
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
        let mut sildenafil =
            parse_smiles("CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12")?;
        generate_2d_coords(&mut sildenafil);
        for a in &sildenafil.atoms {
            let xy = a
                .coords_2d
                .ok_or_else(|| ChemError::Graph("missing 2d coords".to_string()))?;
            assert!(xy[0].is_finite() && xy[1].is_finite());
        }

        let mut aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        generate_2d_coords(&mut aspirin);
        for i in 0..aspirin.atoms.len() {
            for j in (i + 1)..aspirin.atoms.len() {
                if aspirin.bond_between(i, j).is_some() {
                    continue;
                }
                let ci = aspirin.atoms[i]
                    .coords_2d
                    .ok_or_else(|| ChemError::Graph("missing".to_string()))?;
                let cj = aspirin.atoms[j]
                    .coords_2d
                    .ok_or_else(|| ChemError::Graph("missing".to_string()))?;
                let d = ((ci[0] - cj[0]).powi(2) + (ci[1] - cj[1]).powi(2)).sqrt();
                assert!(d > 1.0, "aspirin atoms {i} and {j} overlap (d={d:.2})");
            }
        }

        // Norbornane bridged bicyclic: no elongated > 2.3 bond lengths
        let mut norbornane = parse_smiles("C1CC2CCC1C2")?;
        generate_2d_coords(&mut norbornane);
        for b in &norbornane.bonds {
            let cu = norbornane.atoms[b.u]
                .coords_2d
                .ok_or_else(|| ChemError::Graph("missing".to_string()))?;
            let cv = norbornane.atoms[b.v]
                .coords_2d
                .ok_or_else(|| ChemError::Graph("missing".to_string()))?;
            let len = ((cu[0] - cv[0]).powi(2) + (cu[1] - cv[1]).powi(2)).sqrt();
            assert!(
                len < 2.3,
                "norbornane bond {}-{} too long ({len:.2})",
                b.u,
                b.v
            );
        }
        Ok(())
    }

    #[test]
    fn desalt_strips_counterion_and_warns() -> Result<(), ChemError> {
        let sodium_acetate = parse_smiles("CC(=O)[O-].[Na+]")?;
        assert_eq!(sodium_acetate.formula(), "C2H3O2");
        assert!(
            sodium_acetate
                .warnings
                .iter()
                .any(|w| w.contains("Desalted"))
        );
        Ok(())
    }

    #[test]
    fn stereochemistry_emits_warning() -> Result<(), ChemError> {
        let l_ala = parse_smiles("C[C@@H](C(=O)O)N")?;
        assert!(l_ala.warnings.iter().any(|w| w.contains("Stereochemistry")));
        let trans_butene = parse_smiles("C/C=C/C")?;
        assert!(
            trans_butene
                .warnings
                .iter()
                .any(|w| w.contains("Stereochemistry"))
        );
        Ok(())
    }

    #[test]
    fn huge_chain_uses_iterative_dfs_and_limit_protects() -> Result<(), ChemError> {
        let giant = "C".repeat(28_000);
        let err = parse_smiles(&giant);
        assert!(matches!(err, Err(ChemError::TooManyAtoms { .. })));

        // Iterative DFS handles 3,000-atom chain without stack overflow
        let long_chain = "C".repeat(3_000);
        let mol = parse_smiles_with_limit(&long_chain, 4_000)?;
        let can = mol.to_canonical_smiles();
        assert_eq!(can.len(), 3_000);
        Ok(())
    }

    #[test]
    fn polycycle_over_100_rings_reuses_ring_numbers_and_reparses() -> Result<(), ChemError> {
        // Build a linear chain of 105 spiro/linked cyclopropyl rings (105 ring closures!)
        let smi = format!("C1CC1{}", "C2CC2".repeat(104));
        let mol = parse_smiles(&smi)?;
        assert_eq!(mol.rings.len(), 105);
        let can = mol.to_canonical_smiles();
        let reparsed = parse_smiles(&can)?;
        assert_eq!(reparsed.atoms.len(), mol.atoms.len());
        assert_eq!(reparsed.rings.len(), 105);
        Ok(())
    }

    #[test]
    fn csv_batch_handles_canonical_smiles_and_quoted_commas() -> Result<(), ChemError> {
        let csv_in = "Molecule,Canonical SMILES,MW\n\"Aspirin, acetylsalicylic\",CC(=O)Oc1ccccc1C(=O)O,180.16\n";
        let res = parse_csv_batch(csv_in)?;
        assert_eq!(res.len(), 1);
        let mol = res.into_iter().next().ok_or(ChemError::EmptyInput)??;
        assert_eq!(mol.name.as_deref(), Some("Aspirin, acetylsalicylic"));
        assert_eq!(mol.formula(), "C9H8O4");
        Ok(())
    }

    #[test]
    fn desalt_preserves_mixtures_cocrystals_and_inorganic_salt_pairs() -> Result<(), ChemError> {
        let benzene_dimer = parse_smiles("c1ccccc1.c1ccccc1")?;
        assert_eq!(benzene_dimer.formula(), "C12H12");
        assert!(
            !benzene_dimer
                .warnings
                .iter()
                .any(|w| w.contains("Desalted"))
        );

        let cocrystal = parse_smiles("CCO.CC(=O)O")?;
        assert_eq!(cocrystal.formula(), "C4H10O3");
        assert!(!cocrystal.warnings.iter().any(|w| w.contains("Desalted")));

        let nacl = parse_smiles("[Na+].[Cl-]")?;
        assert_eq!(nacl.formula(), "ClNa");
        assert!(!nacl.warnings.iter().any(|w| w.contains("Desalted")));

        let morphine_hcl_hydrate = parse_smiles("CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5.Cl.O")?;
        assert_eq!(morphine_hcl_hydrate.formula(), "C17H19NO3");
        assert!(
            morphine_hcl_hydrate
                .warnings
                .iter()
                .any(|w| w.contains("Desalted"))
        );

        let na_acetate = parse_smiles("CC(=O)[O-].[Na+]")?;
        assert_eq!(na_acetate.formula(), "C2H3O2");
        assert!(na_acetate.warnings.iter().any(|w| w.contains("Desalted")));
        Ok(())
    }

    #[test]
    fn formal_charge_overflow_and_bad_sdf_are_rejected_cleanly() {
        let plus_128 = format!("[C{}]", "+".repeat(128));
        assert!(matches!(
            parse_smiles(&plus_128),
            Err(ChemError::InvalidSmiles { .. })
        ));
        assert!(matches!(
            parse_smiles("[C+20]"),
            Err(ChemError::InvalidSmiles { .. })
        ));

        let bad_sdf_chg = "bad\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  CHG  1   0  -1\nM  END\n";
        let res_zero = parse_sdf(bad_sdf_chg);
        assert_eq!(res_zero.len(), 1);
        assert!(matches!(&res_zero[0], Ok(m) if m.atoms[0].formal_charge == 0));

        let bad_sdf_min = "bad\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  CHG  1 -2147483648 2147483647\nM  END\n";
        let res_min = parse_sdf(bad_sdf_min);
        assert_eq!(res_min.len(), 1);
        assert!(res_min[0].is_ok());

        let bad_sdf_giant = "bad\n  test\n\n999  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        assert!(matches!(
            parse_sdf(bad_sdf_giant).first(),
            Some(Err(ChemError::InvalidSdf { .. }))
        ));
    }

    #[test]
    fn k21_graph_reports_too_many_open_rings_error() -> Result<(), ChemError> {
        let mut mol = Molecule::new();
        for i in 0..21 {
            mol.add_atom(Atom::new(i, "C")?);
        }
        for i in 0..21 {
            for j in (i + 1)..21 {
                mol.add_bond(i, j, BondOrder::Single, BondStereo::None)?;
            }
        }
        assert!(matches!(
            mol.try_to_canonical_smiles(),
            Err(ChemError::TooManyOpenRings)
        ));
        Ok(())
    }
}
