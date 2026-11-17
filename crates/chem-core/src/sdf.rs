use crate::error::ChemError;
use crate::mol::{Atom, BondOrder, BondStereo, Molecule};

/// Parse a single V2000 MOL block.
pub fn parse_mol_block(block: &str, default_name: &str) -> Result<Molecule, ChemError> {
    let lines: Vec<&str> = block.lines().collect();
    if lines.len() < 4 {
        return Err(ChemError::InvalidSdf {
            line: 1,
            msg: "MOL block has fewer than 4 header/counts lines".to_string(),
        });
    }

    let title = lines[0].trim();
    let counts_line = lines[3];
    if counts_line.len() < 6 {
        return Err(ChemError::InvalidSdf {
            line: 4,
            msg: "counts line too short".to_string(),
        });
    }

    let atom_count: usize = counts_line
        .get(0..3)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| ChemError::InvalidSdf {
            line: 4,
            msg: "invalid atom count".to_string(),
        })?;
    let bond_count: usize = counts_line
        .get(3..6)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| ChemError::InvalidSdf {
            line: 4,
            msg: "invalid bond count".to_string(),
        })?;

    if lines.len() < 4 + atom_count + bond_count {
        return Err(ChemError::InvalidSdf {
            line: 4,
            msg: "truncated atom or bond block".to_string(),
        });
    }

    let mut mol = Molecule::new();
    mol.name = Some(if title.is_empty() {
        default_name.to_string()
    } else {
        title.to_string()
    });

    for idx in 0..atom_count {
        let line_num = 4 + idx;
        let line = lines[line_num];
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 4 {
            return Err(ChemError::InvalidSdf {
                line: line_num + 1,
                msg: "atom line has fewer than 4 columns".to_string(),
            });
        }
        let x: f64 = tokens[0].parse().unwrap_or(0.0);
        let y: f64 = tokens[1].parse().unwrap_or(0.0);
        let symbol = tokens[3];
        let mut atom = Atom::new(idx, symbol)?;
        atom.coords_2d = Some([x, y]);
        // V2000 charge code in column 5 (1=+3, 2=+2, 3=+1, 5=-1, 6=-2, 7=-3)
        if let Some(chg_tok) = tokens.get(5)
            && let Ok(code) = chg_tok.parse::<i8>()
        {
            atom.formal_charge = match code {
                1 => 3,
                2 => 2,
                3 => 1,
                5 => -1,
                6 => -2,
                7 => -3,
                _ => 0,
            };
        }
        mol.add_atom(atom);
    }

    for idx in 0..bond_count {
        let line_num = 4 + atom_count + idx;
        let line = lines[line_num];
        // V2000 fixed-width or whitespace-separated first 3 fields
        let (u_raw, v_raw, t_raw) = if line.len() >= 9 {
            let s1 = line.get(0..3).unwrap_or("").trim();
            let s2 = line.get(3..6).unwrap_or("").trim();
            let s3 = line.get(6..9).unwrap_or("").trim();
            if let (Ok(u), Ok(v), Ok(t)) = (s1.parse::<usize>(), s2.parse::<usize>(), s3.parse::<u8>()) {
                (u, v, t)
            } else {
                parse_bond_tokens(line, line_num + 1)?
            }
        } else {
            parse_bond_tokens(line, line_num + 1)?
        };

        if u_raw == 0 || v_raw == 0 || u_raw > atom_count || v_raw > atom_count {
            return Err(ChemError::InvalidSdf {
                line: line_num + 1,
                msg: "bond atom index out of range".to_string(),
            });
        }
        let order = match t_raw {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            _ => BondOrder::Single,
        };
        mol.add_bond(u_raw - 1, v_raw - 1, order, BondStereo::None)?;
    }

    // Parse M  CHG lines and SDF property tags (> <NAME> or > <TITLE>)
    let mut prop_lines = (4 + atom_count + bond_count)..lines.len();
    while let Some(li) = prop_lines.next() {
        let line = lines[li];
        if let Some(rest) = line.strip_prefix("M  CHG") {
            let nums: Vec<i32> = rest
                .split_whitespace()
                .filter_map(|t| t.parse::<i32>().ok())
                .collect();
            if let Some(&count) = nums.first() {
                for pair in nums[1..].as_chunks::<2>().0.iter().take(count.max(0) as usize) {
                    let a_idx = (pair[0] - 1) as usize;
                    if let Some(atom) = mol.atoms.get_mut(a_idx) {
                        atom.formal_charge = pair[1].clamp(-15, 15) as i8;
                    }
                }
            }
        } else if (line.contains("<NAME>") || line.contains("<Name>") || line.contains("<ID>") || line.contains("<COMMON_NAME>"))
            && title.is_empty()
            && let Some(next_li) = prop_lines.next()
        {
            let val = lines[next_li].trim();
            if !val.is_empty() {
                mol.name = Some(val.to_string());
            }
        }
    }

    mol.suppress_explicit_hydrogens();
    mol.perceive()?;
    mol.input_smiles = Some(mol.to_canonical_smiles());
    Ok(mol)
}

fn parse_bond_tokens(line: &str, line_num: usize) -> Result<(usize, usize, u8), ChemError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(ChemError::InvalidSdf {
            line: line_num,
            msg: "bond line has fewer than 3 columns".to_string(),
        });
    }
    let u = tokens[0].parse().map_err(|_| ChemError::InvalidSdf {
        line: line_num,
        msg: "invalid bond atom 1".to_string(),
    })?;
    let v = tokens[1].parse().map_err(|_| ChemError::InvalidSdf {
        line: line_num,
        msg: "invalid bond atom 2".to_string(),
    })?;
    let t = tokens[2].parse().map_err(|_| ChemError::InvalidSdf {
        line: line_num,
        msg: "invalid bond order".to_string(),
    })?;
    Ok((u, v, t))
}

/// Parse a multi-record SDF file into a vector of `Result<Molecule, ChemError>`.
pub fn parse_sdf(content: &str) -> Vec<Result<Molecule, ChemError>> {
    content
        .split("$$$$")
        .map(str::trim)
        .filter(|rec| !rec.is_empty())
        .enumerate()
        .map(|(idx, rec)| parse_mol_block(rec, &format!("Molecule {}", idx + 1)))
        .collect()
}
