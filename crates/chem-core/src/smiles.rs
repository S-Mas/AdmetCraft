use std::collections::HashMap;

use crate::element::element_by_symbol;
use crate::error::ChemError;
use crate::mol::{Atom, BondOrder, BondStereo, ChiralTag, Molecule};

/// Default upper bound on atoms accepted from a single SMILES input (DoS protection).
pub const DEFAULT_MAX_ATOMS: usize = 2_000;

/// Parse a SMILES string (enforcing `DEFAULT_MAX_ATOMS`).
pub fn parse_smiles(input: &str) -> Result<Molecule, ChemError> {
    parse_smiles_with_limit(input, DEFAULT_MAX_ATOMS)
}

/// Parse a SMILES string with an explicit maximum atom limit.
pub fn parse_smiles_with_limit(input: &str, max_atoms: usize) -> Result<Molecule, ChemError> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(ChemError::EmptyInput);
    }

    let bytes = trimmed.as_bytes();
    let mut i = 0usize;
    let mut mol = Molecule::new();
    mol.input_smiles = Some(trimmed.to_string());

    let mut stack: Vec<Option<usize>> = Vec::new();
    let mut current: Option<usize> = None;
    let mut pending_bond: Option<BondOrder> = None;
    let mut pending_stereo = BondStereo::None;
    let mut saw_chiral_or_ez = false;
    let mut ring_opens: HashMap<u16, (usize, Option<BondOrder>, BondStereo)> = HashMap::new();

    while i < bytes.len() {
        let ch = bytes[i] as char;
        match ch {
            '(' => {
                stack.push(current);
                i += 1;
            }
            ')' => {
                current = stack.pop().ok_or(ChemError::UnmatchedParen { pos: i })?;
                pending_bond = None;
                pending_stereo = BondStereo::None;
                i += 1;
            }
            '.' => {
                current = None;
                pending_bond = None;
                pending_stereo = BondStereo::None;
                i += 1;
            }
            '-' => {
                pending_bond = Some(BondOrder::Single);
                i += 1;
            }
            '=' => {
                pending_bond = Some(BondOrder::Double);
                i += 1;
            }
            '#' => {
                pending_bond = Some(BondOrder::Triple);
                i += 1;
            }
            ':' => {
                pending_bond = Some(BondOrder::Aromatic);
                i += 1;
            }
            '/' => {
                pending_stereo = BondStereo::Up;
                saw_chiral_or_ez = true;
                if pending_bond.is_none() {
                    pending_bond = Some(BondOrder::Single);
                }
                i += 1;
            }
            '\\' => {
                pending_stereo = BondStereo::Down;
                saw_chiral_or_ez = true;
                if pending_bond.is_none() {
                    pending_bond = Some(BondOrder::Single);
                }
                i += 1;
            }
            '%' => {
                if i + 2 >= bytes.len() {
                    return Err(ChemError::InvalidSmiles {
                        pos: i,
                        msg: "incomplete '%' ring closure".to_string(),
                    });
                }
                let d1 = (bytes[i + 1] as char).to_digit(10);
                let d2 = (bytes[i + 2] as char).to_digit(10);
                let (Some(a), Some(b)) = (d1, d2) else {
                    return Err(ChemError::InvalidSmiles {
                        pos: i,
                        msg: "non-digit in '%' ring closure".to_string(),
                    });
                };
                let ring_num = (a * 10 + b) as u16;
                let anchor = current.ok_or_else(|| ChemError::InvalidSmiles {
                    pos: i,
                    msg: "ring closure with no preceding atom".to_string(),
                })?;
                handle_ring_closure(
                    &mut mol,
                    &mut ring_opens,
                    anchor,
                    ring_num,
                    pending_bond.take(),
                    std::mem::take(&mut pending_stereo),
                )?;
                i += 3;
            }
            '0'..='9' => {
                let ring_num = u16::from(bytes[i] - b'0');
                let anchor = current.ok_or_else(|| ChemError::InvalidSmiles {
                    pos: i,
                    msg: "ring closure with no preceding atom".to_string(),
                })?;
                handle_ring_closure(
                    &mut mol,
                    &mut ring_opens,
                    anchor,
                    ring_num,
                    pending_bond.take(),
                    std::mem::take(&mut pending_stereo),
                )?;
                i += 1;
            }
            '[' => {
                let close = trimmed[i..].find(']').map(|off| i + off).ok_or_else(|| {
                    ChemError::InvalidSmiles {
                        pos: i,
                        msg: "unclosed '[' bracket atom".to_string(),
                    }
                })?;
                let content = &trimmed[(i + 1)..close];
                let atom = parse_bracket_atom(content, i)?;
                if atom.chiral != ChiralTag::None {
                    saw_chiral_or_ez = true;
                }
                if mol.atoms.len() >= max_atoms {
                    return Err(ChemError::TooManyAtoms {
                        count: mol.atoms.len() + 1,
                        limit: max_atoms,
                    });
                }
                let new_id = mol.add_atom(atom);
                if let Some(prev) = current {
                    let order = resolve_bond_order(&mol, prev, new_id, pending_bond.take());
                    mol.add_bond(prev, new_id, order, std::mem::take(&mut pending_stereo))?;
                } else {
                    pending_bond = None;
                    pending_stereo = BondStereo::None;
                }
                current = Some(new_id);
                i = close + 1;
            }
            _ => {
                let (atom, consumed) = parse_organic_atom(trimmed, i)?;
                if mol.atoms.len() >= max_atoms {
                    return Err(ChemError::TooManyAtoms {
                        count: mol.atoms.len() + 1,
                        limit: max_atoms,
                    });
                }
                let new_id = mol.add_atom(atom);
                if let Some(prev) = current {
                    let order = resolve_bond_order(&mol, prev, new_id, pending_bond.take());
                    mol.add_bond(prev, new_id, order, std::mem::take(&mut pending_stereo))?;
                } else {
                    pending_bond = None;
                    pending_stereo = BondStereo::None;
                }
                current = Some(new_id);
                i += consumed;
            }
        }
    }

    if !stack.is_empty() {
        return Err(ChemError::UnmatchedParen { pos: trimmed.len() });
    }
    if let Some((&ring, _)) = ring_opens.iter().next() {
        return Err(ChemError::UnclosedRing { ring });
    }
    if mol.atoms.is_empty() {
        return Err(ChemError::EmptyInput);
    }

    if saw_chiral_or_ez {
        mol.warnings.push(
            "Stereochemistry markers (@/@@ or /\\) present in input; 3D/stereo descriptors are not modeled and canonical SMILES omits stereo."
                .to_string(),
        );
    }

    mol.suppress_explicit_hydrogens();
    mol.desalt();
    mol.perceive()?;
    Ok(mol)
}

/// Parse a line formatted as `SMILES [Name]` or `SMILES\tName`.
pub fn parse_smiles_line(line: &str, default_name: &str) -> Result<Molecule, ChemError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Err(ChemError::EmptyInput);
    }
    let mut parts = trimmed.splitn(2, char::is_whitespace);
    let smi = parts.next().unwrap_or("");
    let name = parts
        .next()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or(default_name);
    let mut mol = parse_smiles(smi)?;
    mol.name = Some(name.to_string());
    Ok(mol)
}

fn handle_ring_closure(
    mol: &mut Molecule,
    ring_opens: &mut HashMap<u16, (usize, Option<BondOrder>, BondStereo)>,
    anchor: usize,
    ring_num: u16,
    bond_opt: Option<BondOrder>,
    stereo: BondStereo,
) -> Result<(), ChemError> {
    if let Some((start_atom, start_bond, start_stereo)) = ring_opens.remove(&ring_num) {
        let chosen_bond = bond_opt.or(start_bond);
        let order = resolve_bond_order(mol, start_atom, anchor, chosen_bond);
        let chosen_stereo = if stereo != BondStereo::None {
            stereo
        } else {
            start_stereo
        };
        mol.add_bond(start_atom, anchor, order, chosen_stereo)?;
    } else {
        ring_opens.insert(ring_num, (anchor, bond_opt, stereo));
    }
    Ok(())
}

fn resolve_bond_order(
    mol: &Molecule,
    u: usize,
    v: usize,
    explicit: Option<BondOrder>,
) -> BondOrder {
    if let Some(bo) = explicit {
        return bo;
    }
    let au = mol.atoms.get(u).is_some_and(|a| a.aromatic);
    let av = mol.atoms.get(v).is_some_and(|a| a.aromatic);
    if au && av {
        BondOrder::Aromatic
    } else {
        BondOrder::Single
    }
}

fn parse_organic_atom(s: &str, pos: usize) -> Result<(Atom, usize), ChemError> {
    let rest = &s[pos..];
    if rest.starts_with("Cl") {
        return Ok((Atom::new(0, "Cl")?, 2));
    }
    if rest.starts_with("Br") {
        return Ok((Atom::new(0, "Br")?, 2));
    }
    let ch = rest.chars().next().ok_or(ChemError::EmptyInput)?;
    let (sym, aromatic) = match ch {
        'B' => ("B", false),
        'C' => ("C", false),
        'N' => ("N", false),
        'O' => ("O", false),
        'P' => ("P", false),
        'S' => ("S", false),
        'F' => ("F", false),
        'I' => ("I", false),
        'b' => ("B", true),
        'c' => ("C", true),
        'n' => ("N", true),
        'o' => ("O", true),
        'p' => ("P", true),
        's' => ("S", true),
        _ => {
            return Err(ChemError::InvalidSmiles {
                pos,
                msg: format!("unexpected character '{ch}'"),
            });
        }
    };
    let mut atom = Atom::new(0, sym)?;
    atom.aromatic = aromatic;
    Ok((atom, 1))
}

fn parse_bracket_atom(content: &str, pos: usize) -> Result<Atom, ChemError> {
    if content.is_empty() {
        return Err(ChemError::InvalidSmiles {
            pos,
            msg: "empty bracket atom".to_string(),
        });
    }
    let bytes = content.as_bytes();
    let mut idx = 0usize;

    // Optional isotope
    let iso_start = idx;
    while idx < bytes.len() && bytes[idx].is_ascii_digit() {
        idx += 1;
    }
    let isotope = if idx > iso_start {
        Some(
            content[iso_start..idx]
                .parse::<u16>()
                .map_err(|_| ChemError::InvalidSmiles {
                    pos,
                    msg: "invalid isotope number".to_string(),
                })?,
        )
    } else {
        None
    };

    if idx >= bytes.len() {
        return Err(ChemError::InvalidSmiles {
            pos,
            msg: "missing element in bracket atom".to_string(),
        });
    }

    // Element symbol (1 or 2 letters, lowercase first = aromatic)
    let c0 = bytes[idx] as char;
    if !c0.is_ascii_alphabetic() {
        return Err(ChemError::InvalidSmiles {
            pos,
            msg: format!("invalid element start '{c0}' in bracket atom"),
        });
    }
    let aromatic = c0.is_ascii_lowercase();
    let sym_end = if idx + 1 < bytes.len()
        && (bytes[idx + 1] as char).is_ascii_lowercase()
        && bytes[idx + 1] != b'h'
    {
        let two = format!("{}{}", c0.to_ascii_uppercase(), bytes[idx + 1] as char);
        if element_by_symbol(&two).is_some() {
            idx + 2
        } else {
            idx + 1
        }
    } else {
        idx + 1
    };

    let raw_sym = &content[idx..sym_end];
    let norm_sym = if raw_sym.len() == 1 {
        raw_sym.to_ascii_uppercase()
    } else {
        let mut chars = raw_sym.chars();
        let first = chars.next().map(|c| c.to_ascii_uppercase()).unwrap_or('?');
        let rest: String = chars.collect();
        format!("{first}{rest}")
    };
    idx = sym_end;

    let mut atom = Atom::new(0, &norm_sym)?;
    atom.isotope = isotope;
    atom.aromatic = aromatic;
    atom.from_bracket = true;

    // Optional chirality @ or @@
    if idx < bytes.len() && bytes[idx] == b'@' {
        if idx + 1 < bytes.len() && bytes[idx + 1] == b'@' {
            atom.chiral = ChiralTag::Clockwise;
            idx += 2;
        } else {
            atom.chiral = ChiralTag::CounterClockwise;
            idx += 1;
        }
    }

    // Optional H count
    if idx < bytes.len() && bytes[idx] == b'H' {
        idx += 1;
        let h_start = idx;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            idx += 1;
        }
        atom.explicit_h = if idx > h_start {
            content[h_start..idx]
                .parse::<u8>()
                .map_err(|_| ChemError::InvalidSmiles {
                    pos,
                    msg: "invalid explicit H count".to_string(),
                })?
        } else {
            1
        };
    }

    // Optional formal charge (+, -, +2, -1, ++, --), bounded strictly to [-15, +15]
    if idx < bytes.len() && (bytes[idx] == b'+' || bytes[idx] == b'-') {
        let sign: i16 = if bytes[idx] == b'+' { 1 } else { -1 };
        let ch_byte = bytes[idx];
        idx += 1;
        let charge_i16 = if idx < bytes.len() && bytes[idx].is_ascii_digit() {
            let d_start = idx;
            while idx < bytes.len() && bytes[idx].is_ascii_digit() {
                idx += 1;
            }
            let mag =
                content[d_start..idx]
                    .parse::<i16>()
                    .map_err(|_| ChemError::InvalidSmiles {
                        pos,
                        msg: "invalid formal charge".to_string(),
                    })?;
            sign * mag
        } else {
            let mut count: i16 = 1;
            while idx < bytes.len() && bytes[idx] == ch_byte {
                count = count.saturating_add(1);
                if count > 15 {
                    return Err(ChemError::InvalidSmiles {
                        pos,
                        msg: "formal charge exceeds allowed [-15, +15] range".to_string(),
                    });
                }
                idx += 1;
            }
            sign * count
        };
        if !(-15..=15).contains(&charge_i16) {
            return Err(ChemError::InvalidSmiles {
                pos,
                msg: "formal charge exceeds allowed [-15, +15] range".to_string(),
            });
        }
        atom.formal_charge = charge_i16 as i8;
    }

    // Optional atom class `:N`
    if idx < bytes.len() && bytes[idx] == b':' {
        idx += 1;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            idx += 1;
        }
    }

    if idx < bytes.len() {
        return Err(ChemError::InvalidSmiles {
            pos,
            msg: format!("unexpected trailing bracket content '{}'", &content[idx..]),
        });
    }

    Ok(atom)
}
