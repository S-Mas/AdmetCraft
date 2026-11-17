use std::collections::HashMap;

use crate::element::element_by_symbol;
use crate::error::ChemError;
use crate::mol::{BondOrder, Molecule};

/// Compiled SMARTS substructure pattern.
#[derive(Debug, Clone)]
pub struct SmartsPattern {
    source: String,
    atoms: Vec<AtomExpr>,
    bonds: Vec<SmartsBond>,
    adj: Vec<Vec<(usize, usize)>>,
}

#[derive(Debug, Clone)]
struct SmartsBond {
    u: usize,
    v: usize,
    expr: BondExpr,
}

#[derive(Debug, Clone)]
enum AtomExpr {
    True,
    Not(Box<AtomExpr>),
    And(Vec<AtomExpr>),
    Or(Vec<AtomExpr>),
    AtomicNum(u8),
    AliphaticElement(u8),
    AromaticElement(u8),
    Aromatic(bool),
    TotalH(u8),
    HasImplicitH(Option<u8>),
    Degree(usize),
    Connectivity(usize),
    Valence(u8),
    Charge(i8),
    InRing(bool),
    RingCount(usize),
    RingSize(usize),
    Recursive(SmartsPattern),
}

#[derive(Debug, Clone)]
enum BondExpr {
    SingleOrAromatic,
    Any,
    Order(BondOrder),
    InRing(bool),
    Not(Box<BondExpr>),
    And(Vec<BondExpr>),
    Or(Vec<BondExpr>),
}

impl SmartsPattern {
    /// Parse and compile a SMARTS pattern string.
    pub fn parse(pattern: &str) -> Result<Self, ChemError> {
        let trimmed = pattern.trim();
        if trimmed.is_empty() {
            return Err(ChemError::InvalidSmarts {
                pattern: pattern.to_string(),
                msg: "empty pattern".to_string(),
            });
        }
        let mut parser = SmartsParser {
            src: trimmed,
            bytes: trimmed.as_bytes(),
            pos: 0,
            atoms: Vec::new(),
            bonds: Vec::new(),
        };
        parser.parse_pattern()?;
        let n = parser.atoms.len();
        let mut adj = vec![Vec::new(); n];
        for (bid, b) in parser.bonds.iter().enumerate() {
            adj[b.u].push((b.v, bid));
            adj[b.v].push((b.u, bid));
        }
        Ok(Self {
            source: trimmed.to_string(),
            atoms: parser.atoms,
            bonds: parser.bonds,
            adj,
        })
    }

    /// Return the original SMARTS pattern string.
    pub fn as_str(&self) -> &str {
        &self.source
    }

    /// Return true if the pattern matches at least once in `mol`.
    pub fn matches(&self, mol: &Molecule) -> bool {
        !self.find_matches(mol).is_empty()
    }

    /// Find all subgraph embeddings (each mapping pattern atom `i` -> molecule atom `map[i]`).
    pub fn find_matches(&self, mol: &Molecule) -> Vec<Vec<usize>> {
        let mut results = Vec::new();
        if self.atoms.is_empty() || mol.atoms.len() < self.atoms.len() {
            return results;
        }
        let mut mapping = vec![usize::MAX; self.atoms.len()];
        let mut used = vec![false; mol.atoms.len()];

        for mol_a in 0..mol.atoms.len() {
            if self.atom_matches(mol, mol_a, &self.atoms[0]) {
                mapping[0] = mol_a;
                used[mol_a] = true;
                self.backtrack(mol, 1, &mut mapping, &mut used, &mut results);
                used[mol_a] = false;
                mapping[0] = usize::MAX;
            }
        }
        results
    }

    /// Find unique matched atom sets (deduplicating symmetry-equivalent embeddings).
    pub fn count_unique_matches(&self, mol: &Molecule) -> usize {
        use std::collections::BTreeSet;
        let mut sets: BTreeSet<Vec<usize>> = BTreeSet::new();
        for m in self.find_matches(mol) {
            let mut s = m;
            s.sort_unstable();
            sets.insert(s);
        }
        sets.len()
    }

    fn backtrack(
        &self,
        mol: &Molecule,
        pat_idx: usize,
        mapping: &mut [usize],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
    ) {
        if pat_idx == self.atoms.len() {
            // Verify all pattern bonds (including ring closures) are satisfied
            for b in &self.bonds {
                let mu = mapping[b.u];
                let mv = mapping[b.v];
                let Some(mb) = mol.bond_between(mu, mv) else {
                    return;
                };
                if !bond_matches(mb.order, mb.in_ring, &b.expr) {
                    return;
                }
            }
            results.push(mapping.to_vec());
            return;
        }

        // Find an already-mapped neighbor of pat_idx in the pattern
        let anchor_pat = self.adj[pat_idx]
            .iter()
            .find(|&&(nbr, _)| nbr < pat_idx && mapping[nbr] != usize::MAX);

        if let Some(&(nbr_pat, bid)) = anchor_pat {
            let nbr_mol = mapping[nbr_pat];
            let bond_expr = &self.bonds[bid].expr;
            for &(cand_mol, mbid) in &mol.adj[nbr_mol] {
                if used[cand_mol] {
                    continue;
                }
                let mb = &mol.bonds[mbid];
                if !bond_matches(mb.order, mb.in_ring, bond_expr) {
                    continue;
                }
                if !self.atom_matches(mol, cand_mol, &self.atoms[pat_idx]) {
                    continue;
                }
                // Check bonds to any other already-mapped neighbors of pat_idx
                let ok_all = self.adj[pat_idx].iter().all(|&(other_pat, other_bid)| {
                    if other_pat >= pat_idx || mapping[other_pat] == usize::MAX {
                        return true;
                    }
                    let om = mapping[other_pat];
                    mol.bond_between(cand_mol, om)
                        .is_some_and(|mb2| bond_matches(mb2.order, mb2.in_ring, &self.bonds[other_bid].expr))
                });
                if !ok_all {
                    continue;
                }
                mapping[pat_idx] = cand_mol;
                used[cand_mol] = true;
                self.backtrack(mol, pat_idx + 1, mapping, used, results);
                used[cand_mol] = false;
                mapping[pat_idx] = usize::MAX;
            }
        } else {
            for cand_mol in 0..mol.atoms.len() {
                if used[cand_mol] {
                    continue;
                }
                if self.atom_matches(mol, cand_mol, &self.atoms[pat_idx]) {
                    mapping[pat_idx] = cand_mol;
                    used[cand_mol] = true;
                    self.backtrack(mol, pat_idx + 1, mapping, used, results);
                    used[cand_mol] = false;
                    mapping[pat_idx] = usize::MAX;
                }
            }
        }
    }

    fn atom_matches(&self, mol: &Molecule, mol_a: usize, expr: &AtomExpr) -> bool {
        let Some(atom) = mol.atoms.get(mol_a) else {
            return false;
        };
        match expr {
            AtomExpr::True => true,
            AtomExpr::Not(inner) => !self.atom_matches(mol, mol_a, inner),
            AtomExpr::And(items) => items.iter().all(|e| self.atom_matches(mol, mol_a, e)),
            AtomExpr::Or(items) => items.iter().any(|e| self.atom_matches(mol, mol_a, e)),
            AtomExpr::AtomicNum(z) => atom.atomic_number == *z,
            AtomExpr::AliphaticElement(z) => atom.atomic_number == *z && !atom.aromatic,
            AtomExpr::AromaticElement(z) => atom.atomic_number == *z && atom.aromatic,
            AtomExpr::Aromatic(flag) => atom.aromatic == *flag,
            AtomExpr::TotalH(h) => atom.total_h() == *h,
            AtomExpr::HasImplicitH(opt) => match opt {
                Some(n) => atom.implicit_h == *n,
                None => atom.implicit_h > 0,
            },
            AtomExpr::Degree(d) => mol.degree(mol_a) == *d,
            AtomExpr::Connectivity(x) => mol.total_degree(mol_a) == *x,
            AtomExpr::Valence(v) => {
                let bos = mol.bond_order_sum(mol_a) + f64::from(atom.total_h());
                (bos.round() as u8) == *v
            }
            AtomExpr::Charge(c) => atom.formal_charge == *c,
            AtomExpr::InRing(r) => atom.in_ring == *r,
            AtomExpr::RingCount(cnt) => atom.ring_sizes.len() == *cnt,
            AtomExpr::RingSize(sz) => atom.ring_sizes.contains(sz),
            AtomExpr::Recursive(sub) => {
                if sub.atoms.is_empty() {
                    return false;
                }
                if !sub.atom_matches(mol, mol_a, &sub.atoms[0]) {
                    return false;
                }
                let mut mapping = vec![usize::MAX; sub.atoms.len()];
                let mut used = vec![false; mol.atoms.len()];
                mapping[0] = mol_a;
                used[mol_a] = true;
                let mut found = Vec::new();
                sub.backtrack_first(mol, 1, &mut mapping, &mut used, &mut found)
            }
        }
    }

    fn backtrack_first(
        &self,
        mol: &Molecule,
        pat_idx: usize,
        mapping: &mut [usize],
        used: &mut [bool],
        _dummy: &mut Vec<usize>,
    ) -> bool {
        if pat_idx == self.atoms.len() {
            return self.bonds.iter().all(|b| {
                mol.bond_between(mapping[b.u], mapping[b.v])
                    .is_some_and(|mb| bond_matches(mb.order, mb.in_ring, &b.expr))
            });
        }
        let anchor_pat = self.adj[pat_idx]
            .iter()
            .find(|&&(nbr, _)| nbr < pat_idx && mapping[nbr] != usize::MAX);

        if let Some(&(nbr_pat, bid)) = anchor_pat {
            let nbr_mol = mapping[nbr_pat];
            let bond_expr = &self.bonds[bid].expr;
            for &(cand_mol, mbid) in &mol.adj[nbr_mol] {
                if used[cand_mol] {
                    continue;
                }
                let mb = &mol.bonds[mbid];
                if !bond_matches(mb.order, mb.in_ring, bond_expr) {
                    continue;
                }
                if !self.atom_matches(mol, cand_mol, &self.atoms[pat_idx]) {
                    continue;
                }
                let ok_all = self.adj[pat_idx].iter().all(|&(other_pat, other_bid)| {
                    if other_pat >= pat_idx || mapping[other_pat] == usize::MAX {
                        return true;
                    }
                    let om = mapping[other_pat];
                    mol.bond_between(cand_mol, om)
                        .is_some_and(|mb2| bond_matches(mb2.order, mb2.in_ring, &self.bonds[other_bid].expr))
                });
                if !ok_all {
                    continue;
                }
                mapping[pat_idx] = cand_mol;
                used[cand_mol] = true;
                let hit = self.backtrack_first(mol, pat_idx + 1, mapping, used, _dummy);
                used[cand_mol] = false;
                mapping[pat_idx] = usize::MAX;
                if hit {
                    return true;
                }
            }
            false
        } else {
            for cand_mol in 0..mol.atoms.len() {
                if used[cand_mol] {
                    continue;
                }
                if self.atom_matches(mol, cand_mol, &self.atoms[pat_idx]) {
                    mapping[pat_idx] = cand_mol;
                    used[cand_mol] = true;
                    let hit = self.backtrack_first(mol, pat_idx + 1, mapping, used, _dummy);
                    used[cand_mol] = false;
                    mapping[pat_idx] = usize::MAX;
                    if hit {
                        return true;
                    }
                }
            }
            false
        }
    }
}

fn bond_matches(order: BondOrder, in_ring: bool, expr: &BondExpr) -> bool {
    match expr {
        BondExpr::SingleOrAromatic => matches!(order, BondOrder::Single | BondOrder::Aromatic),
        BondExpr::Any => true,
        BondExpr::Order(o) => order == *o,
        BondExpr::InRing(r) => in_ring == *r,
        BondExpr::Not(inner) => !bond_matches(order, in_ring, inner),
        BondExpr::And(items) => items.iter().all(|e| bond_matches(order, in_ring, e)),
        BondExpr::Or(items) => items.iter().any(|e| bond_matches(order, in_ring, e)),
    }
}

struct SmartsParser<'a> {
    src: &'a str,
    bytes: &'a [u8],
    pos: usize,
    atoms: Vec<AtomExpr>,
    bonds: Vec<SmartsBond>,
}

impl<'a> SmartsParser<'a> {
    fn err(&self, msg: impl Into<String>) -> ChemError {
        ChemError::InvalidSmarts {
            pattern: self.src.to_string(),
            msg: msg.into(),
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn parse_pattern(&mut self) -> Result<(), ChemError> {
        let mut stack: Vec<Option<usize>> = Vec::new();
        let mut current: Option<usize> = None;
        let mut pending_bond: Option<BondExpr> = None;
        let mut ring_opens: HashMap<u16, (usize, Option<BondExpr>)> = HashMap::new();

        while let Some(b) = self.peek() {
            match b {
                b'(' => {
                    stack.push(current);
                    self.pos += 1;
                }
                b')' => {
                    current = stack.pop().ok_or_else(|| self.err("unmatched ')'"))?;
                    pending_bond = None;
                    self.pos += 1;
                }
                b'.' => {
                    current = None;
                    pending_bond = None;
                    self.pos += 1;
                }
                b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'!' | b'/' | b'\\' => {
                    pending_bond = Some(self.parse_bond_expr()?);
                }
                b'%' => {
                    self.pos += 1;
                    if self.pos + 1 >= self.bytes.len() {
                        return Err(self.err("incomplete '%' ring closure"));
                    }
                    let d1 = (self.bytes[self.pos] as char).to_digit(10);
                    let d2 = (self.bytes[self.pos + 1] as char).to_digit(10);
                    let (Some(a), Some(c)) = (d1, d2) else {
                        return Err(self.err("invalid '%' ring closure"));
                    };
                    self.pos += 2;
                    let rnum = (a * 10 + c) as u16;
                    let anchor = current.ok_or_else(|| self.err("ring closure without atom"))?;
                    self.handle_ring(&mut ring_opens, anchor, rnum, pending_bond.take());
                }
                b'0'..=b'9' => {
                    let rnum = u16::from(b - b'0');
                    self.pos += 1;
                    let anchor = current.ok_or_else(|| self.err("ring closure without atom"))?;
                    self.handle_ring(&mut ring_opens, anchor, rnum, pending_bond.take());
                }
                b'[' => {
                    let expr = self.parse_bracket_atom()?;
                    let idx = self.atoms.len();
                    self.atoms.push(expr);
                    if let Some(prev) = current {
                        let bexpr = pending_bond.take().unwrap_or(BondExpr::SingleOrAromatic);
                        self.bonds.push(SmartsBond { u: prev, v: idx, expr: bexpr });
                    } else {
                        pending_bond = None;
                    }
                    current = Some(idx);
                }
                _ => {
                    let expr = self.parse_bare_atom()?;
                    let idx = self.atoms.len();
                    self.atoms.push(expr);
                    if let Some(prev) = current {
                        let bexpr = pending_bond.take().unwrap_or(BondExpr::SingleOrAromatic);
                        self.bonds.push(SmartsBond { u: prev, v: idx, expr: bexpr });
                    } else {
                        pending_bond = None;
                    }
                    current = Some(idx);
                }
            }
        }
        if !stack.is_empty() {
            return Err(self.err("unclosed '('"));
        }
        if !ring_opens.is_empty() {
            return Err(self.err("unclosed ring closure in SMARTS"));
        }
        Ok(())
    }

    fn handle_ring(
        &mut self,
        opens: &mut HashMap<u16, (usize, Option<BondExpr>)>,
        anchor: usize,
        rnum: u16,
        bexpr: Option<BondExpr>,
    ) {
        if let Some((start, start_bond)) = opens.remove(&rnum) {
            let expr = bexpr.or(start_bond).unwrap_or(BondExpr::SingleOrAromatic);
            self.bonds.push(SmartsBond { u: start, v: anchor, expr });
        } else {
            opens.insert(rnum, (anchor, bexpr));
        }
    }

    fn parse_bond_expr(&mut self) -> Result<BondExpr, ChemError> {
        let mut or_terms = Vec::new();
        loop {
            let mut and_terms = Vec::new();
            while let Some(b) = self.peek() {
                let negated = if b == b'!' {
                    self.pos += 1;
                    true
                } else {
                    false
                };
                let Some(tok) = self.peek() else {
                    return Err(self.err("trailing '!' in bond expression"));
                };
                let atom_b = match tok {
                    b'-' => {
                        self.pos += 1;
                        BondExpr::Order(BondOrder::Single)
                    }
                    b'=' => {
                        self.pos += 1;
                        BondExpr::Order(BondOrder::Double)
                    }
                    b'#' => {
                        self.pos += 1;
                        BondExpr::Order(BondOrder::Triple)
                    }
                    b':' => {
                        self.pos += 1;
                        BondExpr::Order(BondOrder::Aromatic)
                    }
                    b'~' => {
                        self.pos += 1;
                        BondExpr::Any
                    }
                    b'@' => {
                        self.pos += 1;
                        BondExpr::InRing(true)
                    }
                    b'/' | b'\\' => {
                        self.pos += 1;
                        BondExpr::Order(BondOrder::Single)
                    }
                    _ => break,
                };
                and_terms.push(if negated {
                    BondExpr::Not(Box::new(atom_b))
                } else {
                    atom_b
                });
                if self.peek() == Some(b'&') {
                    self.pos += 1;
                } else {
                    break;
                }
            }
            if and_terms.is_empty() {
                break;
            }
            let combined = if and_terms.len() == 1 {
                and_terms.remove(0)
            } else {
                BondExpr::And(and_terms)
            };
            or_terms.push(combined);
            if self.peek() == Some(b',') {
                self.pos += 1;
            } else {
                break;
            }
        }
        if or_terms.is_empty() {
            return Err(self.err("invalid bond expression"));
        }
        if or_terms.len() == 1 {
            Ok(or_terms.remove(0))
        } else {
            Ok(BondExpr::Or(or_terms))
        }
    }

    fn parse_bare_atom(&mut self) -> Result<AtomExpr, ChemError> {
        let rest = &self.src[self.pos..];
        if rest.starts_with("Cl") {
            self.pos += 2;
            return Ok(AtomExpr::AliphaticElement(17));
        }
        if rest.starts_with("Br") {
            self.pos += 2;
            return Ok(AtomExpr::AliphaticElement(35));
        }
        let b = self.peek().ok_or_else(|| self.err("unexpected end of SMARTS"))?;
        self.pos += 1;
        match b {
            b'*' => Ok(AtomExpr::True),
            b'A' => Ok(AtomExpr::Aromatic(false)),
            b'a' => Ok(AtomExpr::Aromatic(true)),
            b'B' => Ok(AtomExpr::AliphaticElement(5)),
            b'C' => Ok(AtomExpr::AliphaticElement(6)),
            b'N' => Ok(AtomExpr::AliphaticElement(7)),
            b'O' => Ok(AtomExpr::AliphaticElement(8)),
            b'P' => Ok(AtomExpr::AliphaticElement(15)),
            b'S' => Ok(AtomExpr::AliphaticElement(16)),
            b'F' => Ok(AtomExpr::AliphaticElement(9)),
            b'I' => Ok(AtomExpr::AliphaticElement(53)),
            b'b' => Ok(AtomExpr::AromaticElement(5)),
            b'c' => Ok(AtomExpr::AromaticElement(6)),
            b'n' => Ok(AtomExpr::AromaticElement(7)),
            b'o' => Ok(AtomExpr::AromaticElement(8)),
            b'p' => Ok(AtomExpr::AromaticElement(15)),
            b's' => Ok(AtomExpr::AromaticElement(16)),
            _ => Err(self.err(format!("unexpected atom token '{}'", b as char))),
        }
    }

    fn parse_bracket_atom(&mut self) -> Result<AtomExpr, ChemError> {
        self.pos += 1; // skip '['
        let expr = self.parse_low_and()?;
        if self.peek() != Some(b']') {
            return Err(self.err("expected ']' closing atom expression"));
        }
        self.pos += 1;
        Ok(expr)
    }

    fn parse_low_and(&mut self) -> Result<AtomExpr, ChemError> {
        let mut items = vec![self.parse_or()?];
        while self.peek() == Some(b';') {
            self.pos += 1;
            items.push(self.parse_or()?);
        }
        Ok(if items.len() == 1 {
            items.remove(0)
        } else {
            AtomExpr::And(items)
        })
    }

    fn parse_or(&mut self) -> Result<AtomExpr, ChemError> {
        let mut items = vec![self.parse_high_and()?];
        while self.peek() == Some(b',') {
            self.pos += 1;
            items.push(self.parse_high_and()?);
        }
        Ok(if items.len() == 1 {
            items.remove(0)
        } else {
            AtomExpr::Or(items)
        })
    }

    fn parse_high_and(&mut self) -> Result<AtomExpr, ChemError> {
        let mut items = Vec::new();
        while let Some(b) = self.peek() {
            if matches!(b, b']' | b';' | b',') {
                break;
            }
            if b == b'&' {
                self.pos += 1;
                continue;
            }
            items.push(self.parse_unary()?);
        }
        if items.is_empty() {
            return Err(self.err("empty bracket atom subexpression"));
        }
        Ok(if items.len() == 1 {
            items.remove(0)
        } else {
            AtomExpr::And(items)
        })
    }

    fn parse_unary(&mut self) -> Result<AtomExpr, ChemError> {
        if self.peek() == Some(b'!') {
            self.pos += 1;
            let inner = self.parse_primitive()?;
            return Ok(AtomExpr::Not(Box::new(inner)));
        }
        self.parse_primitive()
    }

    fn parse_uint(&mut self) -> Option<usize> {
        let start = self.pos;
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        if self.pos == start {
            None
        } else {
            self.src[start..self.pos].parse().ok()
        }
    }

    fn parse_primitive(&mut self) -> Result<AtomExpr, ChemError> {
        let Some(b) = self.peek() else {
            return Err(self.err("unexpected end inside '['"));
        };
        match b {
            b'*' => {
                self.pos += 1;
                Ok(AtomExpr::True)
            }
            b'A' => {
                self.pos += 1;
                Ok(AtomExpr::Aromatic(false))
            }
            b'a' => {
                self.pos += 1;
                Ok(AtomExpr::Aromatic(true))
            }
            b'#' => {
                self.pos += 1;
                let z = self
                    .parse_uint()
                    .ok_or_else(|| self.err("expected atomic number after '#'"))?;
                Ok(AtomExpr::AtomicNum(z as u8))
            }
            b'H' => {
                self.pos += 1;
                let n = self.parse_uint().unwrap_or(1) as u8;
                Ok(AtomExpr::TotalH(n))
            }
            b'h' => {
                self.pos += 1;
                let n = self.parse_uint().map(|v| v as u8);
                Ok(AtomExpr::HasImplicitH(n))
            }
            b'D' => {
                self.pos += 1;
                let n = self.parse_uint().unwrap_or(1);
                Ok(AtomExpr::Degree(n))
            }
            b'X' => {
                self.pos += 1;
                let n = self.parse_uint().unwrap_or(1);
                Ok(AtomExpr::Connectivity(n))
            }
            b'v' => {
                self.pos += 1;
                let n = self.parse_uint().unwrap_or(1) as u8;
                Ok(AtomExpr::Valence(n))
            }
            b'R' => {
                self.pos += 1;
                if let Some(n) = self.parse_uint() {
                    if n == 0 {
                        Ok(AtomExpr::InRing(false))
                    } else {
                        Ok(AtomExpr::RingCount(n))
                    }
                } else {
                    Ok(AtomExpr::InRing(true))
                }
            }
            b'r' => {
                self.pos += 1;
                if let Some(n) = self.parse_uint() {
                    Ok(AtomExpr::RingSize(n))
                } else {
                    Ok(AtomExpr::InRing(true))
                }
            }
            b'+' => {
                self.pos += 1;
                if let Some(n) = self.parse_uint() {
                    Ok(AtomExpr::Charge(n as i8))
                } else {
                    let mut cnt = 1i8;
                    while self.peek() == Some(b'+') {
                        self.pos += 1;
                        cnt += 1;
                    }
                    Ok(AtomExpr::Charge(cnt))
                }
            }
            b'-' => {
                self.pos += 1;
                if let Some(n) = self.parse_uint() {
                    Ok(AtomExpr::Charge(-(n as i8)))
                } else {
                    let mut cnt = 1i8;
                    while self.peek() == Some(b'-') {
                        self.pos += 1;
                        cnt += 1;
                    }
                    Ok(AtomExpr::Charge(-cnt))
                }
            }
            b'$' => {
                self.pos += 1;
                if self.peek() != Some(b'(') {
                    return Err(self.err("expected '(' after '$'"));
                }
                self.pos += 1;
                let start = self.pos;
                let mut depth = 1usize;
                while self.pos < self.bytes.len() && depth > 0 {
                    match self.bytes[self.pos] {
                        b'(' => depth += 1,
                        b')' => depth -= 1,
                        _ => {}
                    }
                    if depth > 0 {
                        self.pos += 1;
                    }
                }
                if depth != 0 {
                    return Err(self.err("unclosed recursive SMARTS '$('"));
                }
                let inner_src = &self.src[start..self.pos];
                self.pos += 1; // skip ')'
                let compiled = SmartsPattern::parse(inner_src)?;
                Ok(AtomExpr::Recursive(compiled))
            }
            _ => {
                // Element symbol inside brackets (e.g. C, c, N, n, Cl, Br, Si, Se, etc.)
                let c0 = b as char;
                if !c0.is_ascii_alphabetic() {
                    return Err(self.err(format!("unexpected primitive '{c0}'")));
                }
                let is_arom = c0.is_ascii_lowercase();
                let end = if self.pos + 1 < self.bytes.len()
                    && (self.bytes[self.pos + 1] as char).is_ascii_lowercase()
                    && !matches!(self.bytes[self.pos + 1], b'h' | b'r' | b'v')
                {
                    let two = format!(
                        "{}{}",
                        c0.to_ascii_uppercase(),
                        self.bytes[self.pos + 1] as char
                    );
                    if element_by_symbol(&two).is_some() {
                        self.pos + 2
                    } else {
                        self.pos + 1
                    }
                } else {
                    self.pos + 1
                };
                let sym = &self.src[self.pos..end];
                self.pos = end;
                let norm = if sym.len() == 1 {
                    sym.to_ascii_uppercase()
                } else {
                    let mut chars = sym.chars();
                    let first = chars.next().map(|c| c.to_ascii_uppercase()).unwrap_or('?');
                    format!("{first}{}", chars.as_str())
                };
                let info = element_by_symbol(&norm)
                    .ok_or_else(|| self.err(format!("unknown element '{norm}'")))?;
                if is_arom {
                    Ok(AtomExpr::AromaticElement(info.atomic_number))
                } else {
                    Ok(AtomExpr::AliphaticElement(info.atomic_number))
                }
            }
        }
    }
}
