use std::collections::{BTreeMap, BTreeSet, VecDeque};

use serde::{Deserialize, Serialize};

use crate::element::{element_by_number, element_by_symbol};
use crate::error::ChemError;

/// Covalent bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondOrder {
    /// Single bond (`-`).
    Single,
    /// Double bond (`=`).
    Double,
    /// Triple bond (`#`).
    Triple,
    /// Aromatic bond (`:`).
    Aromatic,
}

impl BondOrder {
    /// Numeric bond order contribution for valence accounting.
    pub fn valence_units(self) -> f64 {
        match self {
            Self::Single => 1.0,
            Self::Double => 2.0,
            Self::Triple => 3.0,
            Self::Aromatic => 1.5,
        }
    }

    /// Integer bond order contribution for non-aromatic valence.
    pub fn int_order(self) -> u8 {
        match self {
            Self::Single => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Aromatic => 1,
        }
    }
}

/// Bond stereochemistry marker.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum BondStereo {
    #[default]
    None,
    Up,
    Down,
    Cis,
    Trans,
}

/// Atom tetrahedral chirality marker.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ChiralTag {
    #[default]
    None,
    CounterClockwise,
    Clockwise,
}

/// Perceived atomic hybridization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum Hybridization {
    #[default]
    Other,
    Sp,
    Sp2,
    Sp3,
}

/// Single atom node in a molecular graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// Zero-based index inside `Molecule::atoms`.
    pub id: usize,
    /// Atomic number (Z).
    pub atomic_number: u8,
    /// Chemical symbol (e.g., `"C"`).
    pub symbol: String,
    /// Formal charge (e.g., `-1`, `0`, `+1`).
    pub formal_charge: i8,
    /// Optional isotopic mass number.
    pub isotope: Option<u16>,
    /// True if atom participates in an aromatic ring.
    pub aromatic: bool,
    /// Explicit hydrogen count from bracket SMILES or explicit H nodes.
    pub explicit_h: u8,
    /// Implicit hydrogen count computed from valence.
    pub implicit_h: u8,
    /// Tetrahedral chirality tag.
    pub chiral: ChiralTag,
    /// Perceived hybridization.
    pub hybridization: Hybridization,
    /// Whether this atom belongs to at least one ring.
    pub in_ring: bool,
    /// Sizes of minimal rings containing this atom.
    pub ring_sizes: Vec<usize>,
    /// Optional 2D coordinates in Angstrom-like units.
    pub coords_2d: Option<[f64; 2]>,
    /// True if parsed from a SMILES bracket `[...]`.
    pub from_bracket: bool,
}

impl Atom {
    /// Construct a neutral atom with default attributes.
    pub fn new(id: usize, symbol: &str) -> Result<Self, ChemError> {
        let info = element_by_symbol(symbol).ok_or_else(|| ChemError::UnknownElement {
            symbol: symbol.to_string(),
        })?;
        Ok(Self {
            id,
            atomic_number: info.atomic_number,
            symbol: info.symbol.to_string(),
            formal_charge: 0,
            isotope: None,
            aromatic: false,
            explicit_h: 0,
            implicit_h: 0,
            chiral: ChiralTag::None,
            hybridization: Hybridization::Other,
            in_ring: false,
            ring_sizes: Vec::new(),
            coords_2d: None,
            from_bracket: false,
        })
    }

    /// Total hydrogens attached to this atom (implicit + explicit).
    pub fn total_h(&self) -> u8 {
        self.implicit_h.saturating_add(self.explicit_h)
    }

    /// Standard IUPAC atomic weight of this atom (excluding attached implicit H).
    pub fn atomic_weight(&self) -> f64 {
        if let Some(iso) = self.isotope {
            return f64::from(iso);
        }
        element_by_number(self.atomic_number)
            .map(|e| e.atomic_weight)
            .unwrap_or(0.0)
    }
}

/// Edge in a molecular graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    /// Zero-based index inside `Molecule::bonds`.
    pub id: usize,
    /// Endpoint atom index 1.
    pub u: usize,
    /// Endpoint atom index 2.
    pub v: usize,
    /// Covalent bond order.
    pub order: BondOrder,
    /// Optional stereo indicator (`/` or `\`).
    pub stereo: BondStereo,
    /// True if edge belongs to at least one cycle.
    pub in_ring: bool,
}

/// Complete molecular graph representation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    /// Optional molecule name/identifier.
    pub name: Option<String>,
    /// Input SMILES (if created from SMILES).
    pub input_smiles: Option<String>,
    /// Atom list (heavy atoms, plus explicit H only if standalone or bracketed).
    pub atoms: Vec<Atom>,
    /// Bond list.
    pub bonds: Vec<Bond>,
    /// Adjacency list: `adj[u]` holds `(neighbor_v, bond_id)`.
    pub adj: Vec<Vec<(usize, usize)>>,
    /// Minimal rings (each ring is an ordered list of atom indices).
    pub rings: Vec<Vec<usize>>,
    /// Non-fatal warnings emitted during parsing, desalting, or standardization.
    #[serde(default)]
    pub warnings: Vec<String>,
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl Molecule {
    /// Create an empty molecule.
    pub fn new() -> Self {
        Self {
            name: None,
            input_smiles: None,
            atoms: Vec::new(),
            bonds: Vec::new(),
            adj: Vec::new(),
            rings: Vec::new(),
            warnings: Vec::new(),
        }
    }

    /// Add an atom and return its index.
    pub fn add_atom(&mut self, mut atom: Atom) -> usize {
        let id = self.atoms.len();
        atom.id = id;
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        id
    }

    /// Add a bond between `u` and `v`. If a bond already exists, returns its id.
    pub fn add_bond(
        &mut self,
        u: usize,
        v: usize,
        order: BondOrder,
        stereo: BondStereo,
    ) -> Result<usize, ChemError> {
        if u >= self.atoms.len() || v >= self.atoms.len() {
            return Err(ChemError::Graph(format!(
                "bond atom index out of bounds: ({u}, {v})"
            )));
        }
        if u == v {
            return Err(ChemError::Graph(
                "self-loop bond is not allowed".to_string(),
            ));
        }
        if let Some((_, existing_id)) = self.adj[u].iter().find(|&&(nbr, _)| nbr == v) {
            return Ok(*existing_id);
        }
        let id = self.bonds.len();
        self.bonds.push(Bond {
            id,
            u,
            v,
            order,
            stereo,
            in_ring: false,
        });
        self.adj[u].push((v, id));
        self.adj[v].push((u, id));
        Ok(id)
    }

    /// Find bond between two atoms if present.
    pub fn bond_between(&self, u: usize, v: usize) -> Option<&Bond> {
        let adj_u = self.adj.get(u)?;
        let (_, bid) = adj_u.iter().find(|&&(nbr, _)| nbr == v)?;
        self.bonds.get(*bid)
    }

    /// Number of non-hydrogen atoms.
    pub fn heavy_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| a.atomic_number != 1).count()
    }

    /// Number of aromatic non-hydrogen atoms.
    pub fn aromatic_heavy_atom_count(&self) -> usize {
        self.atoms
            .iter()
            .filter(|a| a.atomic_number != 1 && a.aromatic)
            .count()
    }

    /// Neighbor atom indices for `u`.
    pub fn neighbors(&self, u: usize) -> impl Iterator<Item = usize> + '_ {
        self.adj
            .get(u)
            .into_iter()
            .flat_map(|v| v.iter().map(|&(nbr, _)| nbr))
    }

    /// Heavy-atom degree (number of bonded heavy atoms).
    pub fn heavy_degree(&self, u: usize) -> usize {
        self.neighbors(u)
            .filter(|&nbr| self.atoms.get(nbr).is_some_and(|a| a.atomic_number != 1))
            .count()
    }

    /// Total explicit degree (bonded atoms in graph, excluding implicit H).
    pub fn degree(&self, u: usize) -> usize {
        self.adj.get(u).map_or(0, Vec::len)
    }

    /// Total connectivity including implicit and explicit hydrogens.
    pub fn total_degree(&self, u: usize) -> usize {
        let h_imp = self.atoms.get(u).map_or(0, |a| usize::from(a.implicit_h));
        let h_exp = self.atoms.get(u).map_or(0, |a| usize::from(a.explicit_h));
        self.degree(u) + h_imp + h_exp
    }

    /// Sum of bond orders on graph edges incident to `u` (aromatic = 1.5).
    pub fn bond_order_sum(&self, u: usize) -> f64 {
        let Some(edges) = self.adj.get(u) else {
            return 0.0;
        };
        edges
            .iter()
            .filter_map(|&(_, bid)| self.bonds.get(bid))
            .map(|b| b.order.valence_units())
            .sum()
    }

    /// Recompute rings, aromaticity, implicit H, and hybridization.
    pub fn perceive(&mut self) -> Result<(), ChemError> {
        self.perceive_rings();
        self.perceive_aromaticity();
        self.recompute_implicit_h();
        self.perceive_hybridization();
        Ok(())
    }

    /// Strip recognized counterions, inorganic salts, and hydrates from multi-component
    /// structures when at least one non-salt parent fragment remains.
    /// Neutral mixtures and cocrystals (e.g. `c1ccccc1.c1ccccc1`, `CCO.CC(=O)O`) and
    /// purely inorganic salt pairs (e.g. `[Na+].[Cl-]`) are preserved untouched.
    pub fn desalt(&mut self) {
        let n = self.atoms.len();
        if n <= 1 {
            return;
        }
        let mut seen = vec![false; n];
        let mut comps: Vec<Vec<usize>> = Vec::new();
        for start in 0..n {
            if seen[start] {
                continue;
            }
            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            seen[start] = true;
            q.push_back(start);
            while let Some(u) = q.pop_front() {
                comp.push(u);
                for &(v, _) in &self.adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        q.push_back(v);
                    }
                }
            }
            comps.push(comp);
        }
        if comps.len() <= 1 {
            return;
        }

        let is_salt: Vec<bool> = comps
            .iter()
            .map(|c| self.is_recognized_counterion(c))
            .collect();

        // Keep everything unless there is both at least one recognized counterion/salt
        // and at least one non-salt parent fragment.
        if !is_salt.iter().any(|&s| s) || !is_salt.iter().any(|&s| !s) {
            return;
        }

        let mut stripped_labels = Vec::new();
        let mut keep_set: BTreeSet<usize> = BTreeSet::new();
        for (idx, comp) in comps.iter().enumerate() {
            if is_salt[idx] {
                let syms: Vec<String> = comp
                    .iter()
                    .take(4)
                    .map(|&a| self.atoms[a].symbol.clone())
                    .collect();
                let label = if comp.len() > 4 {
                    format!("{}...", syms.join(""))
                } else {
                    syms.join("")
                };
                stripped_labels.push(label);
            } else {
                keep_set.extend(comp.iter().copied());
            }
        }

        let mut new_idx_map: Vec<Option<usize>> = vec![None; n];
        let mut new_atoms = Vec::new();
        for (old_idx, mut atom) in self.atoms.drain(..).enumerate() {
            if !keep_set.contains(&old_idx) {
                continue;
            }
            let nid = new_atoms.len();
            atom.id = nid;
            new_idx_map[old_idx] = Some(nid);
            new_atoms.push(atom);
        }
        let mut new_bonds = Vec::new();
        let mut new_adj = vec![Vec::new(); new_atoms.len()];
        for mut b in self.bonds.drain(..) {
            if let (Some(nu), Some(nv)) = (new_idx_map[b.u], new_idx_map[b.v]) {
                let id = new_bonds.len();
                b.id = id;
                b.u = nu;
                b.v = nv;
                new_bonds.push(b);
                new_adj[nu].push((nv, id));
                new_adj[nv].push((nu, id));
            }
        }
        self.atoms = new_atoms;
        self.bonds = new_bonds;
        self.adj = new_adj;
        self.warnings.push(format!(
            "Desalted: stripped {} counterion/salt fragment(s) ({}); preserved parent structure.",
            stripped_labels.len(),
            stripped_labels.join(", ")
        ));
    }

    fn is_recognized_counterion(&self, comp: &[usize]) -> bool {
        let heavy: Vec<usize> = comp
            .iter()
            .copied()
            .filter(|&a| self.atoms[a].atomic_number != 1)
            .collect();
        if heavy.is_empty() {
            // Standalone H+/H2/H- fragment in salt
            return true;
        }
        if heavy.len() == 1 {
            let a = &self.atoms[heavy[0]];
            // Water hydrate (O, OH-, H3O+) or ammonia/ammonium (NH3, NH4+) or halide/HX or common metal counterions
            return match a.atomic_number {
                8 => a.formal_charge >= -1 && a.formal_charge <= 1, // O / [OH-] / [OH3+]
                7 => a.formal_charge == 1,                          // [NH4+]
                9 | 17 | 35 | 53 => true,                           // F, Cl, Br, I and HX
                3 | 11 | 12 | 13 | 19 | 20 | 25 | 26 | 29 | 30 | 37 | 38 | 47 | 55 | 56 | 83 => {
                    true
                } // Li, Na, Mg, Al, K, Ca, Mn, Fe, Cu, Zn, Rb, Sr, Ag, Cs, Ba, Bi
                _ => false,
            };
        }

        let total_charge: i32 = comp
            .iter()
            .map(|&a| i32::from(self.atoms[a].formal_charge))
            .sum();
        let formula = self.component_formula(comp);
        match formula.as_str() {
            // Mineral acids / inorganic oxyanions & fluoroborate/fluorophosphate
            "H2O4S" | "HO4S" | "O4S" => (-2..=0).contains(&total_charge),
            "H3O4P" | "H2O4P" | "HO4P" | "O4P" => (-3..=0).contains(&total_charge),
            "HNO3" | "NO3" => (-1..=0).contains(&total_charge),
            "CH2O3" | "CHO3" | "CO3" => (-2..=0).contains(&total_charge),
            "HClO4" | "ClO4" => (-1..=0).contains(&total_charge),
            "BF4" | "F6P" => total_charge == -1,
            // Common organic pharmaceutical counterions (excluding neutral acetic acid / acetate-Na+)
            "C2HF3O2" | "C2F3O2" => (-1..=0).contains(&total_charge), // TFA
            "CH4O3S" | "CH3O3S" => (-1..=0).contains(&total_charge),  // mesylate
            "C7H8O3S" | "C7H7O3S" => (-1..=0).contains(&total_charge), // tosylate
            "C6H6O3S" | "C6H5O3S" => (-1..=0).contains(&total_charge), // besylate
            "CHO2" => total_charge == -1,                             // formate anion
            "C2H2O4" | "C2HO4" | "C2O4" => (-2..=0).contains(&total_charge), // oxalate
            "C4H4O4" | "C4H3O4" | "C4H2O4" => (-2..=0).contains(&total_charge), // maleate/fumarate
            "C4H6O4" | "C4H5O4" => (-2..=0).contains(&total_charge),  // succinate
            "C4H6O6" | "C4H5O6" | "C4H4O6" => (-2..=0).contains(&total_charge), // tartrate
            "C6H8O7" | "C6H7O7" | "C6H6O7" | "C6H5O7" => (-3..=0).contains(&total_charge), // citrate
            _ => false,
        }
    }

    fn component_formula(&self, comp: &[usize]) -> String {
        let mut counts: BTreeMap<String, u32> = BTreeMap::new();
        for &idx in comp {
            let a = &self.atoms[idx];
            *counts.entry(a.symbol.clone()).or_default() += 1;
            let h = u32::from(a.explicit_h) + u32::from(a.implicit_h);
            if h > 0 && a.atomic_number != 1 {
                *counts.entry("H".to_string()).or_default() += h;
            }
        }
        let mut out = String::new();
        if let Some(c) = counts.remove("C") {
            out.push('C');
            if c > 1 {
                out.push_str(&c.to_string());
            }
            if let Some(h) = counts.remove("H") {
                out.push('H');
                if h > 1 {
                    out.push_str(&h.to_string());
                }
            }
        }
        for (sym, cnt) in counts {
            out.push_str(&sym);
            if cnt > 1 {
                out.push_str(&cnt.to_string());
            }
        }
        out
    }

    /// Strip explicit hydrogen nodes (merging their count onto heavy parent `explicit_h`)
    /// so graphs stay H-suppressed like standard cheminformatics pipelines.
    pub fn suppress_explicit_hydrogens(&mut self) {
        let mut h_parent: Vec<Option<usize>> = vec![None; self.atoms.len()];
        let mut is_remapped_h = vec![false; self.atoms.len()];
        for (idx, atom) in self.atoms.iter().enumerate() {
            if atom.atomic_number == 1 && atom.isotope.is_none() && self.adj[idx].len() == 1 {
                let (nbr, _) = self.adj[idx][0];
                if self.atoms[nbr].atomic_number != 1 {
                    h_parent[idx] = Some(nbr);
                    is_remapped_h[idx] = true;
                }
            }
        }
        if !is_remapped_h.iter().any(|&x| x) {
            return;
        }
        let mut new_idx_map: Vec<Option<usize>> = vec![None; self.atoms.len()];
        let mut new_atoms = Vec::new();
        for (old_idx, mut atom) in self.atoms.drain(..).enumerate() {
            if is_remapped_h[old_idx] {
                continue;
            }
            let new_idx = new_atoms.len();
            atom.id = new_idx;
            new_idx_map[old_idx] = Some(new_idx);
            new_atoms.push(atom);
        }
        for (old_h, parent_opt) in h_parent.into_iter().enumerate() {
            if is_remapped_h[old_h]
                && let Some(old_p) = parent_opt
                && let Some(new_p) = new_idx_map[old_p]
            {
                new_atoms[new_p].explicit_h = new_atoms[new_p].explicit_h.saturating_add(1);
            }
        }
        let mut new_bonds = Vec::new();
        let mut new_adj = vec![Vec::new(); new_atoms.len()];
        for mut b in self.bonds.drain(..) {
            if let (Some(nu), Some(nv)) = (new_idx_map[b.u], new_idx_map[b.v]) {
                let id = new_bonds.len();
                b.id = id;
                b.u = nu;
                b.v = nv;
                new_bonds.push(b);
                new_adj[nu].push((nv, id));
                new_adj[nv].push((nu, id));
            }
        }
        self.atoms = new_atoms;
        self.bonds = new_bonds;
        self.adj = new_adj;
    }

    /// Perceive minimal rings (SSSR-like cycle basis via shortest path on each edge).
    pub fn perceive_rings(&mut self) {
        self.rings.clear();
        for atom in &mut self.atoms {
            atom.in_ring = false;
            atom.ring_sizes.clear();
        }
        for bond in &mut self.bonds {
            bond.in_ring = false;
        }

        let n = self.atoms.len();
        if n < 3 || self.bonds.is_empty() {
            return;
        }

        // Compute connected components & cyclomatic number (M - N + C)
        let mut visited_comp = vec![false; n];
        let mut components = 0usize;
        for start in 0..n {
            if !visited_comp[start] {
                components += 1;
                let mut q = VecDeque::new();
                q.push_back(start);
                visited_comp[start] = true;
                while let Some(u) = q.pop_front() {
                    for &(v, _) in &self.adj[u] {
                        if !visited_comp[v] {
                            visited_comp[v] = true;
                            q.push_back(v);
                        }
                    }
                }
            }
        }
        let target_sssr = self.bonds.len() + components - n;
        if target_sssr == 0 {
            return;
        }

        // Collect candidate shortest cycles for every edge (Horton-style)
        let mut candidates: Vec<(usize, BTreeSet<usize>, Vec<usize>)> = Vec::new();
        for b in &self.bonds {
            if let Some(cycle) = self.shortest_cycle_through_edge(b.u, b.v, b.id) {
                let edge_set: BTreeSet<usize> = cycle
                    .windows(2)
                    .filter_map(|w| self.bond_between(w[0], w[1]).map(|bd| bd.id))
                    .chain(
                        cycle
                            .last()
                            .and_then(|&last| {
                                cycle
                                    .first()
                                    .and_then(|&first| self.bond_between(last, first))
                            })
                            .map(|bd| bd.id),
                    )
                    .collect();
                candidates.push((cycle.len(), edge_set, cycle));
            }
        }

        // Also find cycles passing through vertices of degree >= 3 to catch odd/bridged rings
        for root in 0..n {
            if self.adj[root].len() < 2 {
                continue;
            }
            for i in 0..self.adj[root].len() {
                for j in (i + 1)..self.adj[root].len() {
                    let (n1, b1) = self.adj[root][i];
                    let (n2, b2) = self.adj[root][j];
                    if let Some(path) = self.shortest_path_avoiding_node(n1, n2, root, &[b1, b2]) {
                        let mut cycle = Vec::with_capacity(path.len() + 1);
                        cycle.push(root);
                        cycle.extend(path);
                        let edge_set: BTreeSet<usize> = cycle
                            .windows(2)
                            .filter_map(|w| self.bond_between(w[0], w[1]).map(|bd| bd.id))
                            .chain(
                                cycle
                                    .last()
                                    .and_then(|&last| {
                                        cycle
                                            .first()
                                            .and_then(|&first| self.bond_between(last, first))
                                    })
                                    .map(|bd| bd.id),
                            )
                            .collect();
                        if edge_set.len() == cycle.len() && cycle.len() <= 16 {
                            candidates.push((cycle.len(), edge_set, cycle));
                        }
                    }
                }
            }
        }

        candidates.sort_by_key(|c| c.0);

        // Deduplicate by sorted atom set and select linearly independent basis over GF(2)
        let mut seen_sets: BTreeSet<Vec<usize>> = BTreeSet::new();
        let mut echelon_basis: Vec<(usize, Vec<u64>)> = Vec::new();
        let num_edges = self.bonds.len();

        for (_, edge_set, cycle) in candidates {
            let mut sorted_atoms = cycle.clone();
            sorted_atoms.sort_unstable();
            if !seen_sets.insert(sorted_atoms) {
                continue;
            }
            if insert_gf2_basis(&mut echelon_basis, &edge_set, num_edges) {
                self.rings.push(cycle);
                if self.rings.len() >= target_sssr {
                    break;
                }
            }
        }

        for ring in &self.rings {
            let rlen = ring.len();
            for &aid in ring {
                if let Some(atom) = self.atoms.get_mut(aid) {
                    atom.in_ring = true;
                    atom.ring_sizes.push(rlen);
                }
            }
            for w in ring.windows(2) {
                if let Some(bid) = self.bond_between(w[0], w[1]).map(|b| b.id) {
                    self.bonds[bid].in_ring = true;
                }
            }
            if let (Some(&last), Some(&first)) = (ring.last(), ring.first())
                && let Some(bid) = self.bond_between(last, first).map(|b| b.id)
            {
                self.bonds[bid].in_ring = true;
            }
        }
        for atom in &mut self.atoms {
            atom.ring_sizes.sort_unstable();
        }
    }

    fn shortest_cycle_through_edge(
        &self,
        u: usize,
        v: usize,
        skip_bond: usize,
    ) -> Option<Vec<usize>> {
        let n = self.atoms.len();
        let mut parent = vec![usize::MAX; n];
        let mut dist = vec![usize::MAX; n];
        let mut q = VecDeque::new();
        dist[u] = 0;
        q.push_back(u);
        while let Some(curr) = q.pop_front() {
            if curr == v {
                break;
            }
            if dist[curr] >= 15 {
                continue;
            }
            for &(nbr, bid) in &self.adj[curr] {
                if bid == skip_bond {
                    continue;
                }
                if dist[nbr] == usize::MAX {
                    dist[nbr] = dist[curr] + 1;
                    parent[nbr] = curr;
                    q.push_back(nbr);
                }
            }
        }
        if dist[v] == usize::MAX {
            return None;
        }
        let mut path = Vec::new();
        let mut cur = v;
        while cur != usize::MAX {
            path.push(cur);
            if cur == u {
                break;
            }
            cur = parent[cur];
        }
        path.reverse();
        Some(path)
    }

    fn shortest_path_avoiding_node(
        &self,
        start: usize,
        goal: usize,
        avoid: usize,
        skip_bonds: &[usize],
    ) -> Option<Vec<usize>> {
        let n = self.atoms.len();
        let mut parent = vec![usize::MAX; n];
        let mut dist = vec![usize::MAX; n];
        let mut q = VecDeque::new();
        dist[start] = 0;
        q.push_back(start);
        while let Some(curr) = q.pop_front() {
            if curr == goal {
                break;
            }
            if dist[curr] >= 14 {
                continue;
            }
            for &(nbr, bid) in &self.adj[curr] {
                if nbr == avoid || skip_bonds.contains(&bid) {
                    continue;
                }
                if dist[nbr] == usize::MAX {
                    dist[nbr] = dist[curr] + 1;
                    parent[nbr] = curr;
                    q.push_back(nbr);
                }
            }
        }
        if dist[goal] == usize::MAX {
            return None;
        }
        let mut path = Vec::new();
        let mut cur = goal;
        while cur != usize::MAX {
            path.push(cur);
            if cur == start {
                break;
            }
            cur = parent[cur];
        }
        path.reverse();
        Some(path)
    }

    /// Perceive Hückel 4n+2 aromaticity for rings (both from SMILES aromatic flags and Kekule rings).
    pub fn perceive_aromaticity(&mut self) {
        if self.rings.is_empty() {
            return;
        }

        // First mark any ring whose bonds are already marked Aromatic in SMILES
        let mut aromatic_atoms = vec![false; self.atoms.len()];
        let mut aromatic_bonds = vec![false; self.bonds.len()];
        for b in &self.bonds {
            if b.order == BondOrder::Aromatic {
                aromatic_bonds[b.id] = true;
                aromatic_atoms[b.u] = true;
                aromatic_atoms[b.v] = true;
            }
        }
        for (idx, a) in self.atoms.iter().enumerate() {
            if a.aromatic {
                aromatic_atoms[idx] = true;
            }
        }

        // Check each minimal ring (5..=7) for Hückel 4n+2 pi electrons on planar sp2 atoms
        for ring in &self.rings {
            if ring.len() < 5 || ring.len() > 7 {
                continue;
            }
            if let Some(pi) = self.ring_pi_electrons(ring)
                && pi >= 2
                && (pi - 2) % 4 == 0
            {
                for &aid in ring {
                    aromatic_atoms[aid] = true;
                }
                for w in ring.windows(2) {
                    if let Some(bid) = self.bond_between(w[0], w[1]).map(|b| b.id) {
                        aromatic_bonds[bid] = true;
                    }
                }
                if let (Some(&last), Some(&first)) = (ring.last(), ring.first())
                    && let Some(bid) = self.bond_between(last, first).map(|b| b.id)
                {
                    aromatic_bonds[bid] = true;
                }
            }
        }

        // Also check fused 2-ring aromatic systems (like naphthalene/indole/quinoline/purine/coumarin)
        // where individual rings were already marked or share a double bond
        for i in 0..self.rings.len() {
            for j in (i + 1)..self.rings.len() {
                let r1 = &self.rings[i];
                let r2 = &self.rings[j];
                if r1.len() < 5 || r2.len() < 5 || r1.len() > 7 || r2.len() > 7 {
                    continue;
                }
                let set1: BTreeSet<usize> = r1.iter().copied().collect();
                let set2: BTreeSet<usize> = r2.iter().copied().collect();
                let shared: Vec<usize> = set1.intersection(&set2).copied().collect();
                if shared.len() == 2 {
                    let union: Vec<usize> = set1.union(&set2).copied().collect();
                    if let Some(pi) = self.ring_pi_electrons(&union)
                        && pi >= 2
                        && (pi - 2) % 4 == 0
                    {
                        for &aid in &union {
                            aromatic_atoms[aid] = true;
                        }
                        for &aid in &union {
                            for &(nbr, bid) in &self.adj[aid] {
                                if union.contains(&nbr) {
                                    aromatic_bonds[bid] = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (idx, a) in self.atoms.iter_mut().enumerate() {
            a.aromatic = aromatic_atoms[idx];
        }
        for (idx, b) in self.bonds.iter_mut().enumerate() {
            if aromatic_bonds[idx] {
                b.order = BondOrder::Aromatic;
            }
        }
    }

    fn ring_pi_electrons(&self, ring: &[usize]) -> Option<usize> {
        let ring_set: BTreeSet<usize> = ring.iter().copied().collect();
        let mut total_pi = 0usize;
        for &aid in ring {
            let atom = self.atoms.get(aid)?;
            let mut has_endo_double = false;
            let mut has_endo_aromatic = false;
            let mut has_exo_double_hetero = false;
            let mut total_bonds = 0usize;
            for &(nbr, bid) in &self.adj[aid] {
                total_bonds += 1;
                let bond = self.bonds.get(bid)?;
                if ring_set.contains(&nbr) {
                    if bond.order == BondOrder::Double {
                        has_endo_double = true;
                    } else if bond.order == BondOrder::Aromatic {
                        has_endo_aromatic = true;
                    }
                } else if bond.order == BondOrder::Double
                    && self
                        .atoms
                        .get(nbr)
                        .is_some_and(|n| matches!(n.atomic_number, 7 | 8 | 16))
                {
                    has_exo_double_hetero = true;
                }
            }

            let contrib = match atom.atomic_number {
                6 => {
                    if atom.formal_charge == 1 && total_bonds <= 3 {
                        0
                    } else if atom.formal_charge == -1 && total_bonds <= 3 {
                        2
                    } else if has_exo_double_hetero {
                        // Carbonyl C=O in ring (e.g. coumarin, quinolone, caffeine)
                        0
                    } else if has_endo_double || has_endo_aromatic {
                        1
                    } else {
                        return None;
                    }
                }
                7 => {
                    if atom.formal_charge == 1 && (has_endo_double || has_endo_aromatic) {
                        1
                    } else if has_endo_double {
                        // Pyridine-like N
                        1
                    } else if atom.aromatic
                        && !has_endo_double
                        && atom.explicit_h == 0
                        && total_bonds == 2
                        && !atom.from_bracket
                    {
                        // Aromatic pyridine-like `n` in SMILES
                        1
                    } else if atom.explicit_h > 0
                        || total_bonds == 3
                        || (atom.from_bracket && atom.explicit_h > 0)
                    {
                        // Pyrrole-like [nH] or N-R
                        2
                    } else if has_endo_aromatic {
                        1
                    } else {
                        return None;
                    }
                }
                8 | 16 if total_bonds == 2 && !has_endo_double => 2,
                _ => return None,
            };
            total_pi += contrib;
        }
        Some(total_pi)
    }

    /// Recompute `implicit_h` for non-bracket organic atoms and bracket atoms as appropriate.
    pub fn recompute_implicit_h(&mut self) {
        let mut h_counts = vec![0u8; self.atoms.len()];
        for (idx, atom) in self.atoms.iter().enumerate() {
            if atom.from_bracket {
                // Bracket atoms explicitly encode their H count in SMILES (`explicit_h`).
                h_counts[idx] = 0;
                continue;
            }
            if atom.atomic_number == 1 {
                h_counts[idx] = 0;
                continue;
            }
            let Some(info) = element_by_number(atom.atomic_number) else {
                continue;
            };
            // Compute bond order sum carefully for aromatic vs kekule atoms
            let mut bos_int = 0u8;
            let mut arom_bonds = 0u8;
            let mut exo_double = 0u8;
            for &(_, bid) in &self.adj[idx] {
                if let Some(b) = self.bonds.get(bid) {
                    match b.order {
                        BondOrder::Single => bos_int += 1,
                        BondOrder::Double => {
                            bos_int += 2;
                            exo_double += 1;
                        }
                        BondOrder::Triple => bos_int += 3,
                        BondOrder::Aromatic => {
                            arom_bonds += 1;
                        }
                    }
                }
            }

            let used_valence = if arom_bonds > 0 {
                // In aromatic rings:
                // - Carbon: if it has an exocyclic double bond (e.g. c=O in caffeine/quinolone),
                //   the pi electron goes to the exocyclic double bond -> bos_int + arom_bonds;
                //   otherwise it uses 1 pi electron in the aromatic ring -> bos_int + arom_bonds + 1.
                // - Nitrogen: pyrrole-like substituted N (degree >= 3) uses bos_int + arom_bonds;
                //   pyridine-like `n` uses valence 3 (0 implicit H).
                if atom.atomic_number == 6 {
                    if exo_double > 0 {
                        bos_int + arom_bonds
                    } else {
                        bos_int + arom_bonds + 1
                    }
                } else if atom.atomic_number == 7 {
                    if exo_double > 0 || self.adj[idx].len() >= 3 {
                        bos_int + arom_bonds
                    } else {
                        3
                    }
                } else if atom.atomic_number == 8 || atom.atomic_number == 16 {
                    2 + bos_int
                } else {
                    bos_int + arom_bonds
                }
            } else {
                bos_int
            };

            let target_val = info
                .default_valences
                .iter()
                .copied()
                .find(|&v| v >= used_valence.saturating_add(atom.explicit_h))
                .unwrap_or(used_valence);
            let rem = target_val.saturating_sub(used_valence.saturating_add(atom.explicit_h));
            h_counts[idx] = rem;
        }

        for (idx, h) in h_counts.into_iter().enumerate() {
            self.atoms[idx].implicit_h = h;
        }
    }

    /// Perceive sp / sp2 / sp3 hybridization for every atom.
    pub fn perceive_hybridization(&mut self) {
        let mut hybs = vec![Hybridization::Other; self.atoms.len()];
        for (idx, atom) in self.atoms.iter().enumerate() {
            if atom.atomic_number == 1 {
                hybs[idx] = Hybridization::Other;
                continue;
            }
            if atom.aromatic {
                hybs[idx] = Hybridization::Sp2;
                continue;
            }
            let mut has_triple = false;
            let mut double_count = 0usize;
            for &(_, bid) in &self.adj[idx] {
                if let Some(b) = self.bonds.get(bid) {
                    match b.order {
                        BondOrder::Triple => has_triple = true,
                        BondOrder::Double => double_count += 1,
                        BondOrder::Single | BondOrder::Aromatic => {}
                    }
                }
            }
            hybs[idx] = if has_triple || (double_count >= 2 && matches!(atom.atomic_number, 6 | 7))
            {
                Hybridization::Sp
            } else if double_count == 1 && matches!(atom.atomic_number, 6..=8) {
                Hybridization::Sp2
            } else {
                Hybridization::Sp3
            };
        }
        for (idx, h) in hybs.into_iter().enumerate() {
            self.atoms[idx].hybridization = h;
        }
    }

    /// Hill system molecular formula (e.g. `"C9H8O4"`).
    pub fn formula(&self) -> String {
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for atom in &self.atoms {
            if atom.atomic_number == 1 {
                *counts.entry("H".to_string()).or_insert(0) += 1;
            } else {
                *counts.entry(atom.symbol.clone()).or_insert(0) += 1;
                let h = usize::from(atom.total_h());
                if h > 0 {
                    *counts.entry("H".to_string()).or_insert(0) += h;
                }
            }
        }
        let mut out = String::new();
        if let Some(c) = counts.remove("C") {
            out.push('C');
            if c > 1 {
                out.push_str(&c.to_string());
            }
            if let Some(h) = counts.remove("H") {
                out.push('H');
                if h > 1 {
                    out.push_str(&h.to_string());
                }
            }
        }
        for (sym, cnt) in counts {
            out.push_str(&sym);
            if cnt > 1 {
                out.push_str(&cnt.to_string());
            }
        }
        out
    }

    /// Standard IUPAC average molecular weight (g/mol).
    pub fn molecular_weight(&self) -> f64 {
        let h_weight = element_by_number(1).map_or(1.008, |e| e.atomic_weight);
        self.atoms
            .iter()
            .map(|a| a.atomic_weight() + f64::from(a.total_h()) * h_weight)
            .sum()
    }

    /// Deterministic canonical SMILES representation.
    /// Returns an empty string if ring-closure serialization fails (>99 open rings);
    /// use [`Self::try_to_canonical_smiles`] when fallibility should be reported.
    pub fn to_canonical_smiles(&self) -> String {
        self.try_to_canonical_smiles().unwrap_or_default()
    }

    /// Deterministic canonical SMILES representation, failing if more than 99 rings
    /// are simultaneously open during traversal.
    pub fn try_to_canonical_smiles(&self) -> Result<String, ChemError> {
        if self.atoms.is_empty() {
            return Ok(String::new());
        }
        let n = self.atoms.len();
        // Morgan-like symmetry invariant refinement
        let mut ranks: Vec<u64> = self
            .atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let deg = self.adj[i].len() as u64;
                let z = u64::from(a.atomic_number);
                let h = u64::from(a.total_h());
                let ch = (i64::from(a.formal_charge) + 16) as u64;
                let arom = u64::from(a.aromatic);
                (z << 24) | (deg << 16) | (h << 10) | (ch << 4) | arom
            })
            .collect();

        for _ in 0..6 {
            let mut next_tuples: Vec<(u64, Vec<u64>, usize)> = (0..n)
                .map(|i| {
                    let mut nbr_ranks: Vec<u64> = self.adj[i]
                        .iter()
                        .map(|&(nbr, bid)| {
                            let bo = match self.bonds[bid].order {
                                BondOrder::Single => 1u64,
                                BondOrder::Double => 2,
                                BondOrder::Triple => 3,
                                BondOrder::Aromatic => 4,
                            };
                            (ranks[nbr] << 3) | bo
                        })
                        .collect();
                    nbr_ranks.sort_unstable();
                    (ranks[i], nbr_ranks, i)
                })
                .collect();
            next_tuples.sort();
            let mut new_ranks = vec![0u64; n];
            let mut current_rank = 1u64;
            for idx in 0..next_tuples.len() {
                if idx > 0
                    && (next_tuples[idx].0 != next_tuples[idx - 1].0
                        || next_tuples[idx].1 != next_tuples[idx - 1].1)
                {
                    current_rank += 1;
                }
                new_ranks[next_tuples[idx].2] = current_rank;
            }
            ranks = new_ranks;
        }

        // Break ties deterministically by original index
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| (ranks[i], i));

        let mut visited = vec![false; n];
        let mut component_smiles = Vec::new();
        for &start in &order {
            if visited[start] {
                continue;
            }
            let s = self.serialize_component_smiles(start, &ranks, &mut visited)?;
            component_smiles.push(s);
        }
        Ok(component_smiles.join("."))
    }

    fn serialize_component_smiles(
        &self,
        start: usize,
        ranks: &[u64],
        visited: &mut [bool],
    ) -> Result<String, ChemError> {
        // 1. Iterative DFS spanning tree to identify ring-closure edges (stack-safe for 28k+ atoms)
        let mut tree_visited = vec![false; self.atoms.len()];
        let mut closure_bonds: Vec<usize> = Vec::new();
        let mut closure_set: BTreeSet<usize> = BTreeSet::new();

        struct TreeFrame {
            u: usize,
            parent: usize,
            nbrs: Vec<(usize, usize)>,
            next: usize,
        }

        let mut tree_stack: Vec<TreeFrame> = Vec::new();
        let mut init_nbrs = self.adj[start].clone();
        init_nbrs.sort_by_key(|&(v, _)| (ranks[v], v));
        tree_visited[start] = true;
        tree_stack.push(TreeFrame {
            u: start,
            parent: usize::MAX,
            nbrs: init_nbrs,
            next: 0,
        });

        while let Some(frame) = tree_stack.last_mut() {
            if frame.next >= frame.nbrs.len() {
                tree_stack.pop();
                continue;
            }
            let parent_u = frame.u;
            let (v, bid) = frame.nbrs[frame.next];
            frame.next += 1;
            if v == frame.parent {
                continue;
            }
            if !tree_visited[v] {
                tree_visited[v] = true;
                let mut nbrs = self.adj[v].clone();
                nbrs.sort_by_key(|&(nv, _)| (ranks[nv], nv));
                tree_stack.push(TreeFrame {
                    u: v,
                    parent: parent_u,
                    nbrs,
                    next: 0,
                });
            } else if closure_set.insert(bid) {
                closure_bonds.push(bid);
            }
        }

        // Index closure bonds per endpoint atom
        let mut atom_closure_bonds: Vec<Vec<(usize, usize, BondOrder)>> =
            vec![Vec::new(); self.atoms.len()];
        for &bid in &closure_bonds {
            let b = &self.bonds[bid];
            atom_closure_bonds[b.u].push((bid, b.v, b.order));
            atom_closure_bonds[b.v].push((bid, b.u, b.order));
        }
        for list in &mut atom_closure_bonds {
            list.sort_by_key(|&(_, other, _)| (ranks[other], other));
        }

        // 2. Iterative SMILES emission with reusable ring-closure number pool (1..=99)
        enum EmitOp {
            Visit {
                u: usize,
                p: usize,
                bond_from_parent: Option<BondOrder>,
                parent_aromatic: bool,
            },
            Text(char),
        }

        let mut out = String::new();
        let mut used_nums = [false; 100]; // 1..=99
        let mut active_bond_num: std::collections::HashMap<usize, u16> =
            std::collections::HashMap::new();
        let mut emit_stack = vec![EmitOp::Visit {
            u: start,
            p: usize::MAX,
            bond_from_parent: None,
            parent_aromatic: false,
        }];

        while let Some(op) = emit_stack.pop() {
            match op {
                EmitOp::Text(c) => out.push(c),
                EmitOp::Visit {
                    u,
                    p,
                    bond_from_parent,
                    parent_aromatic,
                } => {
                    if visited[u] {
                        continue;
                    }
                    visited[u] = true;

                    if let Some(border) = bond_from_parent {
                        match border {
                            BondOrder::Double => out.push('='),
                            BondOrder::Triple => out.push('#'),
                            BondOrder::Single => {
                                if parent_aromatic && self.atoms[u].aromatic {
                                    out.push('-');
                                }
                            }
                            BondOrder::Aromatic => {}
                        }
                    }

                    out.push_str(&self.format_smiles_atom(u));

                    // Close already-opened rings first so numbers return to pool immediately
                    let mut to_open: Vec<(usize, BondOrder)> = Vec::new();
                    for &(bid, _, border) in &atom_closure_bonds[u] {
                        if let Some(num) = active_bond_num.remove(&bid) {
                            push_ring_tag(&mut out, num);
                            if usize::from(num) < used_nums.len() {
                                used_nums[usize::from(num)] = false;
                            }
                        } else {
                            to_open.push((bid, border));
                        }
                    }

                    // Open new rings using lowest available ring number 1..=99
                    for (bid, border) in to_open {
                        let Some(num) = (1..=99u16).find(|&n| !used_nums[usize::from(n)]) else {
                            return Err(ChemError::TooManyOpenRings);
                        };
                        used_nums[usize::from(num)] = true;
                        active_bond_num.insert(bid, num);
                        if border == BondOrder::Double {
                            out.push('=');
                        } else if border == BondOrder::Triple {
                            out.push('#');
                        }
                        push_ring_tag(&mut out, num);
                    }

                    let mut children: Vec<(usize, usize)> = self.adj[u]
                        .iter()
                        .copied()
                        .filter(|&(v, bid)| v != p && !closure_set.contains(&bid) && !visited[v])
                        .collect();
                    children.sort_by_key(|&(v, _)| (ranks[v], v));

                    let n_children = children.len();
                    for (idx, (v, bid)) in children.into_iter().enumerate().rev() {
                        let is_last = idx + 1 == n_children;
                        if !is_last {
                            emit_stack.push(EmitOp::Text(')'));
                        }
                        emit_stack.push(EmitOp::Visit {
                            u: v,
                            p: u,
                            bond_from_parent: Some(self.bonds[bid].order),
                            parent_aromatic: self.atoms[u].aromatic,
                        });
                        if !is_last {
                            emit_stack.push(EmitOp::Text('('));
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    fn format_smiles_atom(&self, u: usize) -> String {
        let a = &self.atoms[u];
        let organic = matches!(
            a.symbol.as_str(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
        let needs_brackets = !organic
            || a.formal_charge != 0
            || a.isotope.is_some()
            || (a.aromatic && a.symbol == "N" && a.explicit_h > 0)
            || (a.aromatic && !matches!(a.symbol.as_str(), "B" | "C" | "N" | "O" | "P" | "S"));

        if !needs_brackets {
            if a.aromatic {
                return a.symbol.to_lowercase();
            }
            return a.symbol.clone();
        }
        let mut s = String::from("[");
        if let Some(iso) = a.isotope {
            s.push_str(&iso.to_string());
        }
        if a.aromatic {
            s.push_str(&a.symbol.to_lowercase());
        } else {
            s.push_str(&a.symbol);
        }
        let h = a.total_h();
        if h == 1 {
            s.push('H');
        } else if h > 1 {
            s.push('H');
            s.push_str(&h.to_string());
        }
        if a.formal_charge > 0 {
            s.push('+');
            if a.formal_charge > 1 {
                s.push_str(&a.formal_charge.to_string());
            }
        } else if a.formal_charge < 0 {
            s.push('-');
            let mag = i16::from(a.formal_charge).abs();
            if mag > 1 {
                s.push_str(&mag.to_string());
            }
        }
        s.push(']');
        s
    }
}

fn push_ring_tag(out: &mut String, num: u16) {
    if num < 10 {
        out.push_str(&num.to_string());
    } else {
        out.push('%');
        out.push_str(&format!("{num:02}"));
    }
}

fn insert_gf2_basis(
    echelon: &mut Vec<(usize, Vec<u64>)>,
    candidate: &BTreeSet<usize>,
    num_edges: usize,
) -> bool {
    let mut cand = set_to_words(candidate, num_edges);
    for (pivot, row) in echelon.iter() {
        if get_bit(&cand, *pivot) {
            xor_words(&mut cand, row);
        }
    }
    if let Some(new_pivot) = first_set_bit(&cand) {
        echelon.push((new_pivot, cand));
        echelon.sort_by_key(|(p, _)| *p);
        true
    } else {
        false
    }
}

fn set_to_words(s: &BTreeSet<usize>, num_edges: usize) -> Vec<u64> {
    let words = num_edges.div_ceil(64);
    let mut v = vec![0u64; words];
    for &bit in s {
        let w = bit / 64;
        let b = bit % 64;
        if w < v.len() {
            v[w] |= 1u64 << b;
        }
    }
    v
}

fn first_set_bit(words: &[u64]) -> Option<usize> {
    for (wi, &w) in words.iter().enumerate() {
        if w != 0 {
            return Some(wi * 64 + w.trailing_zeros() as usize);
        }
    }
    None
}

fn get_bit(words: &[u64], bit: usize) -> bool {
    let w = bit / 64;
    let b = bit % 64;
    words.get(w).is_some_and(|&val| (val & (1u64 << b)) != 0)
}

fn xor_words(dst: &mut [u64], src: &[u64]) {
    for (d, s) in dst.iter_mut().zip(src.iter()) {
        *d ^= *s;
    }
}
