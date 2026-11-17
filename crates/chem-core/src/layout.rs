use std::collections::{BTreeSet, VecDeque};
use std::f64::consts::PI;

use crate::mol::{BondOrder, Molecule};

const BOND_LEN: f64 = 1.5;

/// Generate 2D coordinates for all atoms if missing (or regenerate cleanly).
pub fn ensure_2d_coords(mol: &mut Molecule) {
    if mol.atoms.is_empty() {
        return;
    }
    if mol.atoms.iter().all(|a| a.coords_2d.is_some()) {
        return;
    }
    generate_2d_coords(mol);
}

/// Generate clean 2D coordinates for all atoms in `mol`.
pub fn generate_2d_coords(mol: &mut Molecule) {
    let n = mol.atoms.len();
    if n == 0 {
        return;
    }
    let mut coords = vec![[0.0f64; 2]; n];
    let mut placed = vec![false; n];

    // Group rings into fused/spiro ring systems
    let ring_systems = build_ring_systems(mol);

    // Find connected components of the molecule
    let components = connected_components(mol);
    let mut x_offset = 0.0;

    for comp in components {
        let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
        // Place ring systems belonging to this component first, or start at first atom
        let mut comp_systems: Vec<&Vec<usize>> = ring_systems
            .iter()
            .filter(|rs| {
                rs.iter()
                    .any(|&ri| mol.rings[ri].iter().any(|a| comp_set.contains(a)))
            })
            .collect();
        comp_systems.sort_by_key(|rs| std::cmp::Reverse(rs.len()));

        if let Some(first_sys) = comp_systems.first() {
            place_ring_system(mol, first_sys, &mut coords, &mut placed);
        } else if let Some(&root) = comp.first() {
            coords[root] = [0.0, 0.0];
            placed[root] = true;
        }

        // BFS outward to place remaining atoms & ring systems in this component
        let mut q: VecDeque<usize> = comp.iter().copied().filter(|&a| placed[a]).collect();
        while let Some(u) = q.pop_front() {
            // Collect unplaced neighbors
            let unplaced_nbrs: Vec<usize> = mol.neighbors(u).filter(|&v| !placed[v]).collect();
            if unplaced_nbrs.is_empty() {
                continue;
            }

            // Compute average incoming angle from already-placed neighbors
            let placed_nbrs: Vec<usize> = mol.neighbors(u).filter(|&v| placed[v]).collect();

            let base_angle = if placed_nbrs.is_empty() {
                0.0
            } else if placed_nbrs.len() == 1 {
                let p = placed_nbrs[0];
                let dx = coords[u][0] - coords[p][0];
                let dy = coords[u][1] - coords[p][1];
                let incoming = dy.atan2(dx);
                // Linear for sp triple bonds, 120 deg otherwise
                let is_sp = mol
                    .bond_between(u, p)
                    .is_some_and(|b| b.order == BondOrder::Triple)
                    || unplaced_nbrs.iter().any(|&v| {
                        mol.bond_between(u, v)
                            .is_some_and(|b| b.order == BondOrder::Triple)
                    });
                if is_sp && unplaced_nbrs.len() == 1 {
                    incoming
                } else {
                    incoming + PI / 3.0
                }
            } else {
                // Bisect the largest angular gap around u
                best_open_angle(&coords, u, &placed_nbrs)
            };

            let step = if unplaced_nbrs.len() == 1 {
                0.0
            } else {
                2.0 * PI / (unplaced_nbrs.len() + placed_nbrs.len().max(1)) as f64
            };

            for (idx, &v) in unplaced_nbrs.iter().enumerate() {
                // Check if v belongs to an unplaced ring system
                if let Some(sys) = comp_systems
                    .iter()
                    .find(|rs| rs.iter().any(|&ri| mol.rings[ri].contains(&v)))
                {
                    let ang = base_angle + idx as f64 * step;
                    coords[v] = [
                        coords[u][0] + BOND_LEN * ang.cos(),
                        coords[u][1] + BOND_LEN * ang.sin(),
                    ];
                    placed[v] = true;
                    place_ring_system_from_anchor(mol, sys, v, u, &mut coords, &mut placed);
                    for &ri in *sys {
                        for &ra in &mol.rings[ri] {
                            q.push_back(ra);
                        }
                    }
                } else {
                    let ang = if placed_nbrs.len() == 1 && unplaced_nbrs.len() == 1 {
                        let p = placed_nbrs[0];
                        let incoming =
                            (coords[u][1] - coords[p][1]).atan2(coords[u][0] - coords[p][0]);
                        let is_sp = mol
                            .bond_between(u, p)
                            .is_some_and(|b| b.order == BondOrder::Triple)
                            || mol
                                .bond_between(u, v)
                                .is_some_and(|b| b.order == BondOrder::Triple);
                        if is_sp {
                            incoming
                        } else {
                            // Choose +60 deg or -60 deg maximizing clearance from already-placed atoms
                            let a_plus = incoming + PI / 3.0;
                            let a_minus = incoming - PI / 3.0;
                            let p_plus = [
                                coords[u][0] + BOND_LEN * a_plus.cos(),
                                coords[u][1] + BOND_LEN * a_plus.sin(),
                            ];
                            let p_minus = [
                                coords[u][0] + BOND_LEN * a_minus.cos(),
                                coords[u][1] + BOND_LEN * a_minus.sin(),
                            ];
                            if min_dist_to_placed(p_plus, &coords, &placed, u)
                                >= min_dist_to_placed(p_minus, &coords, &placed, u)
                            {
                                a_plus
                            } else {
                                a_minus
                            }
                        }
                    } else if unplaced_nbrs.len() == 2 && placed_nbrs.len() == 1 {
                        let p = placed_nbrs[0];
                        let incoming =
                            (coords[u][1] - coords[p][1]).atan2(coords[u][0] - coords[p][0]);
                        // Test candidate tilts (-25, 0, +25 deg) to avoid colliding with ortho neighbors
                        let mut best_offset = 0.0;
                        let mut best_clearance = -1.0;
                        for &off in &[0.0, 0.35, -0.35, 0.65, -0.65] {
                            let a1 = incoming + off + PI / 3.0;
                            let a2 = incoming + off - PI / 3.0;
                            let c1 = [
                                coords[u][0] + BOND_LEN * a1.cos(),
                                coords[u][1] + BOND_LEN * a1.sin(),
                            ];
                            let c2 = [
                                coords[u][0] + BOND_LEN * a2.cos(),
                                coords[u][1] + BOND_LEN * a2.sin(),
                            ];
                            let cl = min_dist_to_placed(c1, &coords, &placed, u)
                                .min(min_dist_to_placed(c2, &coords, &placed, u));
                            if cl > best_clearance {
                                best_clearance = cl;
                                best_offset = off;
                            }
                        }
                        if idx == 0 {
                            incoming + best_offset + PI / 3.0
                        } else {
                            incoming + best_offset - PI / 3.0
                        }
                    } else {
                        let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
                        let mult = idx.div_ceil(2) as f64;
                        base_angle + sign * mult * (2.0 * PI / 3.0).max(step)
                    };
                    coords[v] = [
                        coords[u][0] + BOND_LEN * ang.cos(),
                        coords[u][1] + BOND_LEN * ang.sin(),
                    ];
                    placed[v] = true;
                    q.push_back(v);
                }
            }
        }

        // Gentle spring-repulsion relaxation for non-bonded overlap cleanup
        relax_overlaps(mol, &comp, &mut coords, 80);

        // Shift component horizontally so multiple fragments don't collide
        let min_x = comp
            .iter()
            .map(|&a| coords[a][0])
            .fold(f64::INFINITY, f64::min);
        let max_x = comp
            .iter()
            .map(|&a| coords[a][0])
            .fold(f64::NEG_INFINITY, f64::max);
        let shift = x_offset - min_x;
        for &a in &comp {
            coords[a][0] += shift;
        }
        x_offset += (max_x - min_x) + 3.0;
    }

    // Center around origin
    let cx = coords.iter().map(|c| c[0]).sum::<f64>() / n as f64;
    let cy = coords.iter().map(|c| c[1]).sum::<f64>() / n as f64;
    for (idx, atom) in mol.atoms.iter_mut().enumerate() {
        atom.coords_2d = Some([coords[idx][0] - cx, coords[idx][1] - cy]);
    }
}

fn connected_components(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut seen = vec![false; n];
    let mut comps = Vec::new();
    for i in 0..n {
        if !seen[i] {
            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            seen[i] = true;
            q.push_back(i);
            while let Some(u) = q.pop_front() {
                comp.push(u);
                for v in mol.neighbors(u) {
                    if !seen[v] {
                        seen[v] = true;
                        q.push_back(v);
                    }
                }
            }
            comps.push(comp);
        }
    }
    comps
}

fn build_ring_systems(mol: &Molecule) -> Vec<Vec<usize>> {
    let nr = mol.rings.len();
    if nr == 0 {
        return Vec::new();
    }
    let ring_sets: Vec<BTreeSet<usize>> = mol
        .rings
        .iter()
        .map(|r| r.iter().copied().collect())
        .collect();
    let mut seen = vec![false; nr];
    let mut systems = Vec::new();
    for i in 0..nr {
        if !seen[i] {
            let mut sys = Vec::new();
            let mut q = VecDeque::new();
            seen[i] = true;
            q.push_back(i);
            while let Some(ri) = q.pop_front() {
                sys.push(ri);
                for (rj, set_j) in ring_sets.iter().enumerate() {
                    if !seen[rj] && ring_sets[ri].intersection(set_j).next().is_some() {
                        seen[rj] = true;
                        q.push_back(rj);
                    }
                }
            }
            systems.push(sys);
        }
    }
    systems
}

fn place_ring_system(mol: &Molecule, sys: &[usize], coords: &mut [[f64; 2]], placed: &mut [bool]) {
    let Some(&first_ri) = sys.first() else {
        return;
    };
    place_single_polygon(&mol.rings[first_ri], [0.0, 0.0], 0.0, coords, placed);

    let mut remaining: Vec<usize> = sys[1..].to_vec();
    for _ in 0..sys.len() {
        let mut next_remaining = Vec::new();
        let mut progressed = false;
        for ri in remaining {
            let ring = &mol.rings[ri];
            let placed_in_ring: Vec<usize> = ring.iter().copied().filter(|&a| placed[a]).collect();
            if placed_in_ring.len() >= 2 {
                // Fuse along first consecutive pair of placed atoms in ring
                place_fused_ring(ring, coords, placed);
                progressed = true;
            } else if placed_in_ring.len() == 1 {
                let anchor = placed_in_ring[0];
                place_spiro_ring(mol, ring, anchor, coords, placed);
                progressed = true;
            } else {
                next_remaining.push(ri);
            }
        }
        remaining = next_remaining;
        if remaining.is_empty() || !progressed {
            break;
        }
    }
    // Fallback for any still unplaced rings in complex cages
    for ri in remaining {
        for &a in &mol.rings[ri] {
            if !placed[a] {
                coords[a] = [0.5 * (a as f64), 0.5 * ((a % 3) as f64)];
                placed[a] = true;
            }
        }
    }
}

fn place_ring_system_from_anchor(
    mol: &Molecule,
    sys: &[usize],
    anchor: usize,
    parent: usize,
    coords: &mut [[f64; 2]],
    placed: &mut [bool],
) {
    if let Some(&ri) = sys.iter().find(|&&ri| mol.rings[ri].contains(&anchor)) {
        let ring = &mol.rings[ri];
        let k = ring.len() as f64;
        let radius = BOND_LEN / (2.0 * (PI / k).sin());
        let dir_x = coords[anchor][0] - coords[parent][0];
        let dir_y = coords[anchor][1] - coords[parent][1];
        let norm = (dir_x * dir_x + dir_y * dir_y).sqrt().max(1e-6);
        let cx = coords[anchor][0] + radius * (dir_x / norm);
        let cy = coords[anchor][1] + radius * (dir_y / norm);
        let anchor_pos = ring.iter().position(|&a| a == anchor).unwrap_or(0);
        let base_ang = (coords[anchor][1] - cy).atan2(coords[anchor][0] - cx);
        let step = 2.0 * PI / k;
        for (offset, &aid) in ring.iter().enumerate() {
            if !placed[aid] {
                let delta = (offset as isize - anchor_pos as isize) as f64;
                let ang = base_ang + delta * step;
                coords[aid] = [cx + radius * ang.cos(), cy + radius * ang.sin()];
                placed[aid] = true;
            }
        }
    }
    place_ring_system_remaining(mol, sys, coords, placed);
}

fn place_ring_system_remaining(
    mol: &Molecule,
    sys: &[usize],
    coords: &mut [[f64; 2]],
    placed: &mut [bool],
) {
    for _ in 0..sys.len() {
        let mut progressed = false;
        for &ri in sys {
            let ring = &mol.rings[ri];
            if ring.iter().all(|&a| placed[a]) {
                continue;
            }
            let placed_cnt = ring.iter().filter(|&&a| placed[a]).count();
            if placed_cnt >= 2 {
                place_fused_ring(ring, coords, placed);
                progressed = true;
            } else if placed_cnt == 1
                && let Some(&anchor) = ring.iter().find(|&&a| placed[a])
            {
                place_spiro_ring(mol, ring, anchor, coords, placed);
                progressed = true;
            }
        }
        if !progressed {
            break;
        }
    }
}

fn place_single_polygon(
    ring: &[usize],
    center: [f64; 2],
    start_angle: f64,
    coords: &mut [[f64; 2]],
    placed: &mut [bool],
) {
    let k = ring.len().max(3) as f64;
    let radius = BOND_LEN / (2.0 * (PI / k).sin());
    for (idx, &aid) in ring.iter().enumerate() {
        let ang = start_angle + (idx as f64) * 2.0 * PI / k;
        coords[aid] = [
            center[0] + radius * ang.cos(),
            center[1] + radius * ang.sin(),
        ];
        placed[aid] = true;
    }
}

fn place_fused_ring(ring: &[usize], coords: &mut [[f64; 2]], placed: &mut [bool]) {
    let m = ring.len();
    if m < 3 {
        return;
    }
    let placed_cnt = ring.iter().filter(|&&a| placed[a]).count();
    if placed_cnt >= 3 {
        place_bridged_ring(ring, coords, placed);
        return;
    }
    // Find an adjacent placed pair (ring[i], ring[(i+1)%m])
    let mut pair_idx = None;
    for i in 0..m {
        if placed[ring[i]] && placed[ring[(i + 1) % m]] {
            pair_idx = Some(i);
            break;
        }
    }
    let Some(i0) = pair_idx else {
        return;
    };
    let a = ring[i0];
    let b = ring[(i0 + 1) % m];
    let mid = [
        0.5 * (coords[a][0] + coords[b][0]),
        0.5 * (coords[a][1] + coords[b][1]),
    ];
    let dx = coords[b][0] - coords[a][0];
    let dy = coords[b][1] - coords[a][1];
    let edge_len = (dx * dx + dy * dy).sqrt().max(1e-6);
    let nx = -dy / edge_len;
    let ny = dx / edge_len;

    let k = m as f64;
    let apothem = (edge_len * 0.5) / (PI / k).tan();
    let radius = (edge_len * 0.5) / (PI / k).sin();

    // Pick normal direction pointing away from already-placed atoms centroid
    let placed_pts: Vec<[f64; 2]> = coords
        .iter()
        .enumerate()
        .filter(|(idx, _)| placed[*idx])
        .map(|(_, c)| *c)
        .collect();
    let c_mean = if placed_pts.is_empty() {
        [0.0, 0.0]
    } else {
        [
            placed_pts.iter().map(|p| p[0]).sum::<f64>() / placed_pts.len() as f64,
            placed_pts.iter().map(|p| p[1]).sum::<f64>() / placed_pts.len() as f64,
        ]
    };
    let c1 = [mid[0] + apothem * nx, mid[1] + apothem * ny];
    let c2 = [mid[0] - apothem * nx, mid[1] - apothem * ny];
    let d1 = (c1[0] - c_mean[0]).powi(2) + (c1[1] - c_mean[1]).powi(2);
    let d2 = (c2[0] - c_mean[0]).powi(2) + (c2[1] - c_mean[1]).powi(2);
    let center = if d1 >= d2 { c1 } else { c2 };

    let ang_a = (coords[a][1] - center[1]).atan2(coords[a][0] - center[0]);
    let ang_b = (coords[b][1] - center[1]).atan2(coords[b][0] - center[0]);
    let mut diff = ang_b - ang_a;
    while diff > PI {
        diff -= 2.0 * PI;
    }
    while diff < -PI {
        diff += 2.0 * PI;
    }
    let step = if diff >= 0.0 {
        2.0 * PI / k
    } else {
        -2.0 * PI / k
    };

    for offset in 0..m {
        let aid = ring[(i0 + offset) % m];
        if !placed[aid] {
            let ang = ang_a + (offset as f64) * step;
            coords[aid] = [
                center[0] + radius * ang.cos(),
                center[1] + radius * ang.sin(),
            ];
            placed[aid] = true;
        }
    }
}

fn place_bridged_ring(ring: &[usize], coords: &mut [[f64; 2]], placed: &mut [bool]) {
    let m = ring.len();
    let placed_pts: Vec<[f64; 2]> = ring
        .iter()
        .filter(|&&a| placed[a])
        .map(|&a| coords[a])
        .collect();
    if placed_pts.is_empty() {
        return;
    }
    let cx = placed_pts.iter().map(|p| p[0]).sum::<f64>() / placed_pts.len() as f64;
    let cy = placed_pts.iter().map(|p| p[1]).sum::<f64>() / placed_pts.len() as f64;

    // Walk around ring and place any contiguous unplaced run between its placed flanking atoms
    let Some(first_placed) = ring.iter().position(|&a| placed[a]) else {
        return;
    };
    let mut i = 0usize;
    while i < m {
        let idx = (first_placed + i) % m;
        if placed[ring[idx]] {
            i += 1;
            continue;
        }
        let start_anchor = ring[(first_placed + i + m - 1) % m];
        let mut run = Vec::new();
        while i < m && !placed[ring[(first_placed + i) % m]] {
            run.push(ring[(first_placed + i) % m]);
            i += 1;
        }
        let end_anchor = ring[(first_placed + i) % m];
        let p0 = coords[start_anchor];
        let p1 = coords[end_anchor];
        let dx = p1[0] - p0[0];
        let dy = p1[1] - p0[1];
        let chord = (dx * dx + dy * dy).sqrt().max(1e-6);
        let nx = -dy / chord;
        let ny = dx / chord;
        let mid = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])];
        // Bow slightly away from the centroid of placed ring atoms (or perpendicular if collinear)
        let dot = (mid[0] - cx) * nx + (mid[1] - cy) * ny;
        let sign = if dot >= 0.0 { 1.0 } else { -1.0 };
        let k = run.len();
        for (j, &aid) in run.iter().enumerate() {
            let t = (j + 1) as f64 / (k + 1) as f64;
            let arch = (PI * t).sin() * BOND_LEN * 0.45 * sign;
            coords[aid] = [p0[0] + dx * t + nx * arch, p0[1] + dy * t + ny * arch];
            placed[aid] = true;
        }
    }
}

fn place_spiro_ring(
    mol: &Molecule,
    ring: &[usize],
    anchor: usize,
    coords: &mut [[f64; 2]],
    placed: &mut [bool],
) {
    let placed_nbrs: Vec<usize> = mol.neighbors(anchor).filter(|&v| placed[v]).collect();
    let dir = if placed_nbrs.is_empty() {
        0.0
    } else {
        best_open_angle(coords, anchor, &placed_nbrs)
    };
    let k = ring.len().max(3) as f64;
    let radius = BOND_LEN / (2.0 * (PI / k).sin());
    let center = [
        coords[anchor][0] + radius * dir.cos(),
        coords[anchor][1] + radius * dir.sin(),
    ];
    let anchor_pos = ring.iter().position(|&a| a == anchor).unwrap_or(0);
    let base_ang = (coords[anchor][1] - center[1]).atan2(coords[anchor][0] - center[0]);
    let step = 2.0 * PI / k;
    for (offset, &aid) in ring.iter().enumerate() {
        if !placed[aid] {
            let delta = (offset as isize - anchor_pos as isize) as f64;
            let ang = base_ang + delta * step;
            coords[aid] = [
                center[0] + radius * ang.cos(),
                center[1] + radius * ang.sin(),
            ];
            placed[aid] = true;
        }
    }
}

fn best_open_angle(coords: &[[f64; 2]], center: usize, nbrs: &[usize]) -> f64 {
    let mut angles: Vec<f64> = nbrs
        .iter()
        .map(|&v| (coords[v][1] - coords[center][1]).atan2(coords[v][0] - coords[center][0]))
        .collect();
    angles.sort_by(|a, b| a.total_cmp(b));
    if angles.is_empty() {
        return 0.0;
    }
    if angles.len() == 1 {
        return angles[0] + PI;
    }
    let mut best_gap = -1.0;
    let mut best_mid = 0.0;
    for i in 0..angles.len() {
        let a1 = angles[i];
        let a2 = if i + 1 < angles.len() {
            angles[i + 1]
        } else {
            angles[0] + 2.0 * PI
        };
        let gap = a2 - a1;
        if gap > best_gap {
            best_gap = gap;
            best_mid = a1 + gap * 0.5;
        }
    }
    best_mid
}

fn min_dist_to_placed(pos: [f64; 2], coords: &[[f64; 2]], placed: &[bool], ignore: usize) -> f64 {
    let mut min_d = f64::INFINITY;
    for (i, &is_p) in placed.iter().enumerate() {
        if !is_p || i == ignore {
            continue;
        }
        let dx = pos[0] - coords[i][0];
        let dy = pos[1] - coords[i][1];
        let d = (dx * dx + dy * dy).sqrt();
        if d < min_d {
            min_d = d;
        }
    }
    min_d
}

fn relax_overlaps(mol: &Molecule, comp: &[usize], coords: &mut [[f64; 2]], iters: usize) {
    if comp.len() < 3 {
        return;
    }
    let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
    for _ in 0..iters {
        let mut forces = vec![[0.0f64; 2]; coords.len()];
        // Non-bonded repulsion for atoms closer than 1.45 units
        for (i, &u) in comp.iter().enumerate() {
            for &v in &comp[(i + 1)..] {
                if mol.bond_between(u, v).is_some() {
                    continue;
                }
                let dx = coords[u][0] - coords[v][0];
                let dy = coords[u][1] - coords[v][1];
                let d2 = dx * dx + dy * dy;
                if d2 <= 1e-6 {
                    let ang = ((u * 7 + v * 13) as f64) * 0.7;
                    forces[u][0] += 0.15 * ang.cos();
                    forces[u][1] += 0.15 * ang.sin();
                    forces[v][0] -= 0.15 * ang.cos();
                    forces[v][1] -= 0.15 * ang.sin();
                } else if d2 < 2.25 {
                    let d = d2.sqrt();
                    let push = 0.14 * (1.50 - d) / d;
                    forces[u][0] += dx * push;
                    forces[u][1] += dy * push;
                    forces[v][0] -= dx * push;
                    forces[v][1] -= dy * push;
                }
            }
        }
        // Bond length spring restoration
        for b in &mol.bonds {
            if !comp_set.contains(&b.u) {
                continue;
            }
            let dx = coords[b.v][0] - coords[b.u][0];
            let dy = coords[b.v][1] - coords[b.u][1];
            let d = (dx * dx + dy * dy).sqrt().max(1e-6);
            let stiffness = if b.in_ring { 0.25 } else { 0.18 };
            let pull = stiffness * (d - BOND_LEN) / d;
            forces[b.u][0] += dx * pull;
            forces[b.u][1] += dy * pull;
            forces[b.v][0] -= dx * pull;
            forces[b.v][1] -= dy * pull;
        }
        // Keep simple planar aromatic rings rigid; allow bridged rings & acyclic atoms to relax
        for &u in comp {
            let is_simple_aromatic_ring =
                mol.atoms[u].in_ring && mol.atoms[u].aromatic && mol.atoms[u].ring_sizes.len() <= 2;
            if !is_simple_aromatic_ring {
                coords[u][0] += forces[u][0].clamp(-0.22, 0.22);
                coords[u][1] += forces[u][1].clamp(-0.22, 0.22);
            }
        }
    }
}
