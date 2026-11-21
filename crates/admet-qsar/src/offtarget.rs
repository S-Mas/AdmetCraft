#![forbid(unsafe_code)]

//! Pure-Rust **SEA (Similarity Ensemble Approach)** off-target safety pharmacology fishing
//! (Keiser et al., *Nature Biotechnology* 2007) over 1024-bit ECFP4 fingerprints (`u64::count_ones`).
//!
//! Ships with a built-in curated **Safety-20 Bowes Reference Panel** covering 20 critical human
//! anti-targets (hERG, 5-HT2B, 5-HT2A, DRD2, ADRB2, ADRA1A, CHRM1, HRH1, OPRM1, COX1, COX2, AChE,
//! PXR, PPARγ, ERα, AR, NR3C1, DAT, SERT, CYP3A4) so off-target prediction works 100 % offline out
//! of the box, and supports loading/building extended ChEMBL indices via `build_sea_index_from_tsv`.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::OnceLock;

use admet_calc::{OffTargetHit, OffTargetReport, OffTargetRisk};
use chem_core::{Molecule, parse_smiles};
use serde::{Deserialize, Serialize};

use crate::data_dir::default_sea_index_path;
use crate::fingerprint::{Ecfp4BitVec, compute_ecfp4};
use crate::pack::QsarError;

/// Tanimoto similarity threshold for SEA raw-score summation (Keiser et al. ECFP4 standard).
pub const SEA_TC_CUTOFF: f64 = 0.28;
const PI_OVER_SQRT6: f64 = 1.282_549_830_161_864;

/// Known active reference ligand inside a SEA target set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaReferenceLigand {
    pub name: String,
    pub smiles: String,
    pub fp: Ecfp4BitVec,
}

/// One human anti-target set with background-calibrated SEA distribution statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaTargetSet {
    pub target_id: String,
    pub target_name: String,
    pub target_class: String,
    pub organ_liability: String,
    pub mu_bg: f64,
    pub sigma_bg: f64,
    pub ligands: Vec<SeaReferenceLigand>,
}

/// Complete serializable SEA Off-Target Fishing index.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaOffTargetIndex {
    pub format: String,
    pub version: String,
    pub tc_cutoff: f64,
    pub targets: Vec<SeaTargetSet>,
}

impl SeaOffTargetIndex {
    pub fn load_from_path(path: &Path) -> Result<Self, QsarError> {
        let raw = std::fs::read_to_string(path)?;
        let idx: Self = serde_json::from_str(&raw)?;
        Ok(idx)
    }

    pub fn save_to_path(&self, path: &Path) -> Result<(), QsarError> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let raw = serde_json::to_string(self)?;
        std::fs::write(path, raw)?;
        Ok(())
    }

    /// Predicts human off-target liabilities for a query ECFP4 fingerprint via SEA Z-score & Gumbel E-value.
    #[must_use]
    pub fn predict_fp(&self, query_fp: &Ecfp4BitVec) -> OffTargetReport {
        let mut hits = Vec::new();
        let mut high_cnt = 0usize;
        let mut mod_cnt = 0usize;
        let n_targets = self.targets.len().max(1) as f64;

        for tgt in &self.targets {
            if tgt.ligands.is_empty() {
                continue;
            }
            let mut raw_score = 0.0_f64;
            let mut max_tc = 0.0_f64;
            let mut nearest = String::new();
            for lig in &tgt.ligands {
                let tc = query_fp.tanimoto(&lig.fp);
                if tc >= self.tc_cutoff {
                    raw_score += tc;
                }
                if tc > max_tc {
                    max_tc = tc;
                    nearest = lig.name.clone();
                }
            }
            let sigma = tgt.sigma_bg.max(0.12);
            let z = (raw_score - tgt.mu_bg) / sigma;
            let p_val = gumbel_tail_p(z);
            let e_val = (p_val * n_targets).min(n_targets);

            let risk = if max_tc >= 0.58 || (z >= 4.0 && max_tc >= 0.40) || e_val <= 1e-4 {
                OffTargetRisk::High
            } else if max_tc >= 0.36 || (z >= 2.2 && max_tc >= 0.28) || e_val <= 0.05 {
                OffTargetRisk::Moderate
            } else {
                OffTargetRisk::Low
            };

            if risk == OffTargetRisk::High {
                high_cnt += 1;
            } else if risk == OffTargetRisk::Moderate {
                mod_cnt += 1;
            }

            if risk != OffTargetRisk::Low || max_tc >= 0.30 {
                hits.push(OffTargetHit {
                    target_id: tgt.target_id.clone(),
                    target_name: tgt.target_name.clone(),
                    target_class: tgt.target_class.clone(),
                    organ_liability: tgt.organ_liability.clone(),
                    max_tanimoto: (max_tc * 1000.0).round() / 1000.0,
                    nearest_ligand: nearest,
                    raw_score: (raw_score * 1000.0).round() / 1000.0,
                    z_score: (z * 100.0).round() / 100.0,
                    e_value: e_val,
                    risk,
                });
            }
        }

        hits.sort_by(|a, b| {
            risk_rank(b.risk)
                .cmp(&risk_rank(a.risk))
                .then_with(|| b.max_tanimoto.total_cmp(&a.max_tanimoto))
                .then_with(|| b.z_score.total_cmp(&a.z_score))
        });
        hits.truncate(12);

        OffTargetReport {
            panel_version: self.version.clone(),
            targets_evaluated: self.targets.len(),
            high_risk_count: high_cnt,
            moderate_risk_count: mod_cnt,
            hits,
        }
    }

    /// Predicts human off-target liabilities for a parsed `chem_core::Molecule`.
    #[must_use]
    pub fn predict_molecule(&self, mol: &Molecule) -> OffTargetReport {
        let (fp, _) = compute_ecfp4(mol);
        self.predict_fp(&fp)
    }
}

fn risk_rank(r: OffTargetRisk) -> u8 {
    match r {
        OffTargetRisk::High => 2,
        OffTargetRisk::Moderate => 1,
        OffTargetRisk::Low => 0,
    }
}

/// Analytical extreme-value distribution (Type-I Gumbel) upper tail p-value (Keiser et al. 2007).
#[inline]
fn gumbel_tail_p(z: f64) -> f64 {
    if z <= -6.0 {
        return 1.0;
    }
    let x = -(z * PI_OVER_SQRT6 + std::f64::consts::EULER_GAMMA);
    if x < -35.0 {
        return x.exp();
    }
    let e = x.exp();
    if e < 1e-6 {
        e
    } else {
        (1.0 - (-e).exp()).clamp(1e-30, 1.0)
    }
}

/// Calibrates empirical background `(mu_bg, sigma_bg)` for each target set by cross-scoring
/// all ligands belonging to *other* targets in the index.
fn calibrate_sea_backgrounds(targets: &mut [SeaTargetSet], tc_cutoff: f64) {
    let all_fps: Vec<(usize, Ecfp4BitVec)> = targets
        .iter()
        .enumerate()
        .flat_map(|(ti, t)| t.ligands.iter().map(move |l| (ti, l.fp)))
        .collect();

    for (ti, tgt) in targets.iter_mut().enumerate() {
        let mut bg_scores = Vec::with_capacity(all_fps.len());
        for &(other_ti, ref qfp) in &all_fps {
            if other_ti == ti {
                continue;
            }
            let mut sum = 0.0_f64;
            for lig in &tgt.ligands {
                let tc = qfp.tanimoto(&lig.fp);
                if tc >= tc_cutoff {
                    sum += tc;
                }
            }
            bg_scores.push(sum);
        }
        if bg_scores.is_empty() {
            tgt.mu_bg = 0.0;
            tgt.sigma_bg = 0.25;
            continue;
        }
        let n = bg_scores.len() as f64;
        let mean = bg_scores.iter().sum::<f64>() / n;
        let var = bg_scores
            .iter()
            .map(|&v| {
                let d = v - mean;
                d * d
            })
            .sum::<f64>()
            / n.max(1.0);
        tgt.mu_bg = (mean * 10000.0).round() / 10000.0;
        tgt.sigma_bg = (var.sqrt().max(0.14) * 10000.0).round() / 10000.0;
    }
}

/// Builds a custom SEA off-target index from a TSV file with header:
/// `target_id\ttarget_name\ttarget_class\torgan_liability\tligand_name\tsmiles`.
pub fn build_sea_index_from_tsv(
    tsv_path: &Path,
    out_json: &Path,
) -> Result<SeaOffTargetIndex, QsarError> {
    let raw = std::fs::read_to_string(tsv_path)?;
    let mut grouped: BTreeMap<String, SeaTargetSet> = BTreeMap::new();
    for (line_idx, line) in raw.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').map(str::trim).collect();
        if cols.len() < 6 {
            continue;
        }
        if line_idx == 0 && cols[0].eq_ignore_ascii_case("target_id") {
            continue;
        }
        let Ok(mol) = parse_smiles(cols[5]) else {
            continue;
        };
        if mol.atoms.is_empty() {
            continue;
        }
        let (fp, _) = compute_ecfp4(&mol);
        let entry = grouped
            .entry(cols[0].to_string())
            .or_insert_with(|| SeaTargetSet {
                target_id: cols[0].to_string(),
                target_name: cols[1].to_string(),
                target_class: cols[2].to_string(),
                organ_liability: cols[3].to_string(),
                mu_bg: 0.0,
                sigma_bg: 0.25,
                ligands: Vec::new(),
            });
        entry.ligands.push(SeaReferenceLigand {
            name: cols[4].to_string(),
            smiles: mol.to_canonical_smiles(),
            fp,
        });
    }

    let mut targets: Vec<SeaTargetSet> = grouped.into_values().collect();
    if targets.is_empty() {
        return Err(QsarError::Dataset {
            endpoint: "sea_offtargets".to_string(),
            message: "no valid target-ligand rows parsed from TSV".to_string(),
        });
    }
    calibrate_sea_backgrounds(&mut targets, SEA_TC_CUTOFF);

    let index = SeaOffTargetIndex {
        format: "admetcraft-sea-index".to_string(),
        version: format!("custom-{}", targets.len()),
        tc_cutoff: SEA_TC_CUTOFF,
        targets,
    };
    index.save_to_path(out_json)?;
    Ok(index)
}

static DEFAULT_SEA_CACHE: OnceLock<SeaOffTargetIndex> = OnceLock::new();

/// Returns the active SEA off-target index: loads user-enriched `CRAFT_DATA_DIR/.../sea_offtargets.json`
/// if present, otherwise falls back cleanly to the embedded Bowes Safety-20 reference panel.
#[must_use]
pub fn load_default_sea_index() -> &'static SeaOffTargetIndex {
    DEFAULT_SEA_CACHE.get_or_init(|| {
        let p = default_sea_index_path();
        if p.is_file()
            && let Ok(idx) = SeaOffTargetIndex::load_from_path(&p)
            && !idx.targets.is_empty()
        {
            return idx;
        }
        builtin_safety20_sea_index()
    })
}

/// Curated offline Bowes Safety-20 human anti-target reference panel (verified clinical drugs & probes).
#[must_use]
pub fn builtin_safety20_sea_index() -> SeaOffTargetIndex {
    type RawTarget<'a> = (&'a str, &'a str, &'a str, &'a str, &'a [(&'a str, &'a str)]);

    const SAFETY20_RAW: &[RawTarget<'static>] = &[
        (
            "hERG",
            "hERG K+ Channel (KCNH2)",
            "Ion Channel",
            "Cardiovascular (QT Prolongation / Torsades de Pointes)",
            &[
                (
                    "Terfenadine",
                    "CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",
                ),
                (
                    "Astemizole",
                    "COc1ccc(CCN2CCC(Nc3nc4ccccc4n3Cc3ccc(F)cc3)CC2)cc1",
                ),
                (
                    "Cisapride",
                    "COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC",
                ),
                (
                    "Dofetilide",
                    "CN(CCc1ccc(NS(C)(=O)=O)cc1)CCOc1ccc(NS(C)(=O)=O)cc1",
                ),
                (
                    "Sertindole",
                    "O=C1NCCN1CCN1CCC(c2cn(-c3ccc(F)cc3)c3ccc(Cl)cc23)CC1",
                ),
                (
                    "Pimozide",
                    "O=c1[nH]c2ccccc2n1C1CCN(CCCC(c2ccc(F)cc2)c2ccc(F)cc2)CC1",
                ),
                (
                    "Haloperidol",
                    "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
                ),
                (
                    "Verapamil",
                    "COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC",
                ),
                (
                    "Amiodarone",
                    "CCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1",
                ),
                ("Quinidine", "COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1"),
                ("Bepridil", "CC(C)COCC(CN(Cc1ccccc1)c1ccccc1)N1CCCC1"),
                (
                    "E4031",
                    "CS(=O)(=O)Nc1ccc(C(=O)CC2CCN(CCc3ccc(NS(C)(=O)=O)cc3)CC2)cc1",
                ),
            ],
        ),
        (
            "5HT2B",
            "Serotonin 5-HT2B Receptor (HTR2B)",
            "GPCR",
            "Cardiovascular (Cardiac Valvulopathy)",
            &[
                ("Fenfluramine", "CCNC(C)Cc1cccc(C(F)(F)F)c1"),
                ("Norfenfluramine", "CC(N)Cc1cccc(C(F)(F)F)c1"),
                ("Pergolide", "CCCN1CC(CSC)CC2c3cccc4[nH]cc(c34)CC21"),
                (
                    "Cabergoline",
                    "CCN(CCCNC(=O)C1CN(CC=C)C2Cc3c[nH]c4cccc(c34)C2C1)C(=O)NC",
                ),
                (
                    "Ergotamine",
                    "CN1CC(C(=O)NC2(C)OC3(O)C4CCCN4C(=O)C(Cc4ccccc4)N3C2=O)C=C2c3cccc4[nH]cc(c34)CC21",
                ),
                (
                    "Methysergide",
                    "CCC(CO)NC(=O)C1CN(C)C2Cc3cn(C)c4cccc(c34)C2=C1",
                ),
                ("MDMA", "CNC(C)Cc1ccc2c(c1)OCO2"),
                (
                    "Bromocriptine",
                    "CC(C)CC1C(=O)N2CCCC2C2(O)OC(C(C)C)(NC(=O)C3CN(C)C4Cc5c(Br)[nH]c6cccc(c56)C4=C3)C(=O)N12",
                ),
                (
                    "Methylergometrine",
                    "CCC(CO)NC(=O)C1CN(C)C2Cc3c[nH]c4cccc(c34)C2=C1",
                ),
                ("Lorcaserin", "CC1CNCCc2ccc(Cl)cc21"),
            ],
        ),
        (
            "5HT2A",
            "Serotonin 5-HT2A Receptor (HTR2A)",
            "GPCR",
            "CNS (Hallucinations / Serotonin Syndrome)",
            &[
                (
                    "Ketanserin",
                    "O=c1[nH]c(=O)n(CCN2CCC(C(=O)c3ccc(F)cc3)CC2)c2ccccc12",
                ),
                (
                    "Risperidone",
                    "Cc1nc2n(c(=O)c1CCN1CCC(c3noc4cc(F)ccc34)CC1)CCCC2",
                ),
                ("Olanzapine", "Cc1cc2c(s1)Nc1ccccc1N=C2N1CCN(C)CC1"),
                ("Clozapine", "CN1CCN(C2=Nc3cc(Cl)ccc3Nc3ccccc32)CC1"),
                ("Trazodone", "O=c1n(CCCN2CCN(c3cccc(Cl)c3)CC2)nc2ccccn12"),
                ("Mirtazapine", "CN1CCN2c3ncccc3Cc3ccccc3C2C1"),
                ("Cyproheptadine", "CN1CCC(=C2c3ccccc3C=Cc3ccccc32)CC1"),
                (
                    "Nefazodone",
                    "CCc1nn(CCCOc2ccccc2)c(=O)n1CCN1CCN(c2cccc(Cl)c2)CC1",
                ),
                ("Quetiapine", "OCCOCCN1CCN(C2=Nc3ccccc3Sc3ccccc32)CC1"),
                (
                    "Ziprasidone",
                    "O=C1Cc2cc(CCN3CCN(c4nsc5ccccc45)CC3)c(Cl)cc2N1",
                ),
            ],
        ),
        (
            "DRD2",
            "Dopamine D2 Receptor (DRD2)",
            "GPCR",
            "CNS (Extrapyramidal Symptoms / Hyperprolactinemia)",
            &[
                (
                    "Haloperidol",
                    "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
                ),
                ("Chlorpromazine", "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21"),
                ("Raclopride", "CCN1CCCC1CNC(=O)c1c(O)c(Cl)cc(Cl)c1OC"),
                ("Sulpiride", "CCN1CCCC1CNC(=O)c1cc(S(N)(=O)=O)ccc1OC"),
                (
                    "Domperidone",
                    "O=c1[nH]c2ccc(Cl)cc2n1CCCN1CCC(n2c(=O)[nH]c3ccccc32)CC1",
                ),
                ("Metoclopramide", "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC"),
                (
                    "Aripiprazole",
                    "O=C1CCc2ccc(OCCCCN3CCN(c4cccc(Cl)c4Cl)CC3)cc2N1",
                ),
                ("Pramipexole", "CCCNc1ccc2nc(N)sc2c1"),
                ("Ropinirole", "CCCN(CCC)CCc1cccc2c1CC(=O)N2"),
                (
                    "Fluphenazine",
                    "OCCN1CCN(CCCN2c3ccccc3Sc3ccc(C(F)(F)F)cc32)CC1",
                ),
            ],
        ),
        (
            "ADRB2",
            "Beta-2 Adrenergic Receptor (ADRB2)",
            "GPCR",
            "Cardiovascular / Pulmonary (Tachycardia / Bronchospasm)",
            &[
                ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12"),
                ("Salbutamol", "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1"),
                ("Terbutaline", "CC(C)(C)NCC(O)c1cc(O)cc(O)c1"),
                ("Salmeterol", "OCc1cc(C(O)CNCCCCCCOCCCCc2ccccc2)ccc1O"),
                ("Formoterol", "COc1ccc(CC(C)NCC(O)c2ccc(O)c(NC=O)c2)cc1"),
                ("Carvedilol", "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12"),
                ("Metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1"),
                ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
                ("Timolol", "CC(C)(C)NCC(O)COc1nsnc1N1CCOCC1"),
                ("Isoprenaline", "CC(C)NCC(O)c1ccc(O)c(O)c1"),
            ],
        ),
        (
            "ADRA1A",
            "Alpha-1A Adrenergic Receptor (ADRA1A)",
            "GPCR",
            "Cardiovascular (Orthostatic Hypotension / Reflex Tachycardia)",
            &[
                ("Prazosin", "COc1cc2nc(N3CCN(C(=O)c4ccco4)CC3)nc(N)c2cc1OC"),
                (
                    "Doxazosin",
                    "COc1cc2nc(N3CCN(C(=O)C4COc5ccccc5O4)CC3)nc(N)c2cc1OC",
                ),
                ("Terazosin", "COc1cc2nc(N3CCN(C(=O)C4CCCO4)CC3)nc(N)c2cc1OC"),
                ("Tamsulosin", "CCOc1ccccc1OCCNC(C)Cc1ccc(OC)c(S(N)(=O)=O)c1"),
                ("Alfuzosin", "COc1cc2nc(N(C)CCCNC(=O)C3CCCO3)nc(N)c2cc1OC"),
                ("Phenylephrine", "CNCC(O)c1cccc(O)c1"),
                ("Oxymetazoline", "Cc1cc(C(C)(C)C)c(O)c(C)c1CC1=NCCN1"),
                ("Phentolamine", "Cc1ccc(N(CC2=NCCN2)c2cccc(O)c2)cc1"),
                (
                    "Silodosin",
                    "CC(Cc1cc(C(N)=O)c2c(c1)CCN2CCCO)NCCOc1ccccc1OCC(F)(F)F",
                ),
            ],
        ),
        (
            "CHRM1",
            "Muscarinic M1 Receptor (CHRM1)",
            "GPCR",
            "CNS / Autonomic (Anticholinergic Delirium / Dry Mouth / Tachycardia)",
            &[
                ("Atropine", "CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2"),
                ("Scopolamine", "CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C1OC21"),
                ("Oxybutynin", "CCN(CC)CC#CCOC(=O)C(O)(c1ccccc1)C1CCCCC1"),
                ("Tolterodine", "Cc1ccc(O)c(C(CCN(C(C)C)C(C)C)c2ccccc2)c1"),
                ("Solifenacin", "O=C(OC1CN2CCC1CC2)N1CCc2ccccc2C1c1ccccc1"),
                (
                    "Ipratropium",
                    "CC(C)[N+]1(C)C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2",
                ),
                ("Pirenzepine", "CN1CCN(CC(=O)N2c3ccccc3C(=O)Nc3cccnc32)CC1"),
                ("Trihexyphenidyl", "OC(CCN1CCCCC1)(c1ccccc1)C1CCCCC1"),
                ("Benztropine", "CN1C2CCC1CC(OC(c1ccccc1)c1ccccc1)C2"),
                ("Amitriptyline", "CN(C)CCC=C1c2ccccc2CCc2ccccc21"),
            ],
        ),
        (
            "HRH1",
            "Histamine H1 Receptor (HRH1)",
            "GPCR",
            "CNS (Sedation / Weight Gain)",
            &[
                ("Diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1"),
                ("Cetirizine", "O=C(O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1"),
                ("Loratadine", "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1"),
                (
                    "Fexofenadine",
                    "CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",
                ),
                ("Chlorpheniramine", "CN(C)CCC(c1ccc(Cl)cc1)c1ccccn1"),
                ("Promethazine", "CC(CN1c2ccccc2Sc2ccccc21)N(C)C"),
                ("Hydroxyzine", "OCCOCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1"),
                ("Doxepin", "CN(C)CCC=C1c2ccccc2COc2ccccc21"),
                ("Desloratadine", "Clc1ccc2c(c1)CCc1cccnc1C2=C1CCNCC1"),
                ("Ketotifen", "CN1CCC(=C2c3ccccc3CC(=O)c3sccc32)CC1"),
            ],
        ),
        (
            "OPRM1",
            "Mu-Opioid Receptor (OPRM1)",
            "GPCR",
            "CNS / Respiratory (Respiratory Depression / Abuse Liability)",
            &[
                ("Morphine", "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"),
                ("Codeine", "COc1ccc2c3c1OC1C(O)C=CC4C(C2)N(C)CCC314"),
                ("Fentanyl", "CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1"),
                ("Methadone", "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1"),
                ("Tramadol", "COc1cccc(C2(O)CCCCC2CN(C)C)c1"),
                ("Oxycodone", "COc1ccc2c3c1OC1C(=O)CCC4(O)C(C2)N(C)CCC314"),
                ("Naloxone", "C=CCN1CCC23c4c5ccc(O)c4OC2C(=O)CCC3(O)C1C5"),
                ("Naltrexone", "O=C1CCC2(O)C3Cc4ccc(O)c5c4C2(CCN3CC2CC2)C1O5"),
                (
                    "Loperamide",
                    "CN(C)C(=O)C(CCN1CCC(O)(c2ccc(Cl)cc2)CC1)(c1ccccc1)c1ccccc1",
                ),
                (
                    "Buprenorphine",
                    "COc1ccc2c3c1OC1C(OC)(C(O)(C)C(C)(C)C)CCC4(C1)C(C2)N(CC1CC1)CCC34",
                ),
            ],
        ),
        (
            "COX1",
            "Cyclooxygenase-1 (PTGS1)",
            "Enzyme",
            "GI / Hematologic (Gastric Ulceration / Bleeding)",
            &[
                ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
                ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
                ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
                (
                    "Indomethacin",
                    "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                ),
                ("Ketoprofen", "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1"),
                ("Flurbiprofen", "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1"),
                ("Piroxicam", "CN1C(=C(O)c2ccccc2S1(=O)=O)C(=O)Nc1ccccn1"),
                ("Ketorolac", "O=C(O)C1CCn2c1ccc2C(=O)c1ccccc1"),
                ("Diclofenac", "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
            ],
        ),
        (
            "COX2",
            "Cyclooxygenase-2 (PTGS2)",
            "Enzyme",
            "Cardiovascular / Renal (Thrombotic Risk / Nephrotoxicity)",
            &[
                (
                    "Celecoxib",
                    "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
                ),
                ("Rofecoxib", "CS(=O)(=O)c1ccc(C2=C(c3ccccc3)C(=O)OC2)cc1"),
                ("Valdecoxib", "Cc1onc(-c2ccccc2)c1-c1ccc(S(N)(=O)=O)cc1"),
                (
                    "Etoricoxib",
                    "Cc1ccc(-c2ncc(Cl)cc2-c2ccc(S(C)(=O)=O)cc2)cn1",
                ),
                ("Meloxicam", "Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1"),
                ("Nimesulide", "CS(=O)(=O)Nc1ccc([N+](=O)[O-])cc1Oc1ccccc1"),
                ("Lumiracoxib", "Cc1ccc(Nc2c(F)cccc2Cl)c(CC(=O)O)c1"),
                ("Diclofenac", "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
            ],
        ),
        (
            "ACHE",
            "Acetylcholinesterase (ACHE)",
            "Enzyme",
            "CNS / Neuromuscular (Cholinergic Crisis / Bradycardia)",
            &[
                ("Donepezil", "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2"),
                ("Rivastigmine", "CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1"),
                ("Galantamine", "COc1ccc2c3c1OC1CC(O)C=CC31CCN(C)C2"),
                ("Tacrine", "Nc1c2c(nc3ccccc13)CCCC2"),
                ("Physostigmine", "CNC(=O)Oc1ccc2c(c1)C1(C)CCN(C)C1N2C"),
                ("Pyridostigmine", "CN(C)C(=O)Oc1ccc[n+](C)c1"),
                ("Neostigmine", "CN(C)C(=O)Oc1cccc([N+](C)(C)C)c1"),
                ("Huperzine_A", "CC=C1C2Cc3ccc(=O)[nH]c3C1(N)CC(C)=C2"),
            ],
        ),
        (
            "PXR",
            "Pregnane X Receptor (NR1I2)",
            "Nuclear Receptor",
            "Hepatic (CYP3A4 Induction / Severe Drug-Drug Interactions)",
            &[
                (
                    "Rifampicin",
                    "COc1c(C)c2c(c(O)c1C=NN1CCN(C)CC1)C(=O)C1(C)OC=CC(OC)C(C)C(O)C(C)C(O)C(C)C(OC(C)=O)C(C)C=CC=C(C)C(=O)NC2=O",
                ),
                (
                    "SR12813",
                    "CC(C)(C)c1cc(C=C(P(=O)(OCC)OCC)P(=O)(OCC)OCC)cc(C(C)(C)C)c1O",
                ),
                ("Carbamazepine", "NC(=O)N1c2ccccc2C=Cc2ccccc21"),
                ("Phenytoin", "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1"),
                ("Phenobarbital", "CCC1(c2ccccc2)C(=O)NC(=O)NC1=O"),
                (
                    "Ritonavir",
                    "CC(C)c1nc(CN(C)C(=O)NC(C(C)C)C(=O)NC(Cc2ccccc2)CC(O)C(Cc2ccccc2)NC(=O)OCc2cncs2)cs1",
                ),
                ("Clotrimazole", "Clc1ccccc1C(c1ccccc1)(c1ccccc1)n1ccnc1"),
                ("Efavirenz", "O=C1Nc2ccc(Cl)cc2C(C#CC2CC2)(C(F)(F)F)O1"),
                (
                    "Bosentan",
                    "COc1ccccc1Oc1c(NS(=O)(=O)c2ccc(C(C)(C)C)cc2)nc(-c2ncccn2)nc1OCCO",
                ),
            ],
        ),
        (
            "PPARG",
            "Peroxisome Proliferator-Activated Receptor Gamma (PPARG)",
            "Nuclear Receptor",
            "Cardiovascular / Metabolic (Fluid Retention / Heart Failure / Weight Gain)",
            &[
                ("Rosiglitazone", "CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1"),
                ("Pioglitazone", "CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1"),
                (
                    "Troglitazone",
                    "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
                ),
                ("Ciglitazone", "O=C1NC(=O)C(Cc2ccc(OCC3(C)CCCCC3)cc2)S1"),
                (
                    "Telmisartan",
                    "CCCc1nc2c(C)cc(-c3nc4ccccc4n3C)cc2n1Cc1ccc(-c2ccccc2C(=O)O)cc1",
                ),
                ("GW9662", "O=[N+]([O-])c1ccc(Cl)c(C(=O)Nc2ccccc2)c1"),
                (
                    "Lobeglitazone",
                    "COc1cc(OCCOc2ccc(CC3SC(=O)NC3=O)cc2)nc(N(C)c2ccccn2)n1",
                ),
                (
                    "Indomethacin",
                    "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                ),
            ],
        ),
        (
            "ESR1",
            "Estrogen Receptor Alpha (ESR1)",
            "Nuclear Receptor",
            "Endocrine (Reproductive / Proliferative Endocrine Disruption)",
            &[
                ("Estradiol", "CC12CCC3c4ccc(O)cc4CCC3C1CCC2O"),
                ("Ethinylestradiol", "C#CC1(O)CCC2C3CCc4cc(O)ccc4C3CCC21C"),
                ("Diethylstilbestrol", "CCC(=C(CC)c1ccc(O)cc1)c1ccc(O)cc1"),
                ("Tamoxifen", "CC/C(=C(/c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
                (
                    "Raloxifene",
                    "O=C(c1ccc(OCCN2CCCCC2)cc1)c1c(-c2ccc(O)cc2)sc2cc(O)ccc12",
                ),
                ("Bisphenol_A", "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"),
                ("Genistein", "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12"),
                (
                    "Fulvestrant",
                    "CC12CCC3c4ccc(O)cc4CC(CCCCCCCCC[S](=O)CCCC(F)(F)C(F)(F)F)C3C1CCC2O",
                ),
                ("Clomifene", "CCN(CC)CCOc1ccc(C(=C(Cl)c2ccccc2)c2ccccc2)cc1"),
            ],
        ),
        (
            "AR",
            "Androgen Receptor (AR)",
            "Nuclear Receptor",
            "Endocrine (Anti-Androgenic / Reproductive Toxicity)",
            &[
                ("Testosterone", "CC12CCC3C(CCC4=CC(=O)CCC43C)C1CCC2O"),
                ("DHT", "CC12CCC3C(CCC4CC(=O)CCC43C)C1CCC2O"),
                (
                    "Bicalutamide",
                    "CC(O)(CS(=O)(=O)c1ccc(F)cc1)C(=O)Nc1ccc(C#N)c(C(F)(F)F)c1",
                ),
                ("Flutamide", "CC(C)C(=O)Nc1ccc([N+](=O)[O-])c(C(F)(F)F)c1"),
                (
                    "Enzalutamide",
                    "CNC(=O)c1ccc(N2C(=S)N(c3ccc(C#N)c(C(F)(F)F)c3)C(=O)C2(C)C)cc1F",
                ),
                (
                    "Nilutamide",
                    "CC1(C)NC(=O)N(c2ccc([N+](=O)[O-])c(C(F)(F)F)c2)C1=O",
                ),
                (
                    "Apalutamide",
                    "CNC(=O)c1ccc(N2C(=S)N(c3ccc(C#N)c(C(F)(F)F)n3)C(=O)C23CCC3)cc1F",
                ),
                (
                    "Spironolactone",
                    "CC(=O)SC1CC2=CC(=O)CCC2(C)C2CCC3(C)C(CCC34CCC(=O)O4)C12",
                ),
            ],
        ),
        (
            "NR3C1",
            "Glucocorticoid Receptor (NR3C1)",
            "Nuclear Receptor",
            "Endocrine / Immune (Immunosuppression / Metabolic Syndrome)",
            &[
                (
                    "Dexamethasone",
                    "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO",
                ),
                (
                    "Prednisolone",
                    "CC12C=CC(=O)C=C1CCC1C2C(O)CC2(C)C1CCC2(O)C(=O)CO",
                ),
                (
                    "Hydrocortisone",
                    "CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C1CCC2(O)C(=O)CO",
                ),
                (
                    "Betamethasone",
                    "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO",
                ),
                (
                    "Triamcinolone",
                    "CC12C=CC(=O)C=C1CCC1C2(F)C(O)CC2(C)C1CC(O)C2(O)C(=O)CO",
                ),
                (
                    "Budesonide",
                    "CCCC1OC2CC3C4CCC5=CC(=O)C=CC5(C)C4C(O)CC3(C)C2(C(=O)CO)O1",
                ),
                (
                    "Mifepristone",
                    "CC#CC1(O)CCC2C3CCC4=CC(=O)CCC4=C3C(c3ccc(N(C)C)cc3)CC21C",
                ),
                (
                    "Fluticasone",
                    "CC1CC2C3CC(F)C4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(OC(=O)CC)C(=O)SCF",
                ),
            ],
        ),
        (
            "SLC6A3_DAT",
            "Dopamine Transporter (DAT / SLC6A3)",
            "Transporter",
            "CNS / Cardiovascular (Stimulant Abuse Liability / Hypertension)",
            &[
                ("Cocaine", "COC(=O)C1C(OC(=O)c2ccccc2)CC2CCC1N2C"),
                ("Methylphenidate", "COC(=O)C(c1ccccc1)C1CCCCN1"),
                ("Bupropion", "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1"),
                ("Amphetamine", "CC(N)Cc1ccccc1"),
                ("Methamphetamine", "CNC(C)Cc1ccccc1"),
                ("Modafinil", "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1"),
                ("Nomifensine", "CN1Cc2c(N)cccc2C(c2ccccc2)C1"),
                ("Mazindol", "OC1(c2ccc(Cl)cc2)c2ccccc2C2=NCCN21"),
                (
                    "GBR12909",
                    "Fc1ccc(C(OCCN2CCN(CCCc3ccccc3)CC2)c2ccc(F)cc2)cc1",
                ),
            ],
        ),
        (
            "SLC6A4_SERT",
            "Serotonin Transporter (SERT / SLC6A4)",
            "Transporter",
            "CNS (Serotonin Syndrome / Drug Interactions)",
            &[
                ("Fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"),
                ("Sertraline", "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21"),
                ("Paroxetine", "Fc1ccc(C2CCNCC2COc2ccc3c(c2)OCO3)cc1"),
                ("Citalopram", "CN(C)CCCC1(c2ccc(F)cc2)OCc2cc(C#N)ccc21"),
                ("Fluvoxamine", "COCCCCC(=NOCCN)c1ccc(C(F)(F)F)cc1"),
                ("Venlafaxine", "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1"),
                ("Duloxetine", "CNCCC(Oc1cccc2ccccc12)c1cccs1"),
                ("Imipramine", "CN(C)CCCN1c2ccccc2CCc2ccccc21"),
                ("Clomipramine", "CN(C)CCCN1c2ccccc2CCc2ccc(Cl)cc21"),
                ("Trazodone", "O=c1n(CCCN2CCN(c3cccc(Cl)c3)CC2)nc2ccccn12"),
            ],
        ),
        (
            "CYP3A4_STRONG",
            "CYP3A4 Strong Inhibitor / Mechanism-Based Inactivator",
            "Metabolic Enzyme",
            "Hepatic (Severe Clinical Drug-Drug Interactions)",
            &[
                (
                    "Ketoconazole",
                    "CC(=O)N1CCN(c2ccc(OCC3COC(Cn4ccnc4)(c4ccc(Cl)cc4Cl)O3)cc2)CC1",
                ),
                (
                    "Itraconazole",
                    "CCC(C)n1ncn(-c2ccc(N3CCN(c4ccc(OCC5COC(Cn6cncn6)(c6ccc(Cl)cc6Cl)O5)cc4)CC3)cc2)c1=O",
                ),
                ("Voriconazole", "CC(c1ncncc1F)C(O)(Cn1cncn1)c1ccc(F)cc1F"),
                ("Fluconazole", "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"),
                (
                    "Clarithromycin",
                    "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(OC)CC(C)C(=O)C(C)C(O)C1(C)O",
                ),
                (
                    "Ritonavir",
                    "CC(C)c1nc(CN(C)C(=O)NC(C(C)C)C(=O)NC(Cc2ccccc2)CC(O)C(Cc2ccccc2)NC(=O)OCc2cncs2)cs1",
                ),
                (
                    "Cobicistat",
                    "CC(C)c1nc(CN(C)C(=O)NC(CCN2CCOCC2)C(=O)NC(Cc2ccccc2)CCC(Cc2ccccc2)NC(=O)OCc2cncs2)cs1",
                ),
                (
                    "Nefazodone",
                    "CCc1nn(CCCOc2ccccc2)c(=O)n1CCN1CCN(c2cccc(Cl)c2)CC1",
                ),
                (
                    "Diltiazem",
                    "COc1ccc(C2Sc3ccccc3N(CCN(C)C)C(=O)C2OC(C)=O)cc1",
                ),
                ("Miconazole", "Clc1ccc(COC(Cn2ccnc2)c2ccc(Cl)cc2Cl)c(Cl)c1"),
            ],
        ),
    ];

    let mut targets = Vec::with_capacity(SAFETY20_RAW.len());
    for &(tid, tname, tclass, organ, ligs) in SAFETY20_RAW {
        let mut parsed = Vec::with_capacity(ligs.len());
        for &(lname, smi) in ligs {
            if let Ok(mol) = parse_smiles(smi)
                && !mol.atoms.is_empty()
            {
                let (fp, _) = compute_ecfp4(&mol);
                parsed.push(SeaReferenceLigand {
                    name: lname.to_string(),
                    smiles: mol.to_canonical_smiles(),
                    fp,
                });
            }
        }
        targets.push(SeaTargetSet {
            target_id: tid.to_string(),
            target_name: tname.to_string(),
            target_class: tclass.to_string(),
            organ_liability: organ.to_string(),
            mu_bg: 0.0,
            sigma_bg: 0.25,
            ligands: parsed,
        });
    }
    calibrate_sea_backgrounds(&mut targets, SEA_TC_CUTOFF);

    SeaOffTargetIndex {
        format: "admetcraft-sea-index".to_string(),
        version: "bowes-safety20-v1.1".to_string(),
        tc_cutoff: SEA_TC_CUTOFF,
        targets,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_safety20_parses_cleanly_and_detects_known_off_targets()
    -> Result<(), Box<dyn std::error::Error>> {
        let idx = builtin_safety20_sea_index();
        assert_eq!(idx.targets.len(), 20);
        for t in &idx.targets {
            assert!(
                t.ligands.len() >= 8,
                "target {} should have >= 8 parsed reference ligands, got {}",
                t.target_id,
                t.ligands.len()
            );
        }

        // Haloperidol analog -> must flag hERG & DRD2 with High risk
        let halo_analog = parse_smiles("O=C(CCCN1CCC(O)(c2ccc(Br)cc2)CC1)c1ccc(F)cc1")?;
        let rep_halo = idx.predict_molecule(&halo_analog);
        assert!(
            rep_halo.high_risk_count >= 1,
            "expected High off-target risk for haloperidol analog, got {rep_halo:?}"
        );
        assert!(
            rep_halo
                .hits
                .iter()
                .any(|h| (h.target_id == "hERG" || h.target_id == "DRD2")
                    && h.risk == OffTargetRisk::High),
            "expected hERG or DRD2 High risk hit, got {:?}",
            rep_halo.hits
        );

        // Simple inert fragment (Urea) -> must NOT raise any High or Moderate off-target alerts
        let urea = parse_smiles("NC(=O)N")?;
        let rep_urea = idx.predict_molecule(&urea);
        assert_eq!(rep_urea.high_risk_count, 0);
        assert_eq!(rep_urea.moderate_risk_count, 0);
        Ok(())
    }
}
