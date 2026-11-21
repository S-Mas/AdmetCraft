#![forbid(unsafe_code)]

use admet_calc::{
    AbsorptionLevel, HergRisk, MoleculeEvaluation, QsarAdmetReport, TppCriterionResult,
    TppMpoReport, TppPreset, TppStatus,
};

fn qsar_val(qsar: Option<&QsarAdmetReport>, id: &str) -> Option<f64> {
    qsar.and_then(|q| {
        q.endpoints
            .iter()
            .find(|e| e.endpoint_id == id)
            .map(|e| e.value)
    })
}

fn qsar_linear(qsar: Option<&QsarAdmetReport>, id: &str) -> Option<f64> {
    qsar.and_then(|q| {
        q.endpoints
            .iter()
            .find(|e| e.endpoint_id == id)
            .and_then(|e| e.linear_value)
    })
}

/// Linear ramp desirability: `1.0` when `x <= good`, `0.0` when `x >= bad`.
fn ramp_low_better(x: f64, good: f64, bad: f64) -> f64 {
    if x <= good {
        1.0
    } else if x >= bad {
        0.0
    } else {
        (bad - x) / (bad - good)
    }
}

/// Linear ramp desirability: `1.0` when `x >= good`, `0.0` when `x <= bad`.
fn ramp_high_better(x: f64, bad: f64, good: f64) -> f64 {
    if x >= good {
        1.0
    } else if x <= bad {
        0.0
    } else {
        (x - bad) / (good - bad)
    }
}

/// Trapezoid band desirability: `1.0` inside `[lo_opt, hi_opt]`, tapering to `0.0` outside `[lo_bad, hi_bad]`.
fn band_desirability(x: f64, lo_bad: f64, lo_opt: f64, hi_opt: f64, hi_bad: f64) -> f64 {
    if x < lo_opt {
        ramp_high_better(x, lo_bad, lo_opt)
    } else if x > hi_opt {
        ramp_low_better(x, hi_opt, hi_bad)
    } else {
        1.0
    }
}

fn status_from_d(d: f64) -> TppStatus {
    if d >= 0.70 {
        TppStatus::Pass
    } else if d >= 0.35 {
        TppStatus::Warn
    } else {
        TppStatus::Fail
    }
}

fn push_crit(
    out: &mut Vec<TppCriterionResult>,
    name: &str,
    observed: String,
    target_spec: &str,
    desirability: f64,
    weight: f64,
) {
    let d = (desirability.clamp(0.0, 1.0) * 1000.0).round() / 1000.0;
    out.push(TppCriterionResult {
        name: name.to_string(),
        observed,
        target_spec: target_spec.to_string(),
        desirability: d,
        weight,
        status: status_from_d(d),
    });
}

/// Computes a configurable pharmaceutical Target Product Profile (`TPP MPO`) scorecard.
#[must_use]
pub fn compute_tpp_mpo(eval: &MoleculeEvaluation, preset: TppPreset) -> TppMpoReport {
    let qsar = eval.qsar.as_ref();
    let mut criteria = Vec::with_capacity(12);
    let mut hard_stops = Vec::new();

    // 1. hERG Cardiotoxicity (QSAR probability preferred, heuristic fallback)
    let herg_p = qsar_val(qsar, "herg_karim").unwrap_or(match eval.safety.herg_risk {
        HergRisk::Low => 0.15,
        HergRisk::Moderate => 0.50,
        HergRisk::High => 0.85,
    });
    let d_herg = ramp_low_better(herg_p, 0.35, 0.75);
    push_crit(
        &mut criteria,
        "hERG Cardiotoxicity",
        format!("P = {:.2}", herg_p),
        "P < 0.35 (hard stop > 0.75)",
        d_herg,
        2.0,
    );
    if herg_p >= 0.75 {
        hard_stops.push(format!("High hERG liability (P = {:.2})", herg_p));
    }

    // 2. AMES Mutagenicity
    let ames_p = qsar_val(qsar, "ames").unwrap_or(if eval.safety.ames_alerts.is_empty() {
        0.15
    } else {
        0.80
    });
    let d_ames = ramp_low_better(ames_p, 0.35, 0.75);
    push_crit(
        &mut criteria,
        "AMES Mutagenicity",
        format!("P = {:.2}", ames_p),
        "P < 0.35 (hard stop > 0.75)",
        d_ames,
        1.8,
    );
    if ames_p >= 0.75 {
        hard_stops.push(format!("High AMES mutagenicity risk (P = {:.2})", ames_p));
    }

    // 3. DILI / Clinical Toxicity
    let dili_p = qsar_val(qsar, "dili").unwrap_or(if eval.safety.dili_alerts.is_empty() {
        0.25
    } else {
        0.75
    });
    let clintox_p = qsar_val(qsar, "clintox").unwrap_or(0.10);
    let d_tox =
        0.60 * ramp_low_better(dili_p, 0.40, 0.80) + 0.40 * ramp_low_better(clintox_p, 0.25, 0.65);
    push_crit(
        &mut criteria,
        "DILI / ClinTox Safety",
        format!("DILI {:.2} | ClinTox {:.2}", dili_p, clintox_p),
        "DILI < 0.40 & ClinTox < 0.25",
        d_tox,
        1.6,
    );
    if clintox_p >= 0.65 {
        hard_stops.push(format!("High ClinTox alert (P = {:.2})", clintox_p));
    }

    // 4. CYP450 DDI Liability (max / count across 5 isoforms)
    let cyp_ids = [
        "cyp1a2_veith",
        "cyp2c9_veith",
        "cyp2c19_veith",
        "cyp2d6_veith",
        "cyp3a4_veith",
    ];
    let cyp_probs: Vec<f64> = cyp_ids.iter().filter_map(|id| qsar_val(qsar, id)).collect();
    let (cyp_high_cnt, cyp_max) = if !cyp_probs.is_empty() {
        let cnt = cyp_probs.iter().filter(|&&p| p >= 0.60).count();
        let mx = cyp_probs.iter().copied().fold(0.0_f64, f64::max);
        (cnt, mx)
    } else {
        let pk = &eval.pharmacokinetics;
        let cnt = usize::from(pk.cyp1a2_inhibitor)
            + usize::from(pk.cyp2c9_inhibitor)
            + usize::from(pk.cyp2c19_inhibitor)
            + usize::from(pk.cyp2d6_inhibitor)
            + usize::from(pk.cyp3a4_inhibitor);
        (cnt, if cnt > 0 { 0.65 } else { 0.20 })
    };
    let d_cyp = ramp_low_better(cyp_high_cnt as f64, 1.0, 4.0) * 0.65
        + ramp_low_better(cyp_max, 0.50, 0.90) * 0.35;
    push_crit(
        &mut criteria,
        "CYP450 DDI Risk (5 isoforms)",
        format!("{cyp_high_cnt}/5 high (max P={cyp_max:.2})"),
        "<= 1 isoform > 0.60",
        d_cyp,
        1.2,
    );

    // 5. Absorption & Permeability (HIA + Caco-2 + P-gp)
    let hia_p = qsar_val(qsar, "hia_hou").unwrap_or(
        if eval.pharmacokinetics.gi_absorption == AbsorptionLevel::High {
            0.85
        } else {
            0.25
        },
    );
    let caco2 = qsar_val(qsar, "caco2_wang").unwrap_or(-5.15);
    let pgp_p =
        qsar_val(qsar, "pgp_broccatelli").unwrap_or(if eval.pharmacokinetics.pgp_substrate {
            0.65
        } else {
            0.25
        });
    let d_abs = 0.50 * ramp_high_better(hia_p, 0.40, 0.80)
        + 0.35 * ramp_high_better(caco2, -6.0, -5.0)
        + 0.15 * ramp_low_better(pgp_p, 0.40, 0.80);
    push_crit(
        &mut criteria,
        "Oral Absorption & Caco-2",
        format!("HIA {:.2} | Caco2 {:.2} | P-gp {:.2}", hia_p, caco2, pgp_p),
        "HIA >= 0.80 & Caco2 >= -5.0",
        d_abs,
        match preset {
            TppPreset::CnsPenetrant => 1.3,
            _ => 1.5,
        },
    );

    // 6. Aqueous Solubility (AqSolDB QSAR if available, else ESOL)
    let log_s = qsar_val(qsar, "solubility_aqsoldb").unwrap_or(eval.water_solubility.esol.log_s);
    let d_sol = ramp_high_better(log_s, -6.0, -3.5);
    push_crit(
        &mut criteria,
        "Aqueous Solubility (LogS)",
        format!("{:.2} log(mol/L)", log_s),
        "LogS >= -3.5 (> -6.0)",
        d_sol,
        1.2,
    );

    // 7. Metabolic Clearance (Microsomal CLint) & Half-Life
    let clint = qsar_linear(qsar, "clearance_microsome_az").unwrap_or(25.0);
    let t_half = qsar_linear(qsar, "half_life_obach").unwrap_or(4.5);
    let d_pk =
        0.65 * ramp_low_better(clint, 20.0, 90.0) + 0.35 * ramp_high_better(t_half, 1.0, 4.0);
    push_crit(
        &mut criteria,
        "Metabolic Stability (CLint / t1/2)",
        format!("CLint {:.1} uL/min/mg | t1/2 {:.1} h", clint, t_half),
        "CLint <= 20 & t1/2 >= 4 h",
        d_pk,
        1.2,
    );

    // 8. Brain Penetration (BBB / P-gp) — Preset-specific direction!
    let bbb_p = qsar_val(qsar, "bbb_martins").unwrap_or(if eval.pharmacokinetics.bbb_permeant {
        0.80
    } else {
        0.25
    });
    match preset {
        TppPreset::CnsPenetrant => {
            let d_bbb = 0.65 * ramp_high_better(bbb_p, 0.40, 0.80)
                + 0.35 * ramp_low_better(pgp_p, 0.30, 0.70);
            push_crit(
                &mut criteria,
                "CNS BBB+ & Low P-gp Efflux",
                format!("BBB P={:.2} | P-gp P={:.2}", bbb_p, pgp_p),
                "BBB >= 0.80 & P-gp <= 0.30",
                d_bbb,
                2.2,
            );
        }
        TppPreset::PeripheralSelective => {
            let d_periph = ramp_low_better(bbb_p, 0.30, 0.75);
            push_crit(
                &mut criteria,
                "Peripheral Restriction (BBB-)",
                format!("BBB P={:.2}", bbb_p),
                "BBB <= 0.30 (avoid CNS exposure)",
                d_periph,
                1.8,
            );
        }
        _ => {}
    }

    // 9. Physicochemical Envelope (MW, LogD7.4, TPSA) — Preset-specific!
    let mw = eval.physicochemical.molecular_weight;
    let logd = eval.ionization.logd_ph74;
    let tpsa = eval.physicochemical.tpsa;
    let (mw_hi, mw_bad, logd_lo, logd_hi, tpsa_hi, tpsa_bad, spec) = match preset {
        TppPreset::OralSmallMolecule => (
            450.0,
            580.0,
            0.5,
            3.5,
            120.0,
            160.0,
            "MW<=450, LogD 0.5-3.5, TPSA<=120",
        ),
        TppPreset::CnsPenetrant => (
            380.0,
            480.0,
            1.0,
            3.5,
            80.0,
            110.0,
            "MW<=380, LogD 1.0-3.5, TPSA<=80",
        ),
        TppPreset::PeripheralSelective => (
            500.0,
            620.0,
            0.0,
            3.5,
            140.0,
            180.0,
            "MW<=500, LogD 0.0-3.5, TPSA<=140",
        ),
        TppPreset::OncologyKinase => (
            520.0,
            650.0,
            1.0,
            4.2,
            135.0,
            175.0,
            "MW<=520, LogD 1.0-4.2, TPSA<=135",
        ),
    };
    let d_physchem = 0.35 * ramp_low_better(mw, mw_hi, mw_bad)
        + 0.35 * band_desirability(logd, logd_lo - 2.0, logd_lo, logd_hi, logd_hi + 2.0)
        + 0.30 * ramp_low_better(tpsa, tpsa_hi, tpsa_bad);
    push_crit(
        &mut criteria,
        "Physicochemical Envelope",
        format!("MW {:.0} | LogD {:.2} | TPSA {:.0}", mw, logd, tpsa),
        spec,
        d_physchem,
        1.1,
    );

    // 10. Structural Alerts (PAINS + Brenk)
    let n_pains = eval.medicinal_chemistry.pains_alerts.len();
    let n_brenk = eval.medicinal_chemistry.brenk_alerts.len();
    let d_alerts = if n_pains == 0 && n_brenk == 0 {
        1.0
    } else if n_pains == 0 && n_brenk == 1 {
        0.55
    } else {
        0.10
    };
    push_crit(
        &mut criteria,
        "MedChem Structural Alerts",
        format!("PAINS={n_pains}, Brenk={n_brenk}"),
        "0 PAINS & 0 Brenk",
        d_alerts,
        1.0,
    );
    if n_pains > 0 {
        hard_stops.push(format!(
            "PAINS structural alert ({})",
            eval.medicinal_chemistry.pains_alerts.join(", ")
        ));
    }

    // 11. Tox21 Cellular Safety (12 pathways: 7 Nuclear Receptors + 5 Stress Responses)
    if let Some(tox21) = &eval.tox21_summary {
        let total_alerts = tox21.endocrine_alerts.len() + tox21.stress_alerts.len();
        let d_tox21 = 0.60 * ramp_low_better(total_alerts as f64, 0.0, 3.0)
            + 0.40 * ramp_low_better(tox21.max_probability, 0.35, 0.80);
        push_crit(
            &mut criteria,
            "Tox21 Cellular & Endocrine Pathways",
            format!(
                "{} alert(s) (NR={}, SR={}, max P={:.2})",
                total_alerts,
                tox21.endocrine_alerts.len(),
                tox21.stress_alerts.len(),
                tox21.max_probability
            ),
            "0 Tox21 alerts & max P < 0.35",
            d_tox21,
            1.2,
        );
        if total_alerts >= 3 || (tox21.stress_alerts.len() >= 2 && tox21.max_probability >= 0.75) {
            hard_stops.push(format!(
                "Multiple Tox21 cellular/endocrine pathway alerts ({total_alerts} flagged, max P={:.2})",
                tox21.max_probability
            ));
        }
    }

    // 12. Off-Target Safety Pharmacology (SEA) & Phase I/II Metabolite Bioactivation Guard
    if eval.off_targets.is_some() || eval.metabolite_guard.is_some() {
        let sea_high = eval.off_targets.as_ref().map_or(0, |o| o.high_risk_count);
        let sea_mod = eval
            .off_targets
            .as_ref()
            .map_or(0, |o| o.moderate_risk_count);
        let met_bio = eval
            .metabolite_guard
            .as_ref()
            .map_or(0, |m| m.bioactivation_alerts);
        let met_esc = eval
            .metabolite_guard
            .as_ref()
            .map_or(0, |m| m.escalation_count);

        let d_sea = ramp_low_better((sea_high as f64) + 0.35 * (sea_mod as f64), 0.0, 2.5);
        let d_met = ramp_low_better((met_bio as f64) + 0.40 * (met_esc as f64), 0.0, 2.5);
        let d_guard = 0.50 * d_sea + 0.50 * d_met;
        push_crit(
            &mut criteria,
            "Off-Target SEA & Metabolite Guard",
            format!("SEA High={sea_high}/Mod={sea_mod} | Bioact={met_bio}/Esc={met_esc}"),
            "0 SEA High & 0 reactive bioactivation",
            d_guard,
            1.3,
        );
        if met_bio >= 2 {
            hard_stops.push(format!(
                "Multiple reactive metabolite bioactivation pathways ({met_bio})"
            ));
        }
        if sea_high >= 2 {
            hard_stops.push(format!(
                "Promiscuous Bowes Safety-20 off-target profile ({sea_high} High-risk targets)"
            ));
        }
    }

    let w_sum: f64 = criteria.iter().map(|c| c.weight).sum();
    let raw_score = if w_sum > 0.0 {
        criteria
            .iter()
            .map(|c| c.desirability * c.weight)
            .sum::<f64>()
            / w_sum
    } else {
        0.0
    };
    let score = (raw_score * 1000.0).round() / 1000.0;
    let decision = if !hard_stops.is_empty() || score < 0.45 {
        "Hold / Red-Flag".to_string()
    } else if score >= 0.72 {
        "Advance".to_string()
    } else {
        "Optimize".to_string()
    };

    TppMpoReport {
        preset,
        score,
        decision,
        hard_stops,
        criteria,
    }
}
