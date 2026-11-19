#![forbid(unsafe_code)]

use admet_calc::{ApplicabilityDomainStatus, QsarTaskKind, TppPreset};
use admet_qsar::{
    AdmetModelPack, default_model_pack_path, evaluate_smiles_full, load_default_model_pack,
};

#[test]
fn tdc_model_pack_validates_metrics_conformal_coverage_and_known_drugs()
-> Result<(), Box<dyn std::error::Error>> {
    let pack_path = default_model_pack_path();
    if !pack_path.is_file() {
        // Clean clone without locally trained TDC pack: verify fallback path runs cleanly
        let ev = evaluate_smiles_full(
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
            Some("Ibuprofen"),
            None,
            TppPreset::OralSmallMolecule,
        )?;
        assert!(ev.tpp_mpo.is_some());
        assert!(!ev.metabolism_sites.is_empty());
        return Ok(());
    }

    let pack =
        load_default_model_pack().ok_or("expected default TDC model pack to load cleanly")?;
    assert!(
        pack.models.len() >= 18,
        "expected >= 18 TDC QSAR endpoints, got {}",
        pack.models.len()
    );

    // Verify held-out Bemis-Murcko scaffold-split acceptance criteria
    let mut cls_aucs = Vec::new();
    let mut covs = Vec::new();
    for m in &pack.models {
        covs.push(m.metrics.conformal_90_coverage);
        if m.task == QsarTaskKind::BinaryClassification {
            cls_aucs.push(m.metrics.roc_auc.unwrap_or(0.5));
        }
    }
    let mean_auc = cls_aucs.iter().sum::<f64>() / (cls_aucs.len().max(1) as f64);
    let mean_cov = covs.iter().sum::<f64>() / (covs.len().max(1) as f64);
    assert!(
        mean_auc >= 0.80,
        "expected mean scaffold ROC AUC >= 0.80 across TDC classification endpoints, got {mean_auc:.3}"
    );
    assert!(
        (0.85..=0.96).contains(&mean_cov),
        "expected mean 90% split-conformal coverage in [0.85, 0.96], got {mean_cov:.3}"
    );

    // Known clinical hERG/CYP3A4/P-gp reference: Verapamil
    let ver = evaluate_smiles_full(
        "COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC",
        Some("Verapamil"),
        Some(pack),
        TppPreset::OralSmallMolecule,
    )?;
    let vq = ver.qsar.as_ref().ok_or("missing qsar on Verapamil")?;
    assert_eq!(vq.overall_ad.status, ApplicabilityDomainStatus::InDomain);
    let herg = vq
        .endpoints
        .iter()
        .find(|e| e.endpoint_id == "herg_karim")
        .ok_or("missing herg_karim")?;
    assert!(
        herg.value >= 0.70,
        "Verapamil must predict high hERG liability, got {:.3}",
        herg.value
    );
    let ppbr = vq
        .endpoints
        .iter()
        .find(|e| e.endpoint_id == "ppbr_az")
        .ok_or("missing ppbr_az")?;
    assert!(
        (0.0..=100.0).contains(&ppbr.value)
            && (0.0..=100.0).contains(&ppbr.ci90_low)
            && (0.0..=100.0).contains(&ppbr.ci90_high),
        "PPBR prediction & 90% CI must stay bounded in [0, 100]%"
    );
    let vtpp = ver.tpp_mpo.as_ref().ok_or("missing tpp_mpo")?;
    assert!(
        !vtpp.hard_stops.is_empty(),
        "Verapamil should trigger hERG hard stop in Oral TPP"
    );

    // Clean oral CNS drug reference: Caffeine
    let caf = evaluate_smiles_full(
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        Some("Caffeine"),
        Some(pack),
        TppPreset::CnsPenetrant,
    )?;
    let cq = caf.qsar.as_ref().ok_or("missing qsar on Caffeine")?;
    let caf_herg = cq
        .endpoints
        .iter()
        .find(|e| e.endpoint_id == "herg_karim")
        .ok_or("missing herg_karim")?;
    assert!(
        caf_herg.value <= 0.45,
        "Caffeine should predict low hERG risk, got {:.3}",
        caf_herg.value
    );
    let caf_bbb = cq
        .endpoints
        .iter()
        .find(|e| e.endpoint_id == "bbb_martins")
        .ok_or("missing bbb_martins")?;
    assert!(
        caf_bbb.value >= 0.70,
        "Caffeine should predict high BBB penetration, got {:.3}",
        caf_bbb.value
    );

    // Ensure reload from disk works identically
    let reloaded = AdmetModelPack::load_from_path(&pack_path)?;
    assert_eq!(reloaded.models.len(), pack.models.len());
    Ok(())
}
