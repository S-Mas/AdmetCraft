#![forbid(unsafe_code)]

use std::path::PathBuf;

/// Returns the centralized Craft bio-data directory (`CRAFT_DATA_DIR`, default
/// `/home/acs/Datos/CraftBio` when `/home/acs/Datos` exists, else `$XDG_DATA_HOME/CraftBio`).
#[must_use]
pub fn craft_data_dir() -> PathBuf {
    if let Ok(raw) = std::env::var("CRAFT_DATA_DIR") {
        let trimmed = raw.trim();
        if !trimmed.is_empty() {
            return PathBuf::from(trimmed);
        }
    }
    let datos = PathBuf::from("/home/acs/Datos");
    if datos.is_dir() {
        return datos.join("CraftBio");
    }
    if let Ok(xdg) = std::env::var("XDG_DATA_HOME") {
        let trimmed = xdg.trim();
        if !trimmed.is_empty() {
            return PathBuf::from(trimmed).join("CraftBio");
        }
    }
    if let Ok(home) = std::env::var("HOME") {
        return PathBuf::from(home).join(".local/share/CraftBio");
    }
    PathBuf::from("CraftBio")
}

/// Directory storing trained AdmetCraft v1.0 QSAR model packs (`.../models/admetcraft_v1`).
#[must_use]
pub fn admet_models_dir() -> PathBuf {
    craft_data_dir().join("models").join("admetcraft_v1")
}

/// Default trained TDC QSAR model pack path (`.../models/admetcraft_v1/admetcraft_tdc_v1.json`).
#[must_use]
pub fn default_model_pack_path() -> PathBuf {
    admet_models_dir().join("admetcraft_tdc_v1.json")
}

/// Optional user-enriched SEA off-target index JSON (`.../models/admetcraft_v1/sea_offtargets.json`).
#[must_use]
pub fn default_sea_index_path() -> PathBuf {
    admet_models_dir().join("sea_offtargets.json")
}

/// Directory storing raw downloaded TDC `.tab` datasets (`.../models/tdc_raw`).
#[must_use]
pub fn tdc_raw_dir() -> PathBuf {
    craft_data_dir().join("models").join("tdc_raw")
}
