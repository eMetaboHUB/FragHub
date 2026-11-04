// src/convertors/mgf_to_dict.rs

use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use rayon::prelude::*;
use std::collections::HashMap;
use crate::convertors::raw_metadata::RawMetadata;

// 1. La structure pour passer les données entre les cœurs CPU sans bloquer Python
/// Structure interne pour le multithreading.
///
/// Pour un développeur Python : En Python, les objets créés dans différents processus
/// doivent être sérialisés (Pickle) pour être renvoyés au processus principal, ce qui est très lent.
/// En Rust, on passe cette structure entre les threads (cœurs CPU) de manière ultra-rapide
/// et sans copie superflue de la mémoire.
struct ParsedSpectrum {
    metadata: RawMetadata, // champs bruts, dans l'ordre du fichier
    peaks: Vec<Vec<f64>>,
}

/// Fonction principale de parsing MGF via Rayon.
///
/// # Arguments
/// * `py` (Python) : Le token PyO3.
/// * `final_mgf` (Vec<String>) : Les spectres MGF bruts.
/// * `keys_dict` (HashMap<String, String>) : Mapping des clés.
/// * `keys_list` (Vec<String>) : Liste des clés à conserver.
/// * `progress_callback` (Option<PyObject>) : Callback de progression.
/// * `total_items_callback` (Option<PyObject>) : Callback de total.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type.
///
/// # Returns
/// * `PyResult<Vec<Spectrum>>` : La liste de spectres parsés.
pub fn mgf_to_dict_processing(
    py: Python,
    final_mgf: Vec<String>,
    keys_dict: HashMap<String, String>,
    keys_list: Vec<String>,
    db_name: String,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {
    let total = final_mgf.len();
    if let Some(cb) = &total_items_callback { cb.call1(py, (total, 0))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Parsing MGF spectrums:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    let mut result_list = Vec::new();
    let mut processed = 0;

    // Découpage par paquets (chunks) comme dans votre version Python
    let chunk_size = 2000;

    for chunk in final_mgf.chunks(chunk_size) {

        // ====================================================================
        // ÉTAPE 1 : MULTITHREADING PUR RUST (Équivalent de ThreadPoolExecutor)
        // ====================================================================
        let parsed_chunk: Vec<ParsedSpectrum> = chunk.par_iter().filter_map(|spectrum| {
            // Attention : Assurez-vous d'avoir bien retiré le "?" dans globals_vars.rs
            // pour METADATA_PEAK_LIST_SPLIT_PATTERN_MGF comme discuté précédemment !
            if let Some(caps) = crate::globals_vars::METADATA_PEAK_LIST_SPLIT_PATTERN_MGF.captures(spectrum) {
                let metadata_str = caps.get(1).map_or("", |m| m.as_str());
                let peaks_str = caps.get(2).map_or("", |m| m.as_str());

                if !metadata_str.is_empty() && !peaks_str.is_empty() {
                    let mut metadata_dict = RawMetadata::default();

                    // Extraction des métadonnées
                    for m in crate::globals_vars::METADATA_PATTERN_MGF.captures_iter(metadata_str) {
                        let k_raw = m.get(1).unwrap().as_str();
                        let v_raw = m.get(2).unwrap().as_str();

                        let k = crate::globals_vars::METADATA_FIELDS_NAME_PATTERN.replace_all(k_raw, "").to_lowercase().trim().to_string();
                        let v = crate::globals_vars::METADATA_STRIP_VALUE_PATTERN.replace_all(v_raw, "").to_string();
                        metadata_dict.insert(k, v);
                    }

                    if !metadata_dict.is_empty() {
                        let mut peaks_array = Vec::new();

                        // Extraction des pics
                        for p in crate::globals_vars::PEAK_LIST_SPLIT_PATTERN.captures_iter(peaks_str) {
                            let mz = p.get(1).unwrap().as_str().replace(",", ".").parse::<f64>().unwrap_or(0.0);
                            let int = p.get(2).unwrap().as_str().replace(",", ".").parse::<f64>().unwrap_or(0.0);
                            peaks_array.push(vec![mz, int]);
                        }

                        // Si on a des métadonnées ET des pics, on valide le spectre
                        if !peaks_array.is_empty() {
                            return Some(ParsedSpectrum { metadata: metadata_dict, peaks: peaks_array });
                        }
                    }
                }
            }
            None // Équivalent de "if not metadata or not peak_list: return None"
        }).collect();

        // ====================================================================
        // ÉTAPE 2 : Mapping et création du dictionnaire (RUST NATIVE)
        // ====================================================================
        for parsed in parsed_chunk {
            let mut spec = Spectrum::default();

            // Renommage vers les clés FragHub : si deux champs donnent la même clé, le premier du fichier gagne
            parsed.metadata.map_into(&keys_dict, |mapped| keys_list.iter().any(|k| k == mapped), &mut spec.metadata);
            spec.metadata.insert("DATABASE_NAME".to_string(), db_name.clone());

            // Ajout de la liste de pics
            spec.peaks = parsed.peaks.into_iter().map(|p| (p[0], p[1])).collect();

            // Complétion des clés manquantes
            for key in &keys_list {
                if !spec.metadata.contains_key(key) && key != "PEAKS_LIST" { 
                    spec.metadata.insert(key.clone(), "".to_string()); 
                }
            }

            result_list.push(spec);
        }

        // ====================================================================
        // ÉTAPE 3 : MISE À JOUR DE LA BARRE DE PROGRESSION VUE.JS
        // ====================================================================
        processed += chunk.len();
        if let Some(cb) = &progress_callback { cb.call1(py, (processed,))?; }
    }

    Ok(result_list)
}