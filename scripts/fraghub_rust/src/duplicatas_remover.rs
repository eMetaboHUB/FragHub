// src/duplicatas_remover.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use csv::WriterBuilder;

/// Supprime les doublons (spectres identiques) en gardant celui qui contient le plus de métadonnées.
///
/// Pour un développeur Python : Retirer des doublons est très long si c'est mal optimisé (`drop_duplicates` etc.).
/// En Rust, on fait un seul passage (`O(N)`) sur les spectres. On calcule la taille (en caractères)
/// des métadonnées de chaque spectre, et on stocke l'index du "gagnant" (le plus gros) dans
/// une table de hachage (`best_indices`) utilisant un Tuple `(SPLASH, INCHIKEY)` comme clé composite.
pub fn remove_duplicatas_processing(
    py: Python,
    spectrum_list: Vec<Spectrum>,
    output_directory: String,
    ordered_columns: Vec<String>,
    update: bool,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<(Vec<Spectrum>, usize)> {

    let total_items = spectrum_list.len();

    // ⚠️ ORDRE CRITIQUE POUR L'INTERFACE VUE.JS
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Removing duplicates:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    let mut best_indices: HashMap<(String, String), (usize, usize)> = HashMap::new();
    let mut empty_inchi_indices: Vec<usize> = Vec::new();

    let mut processed = 0;

    for (i, spec) in spectrum_list.iter().enumerate() {
        let mut row_size = 0;
        for col in &ordered_columns {
            if let Some(val) = spec.metadata.get(col) {
                row_size += val.chars().count();
            }
        }

        let splash = spec.metadata.get("SPLASH").cloned().unwrap_or_default();
        let inchikey = spec.metadata.get("INCHIKEY").map(|s| s.trim().to_string()).unwrap_or_default();

        if inchikey.is_empty() || inchikey.to_lowercase() == "nan" || inchikey.to_lowercase() == "none" {
            empty_inchi_indices.push(i);
        } else {
            let key = (splash, inchikey);
            if let Some(&(_, best_size)) = best_indices.get(&key) {
                if row_size > best_size {
                    best_indices.insert(key, (i, row_size));
                }
            } else {
                best_indices.insert(key, (i, row_size));
            }
        }

        processed += 1;
        if processed % 1000 == 0 {
            if let Some(cb) = &progress_callback { cb.call1(py, (processed,))?; }
        }
    }

    let mut indices_to_keep: HashSet<usize> = HashSet::new();
    for &idx in &empty_inchi_indices { indices_to_keep.insert(idx); }
    for &(idx, _) in best_indices.values() { indices_to_keep.insert(idx); }

    let mut indices_to_delete: Vec<usize> = Vec::new();
    for i in 0..total_items {
        if !indices_to_keep.contains(&i) {
            indices_to_delete.push(i);
        }
    }

    if !indices_to_delete.is_empty() {
        let deleted_dir = Path::new(&output_directory).join("DELETED_SPECTRUMS");
        fs::create_dir_all(&deleted_dir)?;
        let file_path = deleted_dir.join("duplicatas_removed.csv");

        // On preserve l'historique des runs precedents (mode append, sans reset_updates) :
        // si le fichier existe deja, on ajoute les lignes a la suite au lieu de l'ecraser.
        let is_append = update && file_path.exists();
        let file = std::fs::OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(&file_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

        let mut wtr = WriterBuilder::new()
            .delimiter(b'\t')
            .quote(b'"')
            .has_headers(!is_append)
            .from_writer(file);

        if !is_append {
            let mut header = ordered_columns.clone();
            header.push("DELETION_REASON".to_string());
            wtr.write_record(&header).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }

        for &idx in &indices_to_delete {
            let spec = &spectrum_list[idx];
            let mut record: Vec<String> = Vec::with_capacity(ordered_columns.len() + 1);

            for col in &ordered_columns {
                record.push(spec.metadata.get(col).cloned().unwrap_or_default());
            }
            record.push("spectrum deleted because it's a duplicate (SPLASH + INCHIKEY)".to_string());
            wtr.write_record(&record).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        wtr.flush()?;
    }

    let mut final_list = Vec::with_capacity(indices_to_keep.len());
    let mut current_idx = 0;
    for spec in spectrum_list.into_iter() {
        if indices_to_keep.contains(&current_idx) {
            final_list.push(spec);
        }
        current_idx += 1;
    }

    // ⚠️ GARANTIE DU 100% POUR CLÔTURER L'INTERFACE PROPREMENT
    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }

    Ok((final_list, indices_to_delete.len()))
}