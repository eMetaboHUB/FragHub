// src/update_checker.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use csv::WriterBuilder;

/// Vérifie quels spectres ont déjà été traités lors d'une session précédente.
///
/// Pour un développeur Python : Un "Set" (ensemble) en Python utilise des hash tables. En Rust, 
/// on utilise un `HashSet`. La différence est que la recherche `.contains()` est garantie en `O(1)`
/// de manière très stricte et déterministe au niveau de la RAM, là où Python ajoute une couche d'abstraction ("overhead").
pub fn check_for_update_processing(
    py: Python,
    spectrum_list: Vec<Spectrum>,
    output_directory: String,
    ordered_columns: Vec<String>,
    append_mode: bool,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<(Vec<Spectrum>, bool, usize)> {

    let total_items = spectrum_list.len();

    // ⚠️ ORDRE CRITIQUE POUR L'INTERFACE VUE.JS
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Checking for updates:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    // 1. Lire le fichier updates.json existant
    let update_file_path = Path::new(&output_directory).join("updates.json");
    let mut splash_set: HashSet<String> = HashSet::new();
    let mut json_data = serde_json::json!({"SPLASH_LIST": {}});

    if update_file_path.exists() {
        if let Ok(file) = fs::File::open(&update_file_path) {
            if let Ok(parsed) = serde_json::from_reader::<_, serde_json::Value>(file) {
                json_data = parsed;
                if let Some(splash_list) = json_data.get("SPLASH_LIST").and_then(|v| v.as_object()) {
                    for key in splash_list.keys() {
                        splash_set.insert(key.clone()); // On stocke les SPLASH connus dans un Hash ultra-rapide
                    }
                }
            }
        }
    }

    // --- CORRECTION : Utilisation d'un HashSet pour une recherche en O(1) ---
    let mut indices_to_keep: HashSet<usize> = HashSet::new();
    let mut indices_to_delete = Vec::new();
    let mut new_splashes = Vec::new();
    let mut processed = 0;

    // 2. Parcourir les spectres
    for (i, spec) in spectrum_list.iter().enumerate() {
        let splash = spec.metadata.get("SPLASH").cloned().unwrap_or_default();

        if !splash.is_empty() && splash_set.contains(&splash) {
            indices_to_delete.push(i); // Déjà vu -> on supprime
        } else {
            // --- CORRECTION : Utilisation de .insert() au lieu de .push() ---
            indices_to_keep.insert(i); // Nouveau -> on garde
            if !splash.is_empty() {
                new_splashes.push(splash);
            }
        }

        processed += 1;
        if processed % 1000 == 0 {
            if let Some(cb) = &progress_callback { cb.call1(py, (processed,))?; }
        }
    }

    // ⚠️ GARANTIE DU 100% DE L'ÉTAPE 2
    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Saving updates and logs...",))?; }

    let (final_list, update, deleted_count) = py.allow_threads(|| {
        // 3. Écrire les doublons supprimés dans le CSV
        if !indices_to_delete.is_empty() {
            let deleted_dir = Path::new(&output_directory).join("DELETED_SPECTRUMS");
            fs::create_dir_all(&deleted_dir)?;
            let file_path = deleted_dir.join("previously_cleaned.csv");

            // On preserve l'historique des runs precedents (mode append, sans reset_updates) :
            // si le fichier existe deja, on ajoute les lignes a la suite au lieu de l'ecraser.
            let is_append = append_mode && file_path.exists();
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
                record.push("spectrum deleted because already processed in a previous run.".to_string());
                wtr.write_record(&record).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            }
            wtr.flush()?;
        }

        // 4. Mettre à jour le fichier JSON
        let update = !new_splashes.is_empty();
        if update {
            let file = std::fs::File::create(&update_file_path).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            let mut wtr = std::io::BufWriter::new(file);
            use std::io::Write;
            
            wtr.write_all(b"{\"SPLASH_LIST\":{").map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            
            let mut first = true;
            for s in &splash_set {
                if !first { wtr.write_all(b",").map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?; }
                serde_json::to_writer(&mut wtr, s).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
                wtr.write_all(b":true").map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
                first = false;
            }
            for s in &new_splashes {
                if !first { wtr.write_all(b",").map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?; }
                serde_json::to_writer(&mut wtr, s).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
                wtr.write_all(b":true").map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
                first = false;
            }
            
            wtr.write_all(b"}}").map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            wtr.flush().map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }

        // 5. Générer la liste finale
        let mut final_list = Vec::with_capacity(indices_to_keep.len());
        let mut current_idx = 0;
        for spec in spectrum_list.into_iter() {
            if indices_to_keep.contains(&current_idx) {
                final_list.push(spec);
            }
            current_idx += 1;
        }

        Ok::<_, pyo3::PyErr>((final_list, update, indices_to_delete.len()))
    })?;

    // On renvoie la nouvelle liste, le booléen (y a-t-il eu une MAJ ?), et le nombre de supprimés
    Ok((final_list, update, deleted_count))
}