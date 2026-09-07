// src/writers.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::fs::{self, OpenOptions};
use std::io::{Write, Seek, SeekFrom, Read};
use std::path::Path;
use csv::WriterBuilder;
use once_cell::sync::Lazy;
use regex::Regex;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

static RE_3: Lazy<Regex> = Lazy::new(|| Regex::new(r#"\[\n\s*(-?[\d\.eE\+\-]+),\n\s*(-?[\d\.eE\+\-]+),\n\s*"(.*?)"\n\s*\]"#).unwrap());
static RE_2: Lazy<Regex> = Lazy::new(|| Regex::new(r#"\[\n\s*(-?[\d\.eE\+\-]+),\n\s*(-?[\d\.eE\+\-]+)\n\s*\]"#).unwrap());

/// Regroupe les fonctions d'écriture hautement optimisées pour générer les fichiers MSP.
///
/// Pour un développeur Python : Au lieu de garder des gigaoctets de texte en RAM
/// ou de faire des écritures uniques coûteuses, on écrit par "flux" (streams) binaires.
/// L'objet `AtomicUsize` (ex: `global_processed`) permet à de multiples threads concurrents
/// de mettre à jour la barre de progression simultanément, sans aucune corruption mémoire.
///
/// # Arguments
/// * `py` (Python) : Token d'accès au GIL.
/// * `pos_lc` (Vec<String>), `pos_lc_insilico`... : Les différentes listes de spectres formatés selon la polarité/chromato.
/// * `output_directory` (&str) : Dossier cible.
/// * `update` (bool) : Indicateur pour savoir s'il faut écraser ou ajouter (append) au fichier existant.
/// * `progress_callback`, `total_items_callback`, `prefix_callback`, `item_type_callback` : Callbacks UI.
///
/// # Returns
/// * `PyResult<()>` : Succès ou erreur IO.
#[allow(clippy::too_many_arguments)]
pub fn writting_msp_processing(
    py: Python, pos_lc: Vec<String>, pos_lc_insilico: Vec<String>, pos_gc: Vec<String>, pos_gc_insilico: Vec<String>, neg_lc: Vec<String>, neg_lc_insilico: Vec<String>, neg_gc: Vec<String>, neg_gc_insilico: Vec<String>, output_directory: &str, update: bool, progress_callback: Option<PyObject>, total_items_callback: Option<PyObject>, prefix_callback: Option<PyObject>, item_type_callback: Option<PyObject>,
) -> PyResult<()> {

    let tasks = vec![
        (&pos_lc, "POS_LC.msp", "POS"),
        (&pos_lc_insilico, "POS_LC_insilico.msp", "POS"),
        (&pos_gc, "POS_GC.msp", "POS"),
        (&pos_gc_insilico, "POS_GC_insilico.msp", "POS"),
        (&neg_lc, "NEG_LC.msp", "NEG"),
        (&neg_lc_insilico, "NEG_LC_insilico.msp", "NEG"),
        (&neg_gc, "NEG_GC.msp", "NEG"),
        (&neg_gc_insilico, "NEG_GC_insilico.msp", "NEG"),
    ];

    let total_items: usize = tasks.iter().map(|(list, _, _)| list.len()).sum();
    if total_items == 0 { return Ok(()); }

    if let Some(cb) = &prefix_callback { cb.call1(py, ("Writing MSP files:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }

    let global_processed = Arc::new(AtomicUsize::new(0));
    let progress_cb = progress_callback.clone();

    let results: Vec<Result<(), String>> = py.allow_threads(|| {
        tasks.par_iter().map(|(data_list, filename, mode)| {
            if data_list.is_empty() { return Ok(()); }

            let path_dir = Path::new(output_directory).join("MSP").join(mode);
            fs::create_dir_all(&path_dir).map_err(|e| e.to_string())?;
            let file_path = path_dir.join(filename);

            let mut file = OpenOptions::new().write(true).create(true).append(update).truncate(!update).open(&file_path).map_err(|e| e.to_string())?;

            let mut local_processed = 0;
            for spec_str in *data_list {
                file.write_all(spec_str.as_bytes()).map_err(|e| e.to_string())?;
                file.write_all(b"\n\n").map_err(|e| e.to_string())?;

                local_processed += 1;
                if local_processed % 5000 == 0 {
                    let current = global_processed.fetch_add(5000, Ordering::Relaxed) + 5000;
                    if let Some(ref cb) = progress_cb {
                        Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                    }
                }
            }

            let remainder = local_processed % 5000;
            if remainder > 0 {
                let current = global_processed.fetch_add(remainder, Ordering::Relaxed) + remainder;
                if let Some(ref cb) = progress_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                }
            }
            file.flush().map_err(|e| e.to_string())?;
            Ok(())
        }).collect()
    });

    for res in results {
        if let Err(e) = res { return Err(pyo3::exceptions::PyIOError::new_err(e)); }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    Ok(())
}

/// Regroupe les fonctions d'écriture hautement optimisées pour générer les fichiers CSV.
///
/// # Arguments
/// * `py` (Python) : Token d'accès au GIL.
/// * `pos_lc_df` (Vec<Spectrum>), `pos_gc_df`... : Les différentes listes de spectres formatés.
/// * `ordered_columns` (Vec<String>) : Liste ordonnée des colonnes du CSV.
/// * `output_directory` (&str) : Dossier cible.
/// * `update` (bool) : Indicateur pour savoir s'il faut écraser ou ajouter (append) au fichier existant.
/// * `progress_callback`, `total_items_callback`, `prefix_callback`, `item_type_callback` : Callbacks UI.
///
/// # Returns
/// * `PyResult<()>` : Succès ou erreur IO.
#[allow(clippy::too_many_arguments)]
pub fn writting_csv_processing(
    py: Python, pos_lc_df: Vec<Spectrum>, pos_gc_df: Vec<Spectrum>, neg_lc_df: Vec<Spectrum>, neg_gc_df: Vec<Spectrum>, pos_lc_df_insilico: Vec<Spectrum>, pos_gc_df_insilico: Vec<Spectrum>, neg_lc_df_insilico: Vec<Spectrum>, neg_gc_df_insilico: Vec<Spectrum>, ordered_columns: Vec<String>, output_directory: &str, update: bool, progress_callback: Option<PyObject>, total_items_callback: Option<PyObject>, prefix_callback: Option<PyObject>, item_type_callback: Option<PyObject>,
) -> PyResult<()> {

    let tasks = vec![
        (&pos_lc_df, "POS_LC.csv", "POS"),
        (&pos_gc_df, "POS_GC.csv", "POS"),
        (&neg_lc_df, "NEG_LC.csv", "NEG"),
        (&neg_gc_df, "NEG_GC.csv", "NEG"),
        (&pos_lc_df_insilico, "POS_LC_In_Silico.csv", "POS"),
        (&pos_gc_df_insilico, "POS_GC_In_Silico.csv", "POS"),
        (&neg_lc_df_insilico, "NEG_LC_In_Silico.csv", "NEG"),
        (&neg_gc_df_insilico, "NEG_GC_In_Silico.csv", "NEG"),
    ];

    let total_items: usize = tasks.iter().map(|(list, _, _)| list.len()).sum();
    if total_items == 0 { return Ok(()); }

    if let Some(cb) = &prefix_callback { cb.call1(py, ("Writing CSV files:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("rows",))?; }
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }

    let global_processed = Arc::new(AtomicUsize::new(0));
    let progress_cb = progress_callback.clone();
    let prefix_cb = prefix_callback.clone();

    let results: Vec<Result<(), String>> = py.allow_threads(|| {
        tasks.par_iter().map(|(data_list, filename, mode)| {
            if data_list.is_empty() { return Ok(()); }

            let path_dir = Path::new(output_directory).join("CSV").join(mode);
            fs::create_dir_all(&path_dir).map_err(|e| e.to_string())?;
            let file_path = path_dir.join(filename);
            let is_append = update && file_path.exists();

            let file = OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(&file_path).map_err(|e| e.to_string())?;
            let mut wtr = WriterBuilder::new().delimiter(b';').quote(b'"').has_headers(!is_append).from_writer(file);

            if !is_append {
                wtr.write_record(&ordered_columns).map_err(|e| e.to_string())?;
            }

            let mut local_processed = 0;
            for spec in *data_list {
                let mut row = Vec::with_capacity(ordered_columns.len());
                for col in &ordered_columns {
                    let mut cell_val = String::new();
                    if col == "PEAKS_LIST" {
                        // --- NOUVEAU : Récupération des annotations De Novo ---
                        let mut used_metadata = false;
                        if let Some(val) = spec.metadata.get("PEAKS_LIST") {
                            if !val.trim().is_empty() && val != "NOT FOUND" {
                                let sep = if val.contains(';') { ';' } else { '\n' };
                                let lines_count = val.split(sep).filter(|s| !s.trim().is_empty()).count();

                                if lines_count == spec.peaks.len() && val.chars().any(|c| c.is_ascii_alphabetic()) {
                                    cell_val = val.replace("\n", ";"); // Force le format CSV
                                    used_metadata = true;
                                }
                            }
                        }

                        // Fallback rapide
                        if !used_metadata && !spec.peaks.is_empty() {
                            let mut peaks_str = String::with_capacity(spec.peaks.len() * 20);
                            for (i, &(mz, int)) in spec.peaks.iter().enumerate() {
                                if i > 0 { peaks_str.push(';'); }
                                peaks_str.push_str(&format!("{} {}", mz, int));
                            }
                            cell_val = peaks_str;
                        }
                    } else if let Some(val) = spec.metadata.get(col) {
                        if !val.eq_ignore_ascii_case("nan") {
                            cell_val = val.clone();
                        }
                    }
                    row.push(cell_val);
                }
                wtr.write_record(&row).map_err(|e| e.to_string())?;

                local_processed += 1;
                if local_processed % 1000 == 0 {
                    let current = global_processed.fetch_add(1000, Ordering::Relaxed) + 1000;
                    if let Some(ref cb) = progress_cb {
                        Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                    }
                }
            }

            let remainder = local_processed % 1000;
            if remainder > 0 {
                let current = global_processed.fetch_add(remainder, Ordering::Relaxed) + remainder;
                if let Some(ref cb) = progress_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                }
            }

            if let Some(ref cb) = prefix_cb {
                Python::with_gil(|py| { let _ = cb.call1(py, (format!("Writing {} to disk...", filename),)); });
            }

            wtr.flush().map_err(|e| e.to_string())?;
            Ok(())
        }).collect()
    });

    for res in results {
        if let Err(e) = res { return Err(pyo3::exceptions::PyIOError::new_err(e)); }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    Ok(())
}

/// Regroupe les fonctions d'écriture hautement optimisées pour générer les fichiers JSON.
///
/// # Arguments
/// * `py` (Python) : Token d'accès au GIL.
/// * `update` (bool) : Indicateur pour savoir s'il faut écraser ou ajouter (append) au fichier existant.
/// * `pos_lc_df` (Vec<Spectrum>), `pos_gc_df`... : Les différentes listes de spectres formatés.
/// * `output_directory` (&str) : Dossier cible.
/// * `progress_callback`, `total_items_callback`, `prefix_callback`, `item_type_callback` : Callbacks UI.
///
/// # Returns
/// * `PyResult<()>` : Succès ou erreur IO.
#[allow(clippy::too_many_arguments)]
pub fn writting_json_processing(
    py: Python, update: bool, pos_lc_df: Vec<Spectrum>, pos_gc_df: Vec<Spectrum>, neg_lc_df: Vec<Spectrum>, neg_gc_df: Vec<Spectrum>, pos_lc_df_insilico: Vec<Spectrum>, pos_gc_df_insilico: Vec<Spectrum>, neg_lc_df_insilico: Vec<Spectrum>, neg_gc_df_insilico: Vec<Spectrum>, ordered_columns: Vec<String>, output_directory: &str, progress_callback: Option<PyObject>, total_items_callback: Option<PyObject>, prefix_callback: Option<PyObject>, item_type_callback: Option<PyObject>,
) -> PyResult<()> {

    let tasks = vec![
        (&pos_lc_df, "POS_LC.json", "POS"),
        (&pos_gc_df, "POS_GC.json", "POS"),
        (&neg_lc_df, "NEG_LC.json", "NEG"),
        (&neg_gc_df, "NEG_GC.json", "NEG"),
        (&pos_lc_df_insilico, "POS_LC_In_Silico.json", "POS"),
        (&pos_gc_df_insilico, "POS_GC_In_Silico.json", "POS"),
        (&neg_lc_df_insilico, "NEG_LC_In_Silico.json", "NEG"),
        (&neg_gc_df_insilico, "NEG_GC_In_Silico.json", "NEG"),
    ];

    let total_items: usize = tasks.iter().map(|(list, _, _)| list.len()).sum();
    if total_items == 0 { return Ok(()); }

    if let Some(cb) = &prefix_callback { cb.call1(py, ("Writing JSON files:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("rows",))?; }
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }

    let global_processed = Arc::new(AtomicUsize::new(0));
    let progress_cb = progress_callback.clone();
    let prefix_cb = prefix_callback.clone();

    let results: Vec<Result<(), String>> = py.allow_threads(|| {
        tasks.par_iter().map(|(data_list, filename, mode)| {
            if data_list.is_empty() { return Ok(()); }

            let path_dir = Path::new(output_directory).join("JSON").join(mode);
            fs::create_dir_all(&path_dir).map_err(|e| e.to_string())?;
            let file_path = path_dir.join(filename);

            let is_append_mode = update && file_path.exists() && fs::metadata(&file_path).map(|m| m.len()).unwrap_or(0) > 2;

            let file_res = if is_append_mode {
                OpenOptions::new().read(true).write(true).open(&file_path).and_then(|mut f| {
                    let file_len = f.metadata()?.len();
                    if file_len > 0 {
                        let mut buf = [0u8; 1];
                        for offset in 1..=std::cmp::min(10, file_len) {
                            f.seek(SeekFrom::End(-(offset as i64)))?;
                            f.read_exact(&mut buf)?;
                            if buf[0] == b']' {
                                f.set_len(file_len - offset)?;
                                f.seek(SeekFrom::End(0))?;
                                f.write_all(b",\n")?;
                                break;
                            }
                        }
                    }
                    Ok(f)
                })
            } else {
                OpenOptions::new().write(true).create(true).truncate(true).open(&file_path).and_then(|mut f| {
                    f.write_all(b"[\n")?;
                    Ok(f)
                })
            };

            let mut file = file_res.map_err(|e| e.to_string())?;

            let columns: Vec<String> = ordered_columns.clone();

            let len = data_list.len();
            let mut local_processed = 0;

            for (i, spec) in data_list.iter().enumerate() {
                let mut map = serde_json::Map::new();

                for col in &columns {
                    if col == "PEAKS_LIST" || col == "NUM PEAKS" || col.starts_with("TREE_") { continue; }

                    let mut val_str = "NaN".to_string();
                    if let Some(val) = spec.metadata.get(col) {
                        let clean_val = val.replace("\"", "");
                        if !clean_val.is_empty() && !clean_val.eq_ignore_ascii_case("nan") { val_str = clean_val; }
                    }

                    if col == "MSLEVEL" { if let Ok(num) = val_str.parse::<i64>() { map.insert(col.clone(), serde_json::json!(num)); continue; } }
                    if ["PRECURSORMZ", "RT", "ENTROPY"].contains(&col.as_str()) {
                        if let Ok(num) = val_str.parse::<f64>() { map.insert(col.clone(), serde_json::json!(num)); continue; }
                    }
                    map.insert(col.clone(), serde_json::Value::String(val_str));
                }

                let num_peaks_str = spec.metadata.get("NUM PEAKS").cloned().unwrap_or_else(|| "0".to_string());
                map.insert("NUM PEAKS".to_string(), serde_json::json!(num_peaks_str.parse::<i64>().unwrap_or(0)));

                let mut peaks_array = Vec::new();

                // --- NOUVEAU : Récupération intelligente et parsing pour le format JSON ---
                if let Some(val) = spec.metadata.get("PEAKS_LIST") {
                    if !val.trim().is_empty() && val != "NOT FOUND" {
                        let sep = if val.contains(';') { ';' } else { '\n' };
                        let lines_count = val.split(sep).filter(|s| !s.trim().is_empty()).count();

                        if lines_count == spec.peaks.len() && val.chars().any(|c| c.is_ascii_alphabetic()) {
                            for line in val.split(sep).filter(|s| !s.trim().is_empty()) {
                                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                                if parts.len() >= 3 {
                                    let mz: f64 = parts[0].parse().unwrap_or(0.0);
                                    let int: f64 = parts[1].parse().unwrap_or(0.0);
                                    let formula = parts[2..].join(" ");
                                    peaks_array.push(serde_json::json!([mz, int, formula]));
                                } else if parts.len() == 2 {
                                    let mz: f64 = parts[0].parse().unwrap_or(0.0);
                                    let int: f64 = parts[1].parse().unwrap_or(0.0);
                                    peaks_array.push(serde_json::json!([mz, int]));
                                }
                            }
                        }
                    }
                }

                // Si aucune donnée de métadonnée valide, on utilise le tableau brute classique
                if peaks_array.is_empty() {
                    for &(mz, intensity) in &spec.peaks {
                        peaks_array.push(serde_json::json!([mz, intensity]));
                    }
                }

                map.insert("PEAKS_LIST".to_string(), serde_json::Value::Array(peaks_array));

                let item_str_pretty = serde_json::to_string_pretty(&map).unwrap();
                // RE_3 compressera parfaitement tes JSON à 3 arguments !
                let compacted_1 = RE_3.replace_all(&item_str_pretty, "[$1, $2, \"$3\"]").to_string();
                let compacted_2 = RE_2.replace_all(&compacted_1, "[$1, $2]").to_string();
                let indented_str = format!("  {}", compacted_2.replace('\n', "\n  "));

                file.write_all(indented_str.as_bytes()).map_err(|e| e.to_string())?;
                if i < len - 1 {
                    file.write_all(b",\n").map_err(|e| e.to_string())?;
                } else {
                    file.write_all(b"\n").map_err(|e| e.to_string())?;
                }

                local_processed += 1;
                if local_processed % 500 == 0 {
                    let current = global_processed.fetch_add(500, Ordering::Relaxed) + 500;
                    if let Some(ref cb) = progress_cb {
                        Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                    }
                }
            }

            let remainder = local_processed % 500;
            if remainder > 0 {
                let current = global_processed.fetch_add(remainder, Ordering::Relaxed) + remainder;
                if let Some(ref cb) = progress_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                }
            }

            if let Some(ref cb) = prefix_cb {
                Python::with_gil(|py| { let _ = cb.call1(py, (format!("Writing {} to disk...", filename),)); });
            }

            file.write_all(b"]").map_err(|e| e.to_string())?;
            file.flush().map_err(|e| e.to_string())?;
            Ok(())
        }).collect()
    });

    for res in results {
        if let Err(e) = res { return Err(pyo3::exceptions::PyIOError::new_err(e)); }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    Ok(())
}
