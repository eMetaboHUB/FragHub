// src/writers/csv_writer.rs
//! Écriture des fichiers CSV (séparateur « ; »).
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use ::csv::WriterBuilder;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::{WRITE_CHUNK, WRITE_CHUNKS_PER_ROUND};

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
    py: Python, pos_lc_df: &[Spectrum], pos_gc_df: &[Spectrum], neg_lc_df: &[Spectrum], neg_gc_df: &[Spectrum], pos_lc_df_insilico: &[Spectrum], pos_gc_df_insilico: &[Spectrum], neg_lc_df_insilico: &[Spectrum], neg_gc_df_insilico: &[Spectrum], ordered_columns: Vec<String>, output_directory: &str, update: bool, progress_callback: Option<PyObject>, total_items_callback: Option<PyObject>, prefix_callback: Option<PyObject>, item_type_callback: Option<PyObject>,
) -> PyResult<()> {

    let tasks = vec![
        (pos_lc_df, "POS_LC.csv", "POS"),
        (pos_gc_df, "POS_GC.csv", "POS"),
        (neg_lc_df, "NEG_LC.csv", "NEG"),
        (neg_gc_df, "NEG_GC.csv", "NEG"),
        (pos_lc_df_insilico, "POS_LC_In_Silico.csv", "POS"),
        (pos_gc_df_insilico, "POS_GC_In_Silico.csv", "POS"),
        (neg_lc_df_insilico, "NEG_LC_In_Silico.csv", "NEG"),
        (neg_gc_df_insilico, "NEG_GC_In_Silico.csv", "NEG"),
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

            write_csv_file(data_list, &ordered_columns, &file_path, is_append, &|n| {
                let current = global_processed.fetch_add(n, Ordering::Relaxed) + n;
                if let Some(ref cb) = progress_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                }
            })?;

            if let Some(ref cb) = prefix_cb {
                Python::with_gil(|py| { let _ = cb.call1(py, (format!("Writing {} to disk...", filename),)); });
            }
            Ok(())
        }).collect()
    });

    for res in results {
        if let Err(e) = res { return Err(pyo3::exceptions::PyIOError::new_err(e)); }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    Ok(())
}

/// Cellules d'une ligne CSV pour un spectre (code inchangé, sorti de writting_csv_processing).
fn csv_row(spec: &Spectrum, ordered_columns: &[String]) -> Vec<String> {
    let mut row = Vec::with_capacity(ordered_columns.len());
    for col in ordered_columns {
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
    row
}

/// Même réglage du lecteur CSV que l'ancien code (séparateur « ; », guillemets « " »).
fn csv_builder() -> WriterBuilder {
    let mut b = WriterBuilder::new();
    b.delimiter(b';').quote(b'"');
    b
}

/// Écrit un fichier CSV (en-tête sauf en mode ajout, puis une ligne par spectre).
/// `progress(n)` est appelé après l'écriture de n spectres supplémentaires.
pub(crate) fn write_csv_file(
    data_list: &[Spectrum],
    ordered_columns: &[String],
    file_path: &Path,
    is_append: bool,
    progress: &(dyn Fn(usize) + Sync),
) -> Result<(), String> {
    let file = OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(file_path).map_err(|e| e.to_string())?;
    let mut file = std::io::BufWriter::with_capacity(1 << 20, file);

    if !is_append {
        let mut w = csv_builder().has_headers(false).from_writer(Vec::new());
        w.write_record(ordered_columns).map_err(|e| e.to_string())?;
        file.write_all(&w.into_inner().map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    }

    for round in data_list.chunks(WRITE_CHUNK * WRITE_CHUNKS_PER_ROUND) {
        let buffers: Vec<Result<Vec<u8>, String>> = round.par_chunks(WRITE_CHUNK).map(|chunk| {
            let mut w = csv_builder().has_headers(false).from_writer(Vec::new());
            for spec in chunk {
                w.write_record(&csv_row(spec, ordered_columns)).map_err(|e| e.to_string())?;
            }
            w.into_inner().map_err(|e| e.to_string())
        }).collect();
        for buf in buffers {
            file.write_all(&buf?).map_err(|e| e.to_string())?;
        }
        progress(round.len());
    }
    file.flush().map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_data::{jeu, dossier};

    // Ancienne écriture CSV (copie du code d'avant, sans les callbacks), pour comparer octet par octet.
    fn ancien_csv(data_list: &[Spectrum], ordered_columns: &[String], file_path: &Path, is_append: bool) {
        let file = OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(file_path).unwrap();
        let mut wtr = WriterBuilder::new().delimiter(b';').quote(b'"').has_headers(!is_append).from_writer(file);
        if !is_append { wtr.write_record(ordered_columns).unwrap(); }
        for spec in data_list { wtr.write_record(&csv_row(spec, ordered_columns)).unwrap(); }
        wtr.flush().unwrap();
    }

    #[test]
    fn csv_identique_a_l_ancienne_ecriture() {
        let d = dossier("csv");
        // Tailles autour des paquets : vide, 1, juste sous / sur une limite de paquet et de tour
        for &n in &[0usize, 1, 1999, 2000, 2001, 2000 * 64 + 3, 300_001] {
            let (cols, specs) = jeu(n, 42 + n as u64);
            let (a, b) = (d.join(format!("a{n}.csv")), d.join(format!("b{n}.csv")));
            ancien_csv(&specs, &cols, &a, false);
            write_csv_file(&specs, &cols, &b, false, &|_| {}).unwrap();
            // mode ajout (mise à jour)
            ancien_csv(&specs[..n / 2], &cols, &a, true);
            write_csv_file(&specs[..n / 2], &cols, &b, true, &|_| {}).unwrap();
            assert!(fs::read(&a).unwrap() == fs::read(&b).unwrap(), "CSV différent pour n = {n}");
        }
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn progression_compte_tous_les_spectres() {
        let d = dossier("prog");
        let (cols, specs) = jeu(130_005, 3);
        let total = std::sync::atomic::AtomicUsize::new(0);
        write_csv_file(&specs, &cols, &d.join("p.csv"), false, &|n| { total.fetch_add(n, Ordering::Relaxed); }).unwrap();
        assert_eq!(total.load(Ordering::Relaxed), 130_005);
        let _ = fs::remove_dir_all(&d);
    }
}
