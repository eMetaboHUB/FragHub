// src/writers/msp_writer.rs
//! Écriture des fichiers MSP (textes déjà mis en forme par spectra_to_msp.rs).
use pyo3::prelude::*;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

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
    py: Python, pos_lc: &[String], pos_lc_insilico: &[String], pos_gc: &[String], pos_gc_insilico: &[String], neg_lc: &[String], neg_lc_insilico: &[String], neg_gc: &[String], neg_gc_insilico: &[String], output_directory: &str, update: bool, progress_callback: Option<PyObject>, total_items_callback: Option<PyObject>, prefix_callback: Option<PyObject>, item_type_callback: Option<PyObject>,
) -> PyResult<()> {

    let tasks = vec![
        (pos_lc, "POS_LC.msp", "POS"),
        (pos_lc_insilico, "POS_LC_insilico.msp", "POS"),
        (pos_gc, "POS_GC.msp", "POS"),
        (pos_gc_insilico, "POS_GC_insilico.msp", "POS"),
        (neg_lc, "NEG_LC.msp", "NEG"),
        (neg_lc_insilico, "NEG_LC_insilico.msp", "NEG"),
        (neg_gc, "NEG_GC.msp", "NEG"),
        (neg_gc_insilico, "NEG_GC_insilico.msp", "NEG"),
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
