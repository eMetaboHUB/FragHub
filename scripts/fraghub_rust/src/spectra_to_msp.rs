// src/spectra_to_msp
use pyo3::prelude::*;
use rayon::prelude::*;
use crate::spectrum::Spectrum;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::fmt::Write;

// Fonction pour extraire et nettoyer une valeur depuis le dictionnaire Python
fn get_string(spec: &Spectrum, key: &str) -> String {
    if key == "PEAKS_LIST" {
        // --- NOUVEAU : On tente d'abord de récupérer la chaîne avec les formules De Novo ---
        if let Some(val) = spec.metadata.get("PEAKS_LIST") {
            if !val.trim().is_empty() && val != "NOT FOUND" {
                let sep = if val.contains(';') { ';' } else { '\n' };
                let lines_count = val.split(sep).filter(|s| !s.trim().is_empty()).count();

                // Si la chaîne contient des lettres (formules) et correspond au bon nombre de pics, on l'utilise
                if lines_count == spec.peaks.len() && val.chars().any(|c| c.is_ascii_alphabetic()) {
                    return val.replace(";", "\n"); // On s'assure d'avoir le format MSP
                }
            }
        }

        // --- COMPORTEMENT PAR DÉFAUT : Reconstruction rapide depuis les vecteurs ---
        if spec.peaks.is_empty() { return "NOT FOUND".to_string(); }
        let mut peaks_str = String::with_capacity(spec.peaks.len() * 20);
        for (i, &(mz, int)) in spec.peaks.iter().enumerate() {
            if i > 0 { peaks_str.push('\n'); }
            peaks_str.push_str(&format!("{} {}", mz, int));
        }
        return peaks_str;
    }
    if let Some(val) = spec.metadata.get(key) {
        let trimmed = val.trim();
        if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("nan") {
            return "NOT FOUND".to_string();
        }
        return val.clone();
    }
    "NOT FOUND".to_string()
}

/// Convertit tous les spectres en texte (format MSP) pour l'écriture finale, en utilisant le multithreading.
///
/// Pour un développeur Python : Construire une énorme chaîne de texte ligne par ligne
/// en Python (ex: `texte += "\n"` ou `''.join()`) est très coûteux.
/// En Rust, on utilise `String::with_capacity(2048)` pour pré-allouer exactement la RAM
/// nécessaire à un spectre, puis la macro `write!` pousse les valeurs formatées
/// directement dans cette zone mémoire, le tout en parallèle grâce à `.par_iter()`.
#[allow(clippy::too_many_arguments)]
pub fn spectra_to_msp_processing(
    py: Python,
    pos_lc_df: Vec<Spectrum>,
    pos_lc_df_insilico: Vec<Spectrum>,
    pos_gc_df: Vec<Spectrum>,
    pos_gc_df_insilico: Vec<Spectrum>,
    neg_lc_df: Vec<Spectrum>,
    neg_lc_df_insilico: Vec<Spectrum>,
    neg_gc_df: Vec<Spectrum>,
    neg_gc_df_insilico: Vec<Spectrum>,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<(Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>)> {

    let tasks = vec![
        pos_lc_df, pos_lc_df_insilico, pos_gc_df, pos_gc_df_insilico,
        neg_lc_df, neg_lc_df_insilico, neg_gc_df, neg_gc_df_insilico
    ];

    let total_items: usize = tasks.iter().map(|l| l.len()).sum();

    if let Some(cb) = &prefix_callback { cb.call1(py, ("Converting to MSP format:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("rows",))?; }
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }

    if total_items == 0 {
        return Ok((Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()));
    }

    let global_processed = Arc::new(AtomicUsize::new(0));
    let progress_cb = progress_callback.clone();

    let results: Result<Vec<Vec<String>>, PyErr> = py.allow_threads(|| {
        let res: Vec<Vec<String>> = tasks.into_par_iter().map(|data_list| {
            if data_list.is_empty() {
                return Vec::new();
            }

            let mut spectrum_list = Vec::with_capacity(data_list.len());
            let chunk_size = 2500;

            for chunk in data_list.chunks(chunk_size) {
                let chunk_results: Vec<String> = chunk.par_iter().map(|spec| {
                    let mut spectrum = String::with_capacity(2048);

                    let _ = write!(
                        &mut spectrum,
                        "NAME: {}\nPRECURSORMZ: {}\nPRECURSORTYPE: {}\nFORMULA: {}\nINCHIKEY: {}\nINCHI: {}\nSMILES: {}\nRT: {}\nIONMODE: {}\nINSTRUMENTTYPE: {}\nINSTRUMENT: {}\nCOLLISIONENERGY: {}\nEXACTMASS: {}\nIONIZATION: {}\nMSLEVEL: {}\nDISSOCIATIONMETHOD: {}\nAVERAGEMASS: {}\nENTROPY: {}\nCLASSYFIRE_SUPERCLASS: {}\nCLASSYFIRE_CLASS: {}\nCLASSYFIRE_SUBCLASS: {}\nNPCLASS_PATHWAY: {}\nNPCLASS_SUPERCLASS: {}\nNPCLASS_CLASS: {}\nCOMMENT: FILENAME={}; DATABASE_NAME={}; FILEHASH={}; PREDICTED={}; SPLASH={}; SPECTRUMID={}; RESOLUTION={}; SYNON={}\nNUM PEAKS: {}\n",
                        get_string(spec, "NAME"),
                        get_string(spec, "PRECURSORMZ"),
                        get_string(spec, "PRECURSORTYPE"),
                        get_string(spec, "FORMULA"),
                        get_string(spec, "INCHIKEY"),
                        get_string(spec, "INCHI"),
                        get_string(spec, "SMILES"),
                        get_string(spec, "RT"),
                        get_string(spec, "IONMODE"),
                        get_string(spec, "INSTRUMENTTYPE"),
                        get_string(spec, "INSTRUMENT"),
                        get_string(spec, "COLLISIONENERGY"),
                        get_string(spec, "EXACTMASS"),
                        get_string(spec, "IONIZATION"),
                        get_string(spec, "MSLEVEL"),
                        get_string(spec, "FRAGMENTATIONMETHODE"), // Note: internal col is still FRAGMENTATIONMETHODE
                        get_string(spec, "AVERAGEMASS"),
                        get_string(spec, "ENTROPY"),
                        get_string(spec, "CLASSYFIRE_SUPERCLASS"),
                        get_string(spec, "CLASSYFIRE_CLASS"),
                        get_string(spec, "CLASSYFIRE_SUBCLASS"),
                        get_string(spec, "NPCLASS_PATHWAY"),
                        get_string(spec, "NPCLASS_SUPERCLASS"),
                        get_string(spec, "NPCLASS_CLASS"),
                        get_string(spec, "FILENAME"), get_string(spec, "DATABASE_NAME"), get_string(spec, "FILEHASH"), get_string(spec, "PREDICTED"),
                        get_string(spec, "SPLASH"), get_string(spec, "SPECTRUMID"), get_string(spec, "RESOLUTION"), get_string(spec, "SYNON"),
                        get_string(spec, "NUM PEAKS")
                    );

                    let peaks = get_string(spec, "PEAKS_LIST");
                    if peaks != "NOT FOUND" { spectrum.push_str(&peaks); }
                    spectrum.push('\n');

                    spectrum
                }).collect();

                spectrum_list.extend(chunk_results);

                let current = global_processed.fetch_add(chunk.len(), Ordering::Relaxed) + chunk.len();
                if let Some(ref cb) = progress_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                }
            }
            spectrum_list
        }).collect();

        Ok(res)
    });

    let mut final_lists = results?;

    let neg_gc_in = final_lists.pop().unwrap();
    let neg_gc = final_lists.pop().unwrap();
    let neg_lc_in = final_lists.pop().unwrap();
    let neg_lc = final_lists.pop().unwrap();
    let pos_gc_in = final_lists.pop().unwrap();
    let pos_gc = final_lists.pop().unwrap();
    let pos_lc_in = final_lists.pop().unwrap();
    let pos_lc = final_lists.pop().unwrap();

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }

    Ok((pos_lc, pos_lc_in, pos_gc, pos_gc_in, neg_lc, neg_lc_in, neg_gc, neg_gc_in))
}