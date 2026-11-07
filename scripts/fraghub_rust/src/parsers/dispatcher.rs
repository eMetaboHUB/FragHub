// src/parsers/dispatcher.rs
use pyo3::prelude::*;

use crate::loaders::{load_spectrum_list_json, load_spectrum_list_json_2, load_spectrum_list_from_msp, load_spectrum_list_from_mgf};
use crate::parsers::json_parser::json_to_dict_processing;
use crate::parsers::msp_parser::msp_to_dict_processing;
use crate::parsers::mgf_parser::mgf_to_dict_processing;
use crate::parsers::csv_parser::load_and_parse_csv;
use crate::spectrum::Spectrum;

/// Orchestrateur principal de parsing. Route les fichiers vers le bon parseur (JSON, MSP, MGF, CSV).
///
/// Pour un développeur Python : C'est la fonction "Chef d'Orchestre". Elle extrait l'état global 
/// (les clés de dictionnaire), trie les fichiers par extension, puis lance les bonnes routines
/// de parsing multi-threadées. Notez l'utilisation de `py.allow_threads(...)` qui libère le GIL
/// pour permettre à l'interface graphique (Python/Electron) de ne pas se bloquer (freezer) pendant le traitement.
///
/// # Arguments
/// * `py` (Python) : Le token d'accès au GIL Python (fourni par PyO3).
/// * `input_paths` (Vec<String>) : Liste des chemins absolus vers les fichiers à traiter.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression globale.
/// * `total_items_callback` (Option<PyObject>) : Callback pour indiquer le nombre total d'éléments.
/// * `prefix_callback` (Option<PyObject>) : Callback pour mettre à jour le message texte de l'UI.
/// * `item_type_callback` (Option<PyObject>) : Callback pour spécifier le type d'élément traité.
/// * `step_callback` (Option<PyObject>) : Callback pour indiquer l'étape en cours.
///
/// # Returns
/// * `PyResult<(Vec<Spectrum>, Vec<Spectrum>, Vec<Spectrum>, Vec<Spectrum>)>` : Un tuple de 4 listes
///   de spectres correspondant aux 4 formats de fichiers (MSP, CSV, JSON, MGF).
pub fn parsing_to_dict_processing(
    py: Python,
    input_paths: Vec<String>,
    input_db_names: std::collections::HashMap<String, String>,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
    step_callback: Option<PyObject>,
) -> PyResult<(Vec<Spectrum>, Vec<Spectrum>, Vec<Spectrum>, Vec<Spectrum>)> {

    let (keys_dict, keys_list) = {
        let state = crate::global_state::STATE.read().unwrap();
        (state.keys_dict.clone(), state.keys_list.clone())
    };

    let mut final_json = Vec::new();
    let mut final_msp = Vec::new();
    let mut final_mgf = Vec::new();

    let mut json_files = Vec::new();
    let mut msp_files = Vec::new();
    let mut mgf_files = Vec::new();
    let mut csv_files = Vec::new();

    // Tri des fichiers par extension, sans tenir compte des majuscules (".MSP", ".Mgf"... : le script
    // Python les accepte, ils doivent donc être lus ici aussi). ".jsonl" : un spectre JSON par ligne.
    for path in input_paths {
        let lower = path.to_lowercase();
        if lower.ends_with(".json") || lower.ends_with(".jsonl") { json_files.push(path); }
        else if lower.ends_with(".msp") { msp_files.push(path); }
        else if lower.ends_with(".mgf") { mgf_files.push(path); }
        else if lower.ends_with(".csv") { csv_files.push(path); }
    }

    // ======================================================
    // JSON PROCESSING
    // ======================================================
    if !json_files.is_empty() {
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { let _ = cb.call1(py, ("-- PARSING JSON TO DICT --",)); }

        for file in json_files {
            let db_name = input_db_names.get(&file).cloned().unwrap_or_else(|| "Unknown".to_string());
            // .jsonl : lu directement ligne par ligne ; .json : tableau, et ligne par ligne si ce n'est pas un tableau
            let rust_strings = if file.to_lowercase().ends_with(".jsonl") {
                load_spectrum_list_json_2(py, &file, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?
            } else {
                match load_spectrum_list_json(py, &file, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()) {
                    Ok(t) => t,
                    Err(_) => load_spectrum_list_json_2(py, &file, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?,
                }
            };

            let dict_list = json_to_dict_processing(py, rust_strings, keys_dict.clone(), keys_list.clone(), db_name, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;

            final_json.extend(dict_list);
        }
    }

    // ======================================================
    // MSP PROCESSING
    // ======================================================
    if !msp_files.is_empty() {
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { let _ = cb.call1(py, ("-- PARSING MSP TO DICT --",)); }

        for file in msp_files {
            let db_name = input_db_names.get(&file).cloned().unwrap_or_else(|| "Unknown".to_string());
            let rust_strings = load_spectrum_list_from_msp(py, &file, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;
            let dict_list = msp_to_dict_processing(py, rust_strings, keys_dict.clone(), keys_list.clone(), db_name, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;

            final_msp.extend(dict_list);
        }
    }

    // ======================================================
    // MGF PROCESSING
    // ======================================================
    if !mgf_files.is_empty() {
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { let _ = cb.call1(py, ("-- PARSING MGF TO DICT --",)); }

        for file in mgf_files {
            let db_name = input_db_names.get(&file).cloned().unwrap_or_else(|| "Unknown".to_string());
            let spectra = load_spectrum_list_from_mgf(py, &file, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;
            let dict_list = mgf_to_dict_processing(py, spectra, keys_dict.clone(), keys_list.clone(), db_name, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;
            final_mgf.extend(dict_list);
        }
    }

    // ======================================================
    // CSV PROCESSING
    // ======================================================
    let mut final_csv = Vec::new();
    if !csv_files.is_empty() {
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { let _ = cb.call1(py, ("-- PARSING CSV TO DICT --",)); }

        final_csv = load_and_parse_csv(py, csv_files, keys_dict, keys_list, input_db_names, progress_callback, total_items_callback, prefix_callback, item_type_callback)?;
    }

    Ok((final_msp, final_csv, final_json, final_mgf))
}
