// src/main_orchestrator.rs
use pyo3::prelude::*;
use pyo3::exceptions::PyException;
use std::time::SystemTime;
use chrono::Local;

// On importe toutes les fonctions du module pour pouvoir les utiliser
use crate::convertors::parsing_to_dict::parsing_to_dict_processing;
use crate::splash_generator::generate_splash_processing;
use crate::duplicatas_remover::remove_duplicatas_processing;
use crate::update_checker::check_for_update_processing;
use crate::spectrum_cleaning::spectrum_cleaning_processing;
use crate::rdkit_bridge::process_mols;
use crate::complete_from_pubchem_datas::complete_from_pubchem_datas;
use crate::ontologies_completion::ontologies_completion_processing;
use crate::de_novo_calculation::de_novo_calculation_processing;
use crate::normalize_to_not_found::normalize_to_not_found_processing;

use crate::splitter::master_splitter;

use crate::spectra_to_msp::spectra_to_msp_processing;
use crate::writers::{writting_csv_processing, writting_msp_processing, writting_json_processing};
use crate::report::generate_report_processing;
use crate::set_projects::{init_project, reset_updates};
use crate::deletion_report::DeletionReport;
use crate::mz_correction::mz_correction_processing;

/// La "Tour de Contrôle" de l'application Rust.
///
/// Pour un développeur Python : C'est l'équivalent du script principal (votre ancien `main.py`).
/// C'est cette fonction précise que l'interface Vue.js (via Python) appelle lorsqu'on clique sur "Lancer".
/// Observez comment on extrait les arguments Python (`parameters_dict`) et comment on appelle
/// régulièrement `check_stop_flag()` pour permettre à l'utilisateur d'annuler le traitement 
/// depuis l'interface à n'importe quel moment, proprement.
///
/// # Arguments
/// * `py` (Python) : Token d'accès au GIL.
/// * `parameters_dict` (&PyDict) : Les paramètres choisis par l'utilisateur dans l'UI.
/// * `progress_callback` (Option<PyObject>) : Fonction de rappel pour la progression.
/// * `total_items_callback` (Option<PyObject>) : Fonction de rappel pour le nombre total d'éléments.
/// * `prefix_callback` (Option<PyObject>) : Fonction de rappel pour l'étape en cours.
/// * `item_type_callback` (Option<PyObject>) : Fonction de rappel pour l'unité (spectra, MB, etc.).
/// * `step_callback` (Option<PyObject>) : Fonction de rappel pour le log.
/// * `completion_callback` (Option<PyObject>) : Fonction de rappel pour la fin de processus.
/// * `deletion_callback` (Option<PyObject>) : Fonction de rappel pour les rapports de suppression.
/// * `stop_flag` (Option<PyObject>) : Fonction de rappel (bool) indiquant si l'utilisateur a cliqué sur "Stop".
///
/// # Returns
/// * `PyResult<i32>` : 0 si le traitement a réussi.
#[pyfunction]
#[pyo3(signature = (parameters_dict, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None, step_callback=None, completion_callback=None, deletion_callback=None, stop_flag=None))]
pub fn main_orchestrator(
    py: Python,
    parameters_dict: &pyo3::Bound<'_, pyo3::types::PyDict>,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
    step_callback: Option<PyObject>,
    completion_callback: Option<PyObject>,
    deletion_callback: Option<PyObject>,
    stop_flag: Option<PyObject>,
) -> PyResult<i32> {

    let check_stop_flag = || -> PyResult<()> {
        if let Some(cb) = &stop_flag {
            let res = cb.call0(py)?;
            if res.is_truthy(py)? {
                return Err(PyException::new_err("Process stopped by user."));
            }
        }
        Ok(())
    };

    let start_time = SystemTime::now();

    // 1. Initialisation
    let output_directory: String = parameters_dict.get_item("output_directory")?.unwrap().extract()?;
    let input_paths_any = parameters_dict.get_item("input_directory")?.unwrap();
    let input_paths: Vec<String> = input_paths_any.extract()?;

    let input_db_names_any = parameters_dict.get_item("input_db_names")?.unwrap();
    let input_db_names: std::collections::HashMap<String, String> = input_db_names_any.extract().unwrap_or_default();

    let reset_updates_val: f64 = parameters_dict.get_item("reset_updates")?.unwrap().extract()?;
    if reset_updates_val == 1.0 {
        reset_updates(output_directory.clone())?;
    }
    init_project(output_directory.clone())?;

    let mut deletion_report = DeletionReport::default();

    check_stop_flag()?;

    // STEP 1: PARSING
    let tuple_res = parsing_to_dict_processing(
        py, input_paths.clone(), input_db_names.clone(),
        progress_callback.clone(), total_items_callback.clone(),
        prefix_callback.clone(), item_type_callback.clone(), step_callback.clone()
    )?;

    check_stop_flag()?;

    let final_msp = tuple_res.0;
    let final_csv = tuple_res.1;
    let final_json = tuple_res.2;
    let final_mgf = tuple_res.3;

    if final_msp.is_empty() && final_csv.is_empty() && final_json.is_empty() && final_mgf.is_empty() {
        if let Some(cb) = &deletion_callback { cb.call1(py, ("No files to process",))?; }
        if let Some(cb) = &completion_callback { cb.call1(py, ("--- TOTAL TIME ... ---",))?; }
        return Ok(0);
    }

    // STEP 2: SPLASH KEY
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
    if let Some(cb) = &step_callback { cb.call1(py, ("Generating SPLASH IDs",))?; }
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

    let final_msp = if !final_msp.is_empty() { generate_splash_processing(py, final_msp, "MSP".to_string(), progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())? } else { final_msp };
    let final_csv = if !final_csv.is_empty() { generate_splash_processing(py, final_csv, "CSV".to_string(), progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())? } else { final_csv };
    let final_json = if !final_json.is_empty() { generate_splash_processing(py, final_json, "JSON".to_string(), progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())? } else { final_json };
    let final_mgf = if !final_mgf.is_empty() { generate_splash_processing(py, final_mgf, "MGF".to_string(), progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())? } else { final_mgf };
    check_stop_flag()?;

    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(100)); });
    if let Some(cb) = &total_items_callback { cb.call1(py, (1,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Merging files (Please wait)...",))?; }
    if let Some(cb) = &progress_callback { cb.call1(py, (0,))?; }
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(100)); });

    // Merging into one single PyList
    let mut spectrum_list = Vec::with_capacity(final_msp.len() + final_csv.len() + final_json.len() + final_mgf.len());
    spectrum_list.extend(final_msp);
    spectrum_list.extend(final_csv);
    spectrum_list.extend(final_json);
    spectrum_list.extend(final_mgf);

    if let Some(cb) = &progress_callback { cb.call1(py, (1,))?; }
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

    let ordered_columns_str = vec!["FILENAME", "DATABASE_NAME", "FILEHASH", "PREDICTED", "SPLASH", "SPECTRUMID", "RESOLUTION", "SYNON", "IONIZATION", "MSLEVEL", "FRAGMENTATIONMETHODE", "NAME", "PRECURSORMZ", "EXACTMASS", "AVERAGEMASS", "PRECURSORTYPE", "INSTRUMENTTYPE", "INSTRUMENT", "SMILES", "INCHI", "INCHIKEY", "COLLISIONENERGY", "FORMULA", "RT", "IONMODE", "COMMENT", "ENTROPY", "CLASSYFIRE_SUPERCLASS", "CLASSYFIRE_CLASS", "CLASSYFIRE_SUBCLASS", "NPCLASS_PATHWAY", "NPCLASS_SUPERCLASS", "NPCLASS_CLASS", "NUM PEAKS", "PEAKS_LIST"];
    let ordered_columns: Vec<String> = ordered_columns_str.into_iter().map(|s| s.to_string()).collect();

    // STEP 3: DUPLICATAS
    if let Some(cb) = &step_callback { cb.call1(py, ("Removing Duplicates",))?; }
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

    let tuple_dup = remove_duplicatas_processing(
        py, spectrum_list, output_directory.clone(), ordered_columns.clone(),
        progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
    )?;
    let mut spectrum_list = tuple_dup.0;
    let deleted_count_dup = tuple_dup.1;
    {
        
        deletion_report.duplicatas_removed = deleted_count_dup;
    }
    if let Some(cb) = &deletion_callback { cb.call1(py, (format!("duplicatas removed: {}", deleted_count_dup),))?; }
    check_stop_flag()?;

    // STEP 4: UPDATES
    let mut update = false;
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
    if let Some(cb) = &step_callback { cb.call1(py, ("Checking for Updates",))?; }
    py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

    let tuple_up = check_for_update_processing(
        py, spectrum_list, output_directory.clone(), ordered_columns.clone(),
        progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
    )?;
    spectrum_list = tuple_up.0;
    let update_temp: bool = tuple_up.1;
    let deleted_count_up = tuple_up.2;
    {
        
        deletion_report.previously_cleaned = deleted_count_up;
    }
    if let Some(cb) = &deletion_callback { cb.call1(py, (format!("previously cleaned: {}", deleted_count_up),))?; }
    check_stop_flag()?;

    if !spectrum_list.is_empty() {
        if update_temp { update = true; }
        if reset_updates_val == 1.0 { update = false; }

        // STEP 5: CLEANING
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("Cleaning Spectrums",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        // Convert parameters to native map
        let mut params_f64 = std::collections::HashMap::new();
        for (k, v) in parameters_dict.iter() {
            if let Ok(key_str) = k.extract::<String>() {
                if let Ok(val_float) = v.extract::<f64>() {
                    params_f64.insert(key_str, val_float);
                }
            }
        }

        spectrum_list = spectrum_cleaning_processing(
            py, spectrum_list.clone(), output_directory.clone(), ordered_columns.clone(), &mut deletion_report, &params_f64,
            progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
        )?;

        // ---> AJOUT DU BLOC DELETION CALLBACK (COMME EN PYTHON) <---
        if let Some(cb) = &deletion_callback {
            let report = &deletion_report;
            let msg = format!(
                "\n                No peaks list: {}\n                No smiles, no inchi, no inchikey: {}\n                No precursor mz: {}\n                No or bad adduct: {}\n                Low entropy score: {}\n                Minimum peaks not required: {}\n                All peaks above precursor mz: {}\n                No peaks in mz range: {}\n                Minimum high peaks not required: {}\n                ",
                report.no_peaks_list,
                report.no_smiles_no_inchi_no_inchikey,
                report.no_precursor_mz,
                report.no_or_bad_adduct, // Intégré ici
                report.low_entropy_score,
                report.minimum_peaks_not_requiered,
                report.all_peaks_above_precursor_mz,
                report.no_peaks_in_mz_range,
                report.minimum_high_peaks_not_requiered
            );
            cb.call1(py, (msg,))?;
        }

        if spectrum_list.is_empty() {
            if let Some(cb) = &deletion_callback { cb.call1(py, ("No spectrums to process after cleaning",))?; }
            return Ok(0);
        }
        check_stop_flag()?;

        // STEP 6: MOLS DERIVATIONS (RDKit in Rust)
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("RDKit Derivation and Mass Calculation",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        spectrum_list = process_mols(
            py, spectrum_list, &output_directory, &mut deletion_report,
            progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
        )?;

        // ---> AJOUT DE LA MISE À JOUR MOLS (COMME EN PYTHON) <---
        if let Some(cb) = &deletion_callback {
            let report = &deletion_report;
            cb.call1(py, (format!("No smiles, no inchi, no inchikey (updated): {}", report.no_smiles_no_inchi_no_inchikey),))?;
        }

        check_stop_flag()?;

        // STEP 6.5: MASS & ADDUCT CORRECTION
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("Mass & Adduct Correction",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        spectrum_list = mz_correction_processing(
            py, spectrum_list, &mut deletion_report,
            progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
        )?;

        if let Some(cb) = &deletion_callback {
            let report = &deletion_report;
            cb.call1(py, (format!("Deleted by mz mismatch (NO_MATCH): {}", report.no_or_bad_adduct),))?;
        }

        check_stop_flag()?;

        // STEP 7: PUBCHEM
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("Completing from PubChem Datas",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        spectrum_list = complete_from_pubchem_datas(
            py, spectrum_list, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
        )?;
        check_stop_flag()?;

        // STEP 8: ONTOLOGIES
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("Ontologies Completion",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        spectrum_list = ontologies_completion_processing(
            py, spectrum_list, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
        )?;
        check_stop_flag()?;

        // STEP 9: DE NOVO
        let calculate_de_novo: f64 = parameters_dict.get_item("calculate_de_novo")?.map(|v| v.extract().unwrap_or(0.0)).unwrap_or(0.0);
        if calculate_de_novo == 1.0 {
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
            if let Some(cb) = &step_callback { cb.call1(py, ("De Novo Calculations",))?; }
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

            spectrum_list = de_novo_calculation_processing(
                py, spectrum_list, &params_f64, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone()
            )?;

            if let Some(cb) = &deletion_callback {
                let report = &deletion_report;
                let msg = format!(
                    "Low resolution MS2: {}\n                MS2 chemical crash (De Novo): {}",
                    report.low_resolution_ms2,
                    report.ms2_chemical_crash
                );
                cb.call1(py, (msg,))?;
            }

            check_stop_flag()?;
        }

        spectrum_list = normalize_to_not_found_processing(py, spectrum_list)?;

        // STEP 10: SPLITTING
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("Splitting Spectrums",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        let tuple_split = master_splitter(
            py,
            &spectrum_list,
            progress_callback.clone(),
            total_items_callback.clone(),
            prefix_callback.clone(),
            item_type_callback.clone()
        )?;

        let pos_lc_df = tuple_split.0;
        let pos_lc_in_silico_df = tuple_split.1;
        let pos_gc_df = tuple_split.2;
        let pos_gc_in_silico_df = tuple_split.3;
        let neg_lc_df = tuple_split.4;
        let neg_lc_in_silico_df = tuple_split.5;
        let neg_gc_df = tuple_split.6;
        let neg_gc_in_silico_df = tuple_split.7;

        check_stop_flag()?;

        // STEP 11: MSP / CSV / JSON
        let mut pos_lc_msp = Vec::new();
        let mut pos_lc_insilico_msp = Vec::new();
        let mut pos_gc_msp = Vec::new();
        let mut pos_gc_insilico_msp = Vec::new();
        let mut neg_lc_msp = Vec::new();
        let mut neg_lc_insilico_msp = Vec::new();
        let mut neg_gc_msp = Vec::new();
        let mut neg_gc_insilico_msp = Vec::new();

        let msp_val: f64 = parameters_dict.get_item("msp")?.map(|v| v.extract().unwrap_or(0.0)).unwrap_or(0.0);
        if msp_val == 1.0 {
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
            if let Some(cb) = &step_callback { cb.call1(py, ("Converting CSV to MSP",))?; }
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

            let tuple_msp = spectra_to_msp_processing(py, pos_lc_df.clone(), pos_lc_in_silico_df.clone(), pos_gc_df.clone(), pos_gc_in_silico_df.clone(), neg_lc_df.clone(), neg_lc_in_silico_df.clone(), neg_gc_df.clone(), neg_gc_in_silico_df.clone(), progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;
            pos_lc_msp = tuple_msp.0;
            pos_lc_insilico_msp = tuple_msp.1;
            pos_gc_msp = tuple_msp.2;
            pos_gc_insilico_msp = tuple_msp.3;
            neg_lc_msp = tuple_msp.4;
            neg_lc_insilico_msp = tuple_msp.5;
            neg_gc_msp = tuple_msp.6;
            neg_gc_insilico_msp = tuple_msp.7;
        }
        check_stop_flag()?;

        let csv_val: f64 = parameters_dict.get_item("csv")?.map(|v| v.extract().unwrap_or(0.0)).unwrap_or(0.0);
        if csv_val == 1.0 {
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
            if let Some(cb) = &step_callback { cb.call1(py, ("Writing CSV",))?; }
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

            writting_csv_processing(py, pos_lc_df.clone(), pos_gc_df.clone(), neg_lc_df.clone(), neg_gc_df.clone(), pos_lc_in_silico_df.clone(), pos_gc_in_silico_df.clone(), neg_lc_in_silico_df.clone(), neg_gc_in_silico_df.clone(), ordered_columns.clone(), &output_directory, update, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;
        }

        if msp_val == 1.0 {
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
            if let Some(cb) = &step_callback { cb.call1(py, ("Writing MSP",))?; }
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

            writting_msp_processing(py, pos_lc_msp.clone(), pos_lc_insilico_msp.clone(), pos_gc_msp.clone(), pos_gc_insilico_msp.clone(), neg_lc_msp.clone(), neg_lc_insilico_msp.clone(), neg_gc_msp.clone(), neg_gc_insilico_msp.clone(), &output_directory, update, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;
        }

        let json_val: f64 = parameters_dict.get_item("json")?.map(|v| v.extract().unwrap_or(0.0)).unwrap_or(0.0);
        if json_val == 1.0 {
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
            if let Some(cb) = &step_callback { cb.call1(py, ("Writing JSON",))?; }
            py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

            writting_json_processing(py, update, pos_lc_df.clone(), pos_gc_df.clone(), neg_lc_df.clone(), neg_gc_df.clone(), pos_lc_in_silico_df.clone(), pos_gc_in_silico_df.clone(), neg_lc_in_silico_df.clone(), neg_gc_in_silico_df.clone(), ordered_columns.clone(), &output_directory, progress_callback.clone(), total_items_callback.clone(), prefix_callback.clone(), item_type_callback.clone())?;


        }

        if let Some(cb) = &deletion_callback {
            let report = &deletion_report;
            // ---> CORRECTION : J'ai ajouté no_or_bad_adduct dans le total des suppressions <---
            let total_deletions = report.duplicatas_removed + report.previously_cleaned + report.no_peaks_list + report.no_smiles_no_inchi_no_inchikey + report.no_precursor_mz + report.no_or_bad_adduct + report.low_entropy_score + report.minimum_peaks_not_requiered + report.all_peaks_above_precursor_mz + report.no_peaks_in_mz_range + report.minimum_high_peaks_not_requiered;
            cb.call1(py, (format!("Total deletions: {}", total_deletions),))?;
        }

        // STEP 12: REPORT
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });
        if let Some(cb) = &step_callback { cb.call1(py, ("Generating Report",))?; }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(10)); });

        let current_datetime = Local::now().format("%d_%m_%Y__%H_%M_%S").to_string();
        let report_filename = format!("report_{}.html", current_datetime);
        let report_path = std::path::Path::new(&output_directory).join(&report_filename);
        let final_report_path_str = report_path.to_string_lossy().to_string();
        
        let elapsed_for_report = start_time.elapsed().unwrap().as_secs();
        let hours_r = elapsed_for_report / 3600;
        let mins_r = (elapsed_for_report % 3600) / 60;
        let secs_r = elapsed_for_report % 60;
        let total_time_str = format!("{:02}:{:02}:{:02}", hours_r, mins_r, secs_r);

        generate_report_processing(py, output_directory, current_datetime, total_time_str, &params_f64, &input_paths, &deletion_report, &pos_lc_df, &pos_lc_in_silico_df, &pos_gc_df, &pos_gc_in_silico_df, &neg_lc_df, &neg_lc_in_silico_df, &neg_gc_df, &neg_gc_in_silico_df)?;

        if let Some(cb) = &completion_callback {
            let elapsed = start_time.elapsed().unwrap().as_secs();
            let hours = elapsed / 3600;
            let mins = (elapsed % 3600) / 60;
            let secs = elapsed % 60;
            cb.call1(py, (format!("--- TOTAL TIME: {:02}:{:02}:{:02} ---", hours, mins, secs), final_report_path_str))?;
        }

    } else {
        if let Some(cb) = &completion_callback {
            cb.call1(py, ("Process stopped by user", ""))?;
        }
    }

    Ok(0)
}