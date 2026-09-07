// src/spectrum_cleaning.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use csv::WriterBuilder;

/// Moteur principal de nettoyage et de filtrage d'une liste de spectres.
///
/// Pour un développeur Python : Ce module coordonne toutes les étapes de nettoyage
/// (normalisation, filtres, calcul de l'entropie) pour un gros paquet (chunk) de spectres.
/// Observez l'utilisation élégante des `Result<Ok, Err>` (où `Err` capture le spectre rejeté ET son motif).
/// Les spectres rejetés sont ensuite extraits et sauvegardés directement dans des fichiers CSV "DELETED_SPECTRUMS" 
/// à la volée, sans interrompre le flux principal des spectres valides (`Ok`).
///
/// # Arguments
/// * `py` (Python) : Le token d'accès au GIL.
/// * `spectrum_list` (Vec<Spectrum>) : Liste de spectres à nettoyer.
/// * `output_directory` (String) : Répertoire de sortie (pour y sauvegarder les spectres supprimés).
/// * `ordered_columns` (Vec<String>) : Liste des colonnes pour le CSV de sortie.
/// * `deletion_report` (&mut crate::deletion_report::DeletionReport) : Structure mutable pour suivre le compteur des suppressions.
/// * `parameters_dict` (&std::collections::HashMap<String, f64>) : Paramètres utilisateurs (seuils, filtres activés...).
/// * `progress_callback` (Option<PyObject>) : Callback de progression.
/// * `total_items_callback` (Option<PyObject>) : Callback du nombre total.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément.
///
/// # Returns
/// * `PyResult<Vec<Spectrum>>` : La liste filtrée contenant uniquement les spectres valides.
pub fn spectrum_cleaning_processing(
    py: Python,
    spectrum_list: Vec<Spectrum>,
    output_directory: String,
    ordered_columns: Vec<String>,
    deletion_report: &mut crate::deletion_report::DeletionReport,
    parameters_dict: &std::collections::HashMap<String, f64>,
    update: bool,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {

    let total_items = spectrum_list.len();

    // ⚠️ ORDRE CRITIQUE POUR VUE.JS
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("cleaning spectrums:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    let context = {
        let state = crate::global_state::STATE.read().unwrap();
        crate::normalizer::NormalizerContext {
            adduct_pos: state.adduct_dict_pos.clone(),
            adduct_neg: state.adduct_dict_neg.clone(),
            adduct_massdiff_pos: state.adduct_massdiff_dict_pos.clone(),
            adduct_massdiff_neg: state.adduct_massdiff_dict_neg.clone(),
            instrument_tree: state.instrument_tree.clone(),
        }
    };

    // 2. Traitement Multithreadé (Rayon)
    let chunk_size = 2000;
    let mut processed = 0;

    // CORRECTION : On stocke maintenant un tuple (metadata, peaks) pour ne pas perdre les pics !
    let mut kept_spectra = Vec::new();
    let mut deleted_spectra: HashMap<String, Vec<HashMap<String, String>>> = HashMap::new();

    for chunk in spectrum_list.chunks(chunk_size) {
        let results: Vec<_> = py.allow_threads(|| {
            chunk.par_iter().map(|spec| {
                let meta = spec.metadata.clone();

                if spec.peaks.is_empty() {
                    return Err((meta, "spectrum deleted because peaks list is empty".to_string()));
                }

                let mut deletion_reason = None;
                // Appel du Normalizer AVEC le contexte
                let metadata_opt = crate::normalizer::values_normalizer::normalize_values(meta.clone(), &mut deletion_reason, &context);

                if let Some(mut valid_meta) = metadata_opt {
                    let filename = valid_meta.get("FILENAME").cloned().unwrap_or_default();
                    let instrument = valid_meta.get("INSTRUMENTTYPE").cloned().unwrap_or_default();
                    let is_gc = filename.contains("_GC") || crate::globals_vars::GC_PATTERN.is_match(&instrument);

                    let mut float_pmz: Option<f64> = None;

                    if !is_gc {
                        if let Some(pmz) = valid_meta.get("PRECURSORMZ") {
                            if let Some(caps) = crate::globals_vars::FLOAT_CHECK_PATTERN.captures(pmz) {
                                let matched_str = caps.get(1).unwrap().as_str();
                                let replaced = matched_str.replace(',', ".");
                                if let Ok(val) = replaced.parse::<f64>() {
                                    if val <= 0.0 { return Err((valid_meta, "spectrum deleted because precursor mz is less than or equal to zero.".to_string())); }
                                    valid_meta.insert("PRECURSORMZ".to_string(), matched_str.to_string());
                                    float_pmz = Some(val);
                                } else { return Err((valid_meta, "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).".to_string())); }
                            } else { return Err((valid_meta, "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).".to_string())); }
                        } else { return Err((valid_meta, "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).".to_string())); }
                    }

                    // Filtrage des pics
                    let mut peaks = spec.peaks.clone();
                    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

                    let mut peak_del_reason = None;
                    peaks = crate::peaks_filters::filters::apply_filters(peaks, float_pmz, &parameters_dict, &mut peak_del_reason);

                    if peaks.is_empty() {
                        let reason = peak_del_reason.unwrap_or_else(|| "spectrum deleted because peaks list is empty".to_string());
                        return Err((valid_meta, reason));
                    }

                    // Calcul de l'entropie
                    let intensities: Vec<f64> = peaks.iter().map(|p| p.1).collect();
                    let entropy = crate::peaks_filters::entropy_calculation::entropy_calculation(&intensities);
                    valid_meta.insert("ENTROPY".to_string(), format!("{:.8}", entropy));

                    if *parameters_dict.get("remove_spectrum_under_entropy_score").unwrap_or(&0.0) == 1.0 {
                        let threshold = *parameters_dict.get("remove_spectrum_under_entropy_score_value").unwrap_or(&0.0);
                        if entropy < threshold {
                            return Err((valid_meta, "spectrum deleted because it's entropy score is lower than the threshold selected by the user.".to_string()));
                        }
                    }

                    valid_meta.insert("NUM PEAKS".to_string(), peaks.len().to_string());

                    let mut formatted_peaks = String::new();
                    for (i, &(mz, int)) in peaks.iter().enumerate() {
                        if i > 0 { formatted_peaks.push('\n'); }
                        formatted_peaks.push_str(&format!("{:.8} {:.8}", mz, int));
                    }
                    valid_meta.insert("PEAKS_LIST".to_string(), formatted_peaks);

                    // CORRECTION : On retourne les métadonnées ET les pics filtrés
                    Ok((valid_meta, peaks))
                } else {
                    Err((meta, deletion_reason.unwrap_or_else(|| "Unknown deletion reason".to_string())))
                }
            }).collect()
        });

        // Répartition des succès et des erreurs
        for res in results {
            match res {
                Ok((meta, peaks)) => kept_spectra.push((meta, peaks)), // CORRECTION : on sauvegarde le tuple
                Err((mut meta, reason)) => {
                    meta.insert("DELETION_REASON".to_string(), reason.clone());
                    deleted_spectra.entry(reason).or_default().push(meta);
                }
            }
        }

        processed += chunk.len();
        if let Some(cb) = &progress_callback { cb.call1(py, (processed,))?; }
    }

    // 3. Mise à jour de L'OBJET RUST DeletionReport EN DIRECT (Nativement !)
    for (reason, group) in &deleted_spectra {
        let count = group.len();
        match reason.as_str() {
            "spectrum deleted because peaks list is empty" => deletion_report.no_peaks_list += count,
            "spectrum deleted because precursor mz is less than or equal to zero." => deletion_report.no_precursor_mz += count,
            "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number)." => deletion_report.no_precursor_mz += count,
            "spectrum deleted because it's entropy score is lower than the threshold selected by the user." => deletion_report.low_entropy_score += count,
            "spectrum deleted because it has neither inchi nor smiles nor inchikey" => deletion_report.no_smiles_no_inchi_no_inchikey += count,
            "spectrum deleted because its adduct field is empty or the value entered is not an adduct" => deletion_report.no_or_bad_adduct += count,
            "spectrum deleted because the adduct corresponds to the wrong ionization mode (neg adduct in pos ionmode)." => deletion_report.no_or_bad_adduct += count,
            "spectrum deleted because the adduct corresponds to the wrong ionization mode (pos adduct in neg ionmode)." => deletion_report.no_or_bad_adduct += count,
            "spectrum deleted because its number of peaks is below the threshold chosen by the user" => deletion_report.minimum_peaks_not_requiered += count,
            "spectrum deleted because peaks list is empty after removing peaks above precursor m/z" => deletion_report.all_peaks_above_precursor_mz += count,
            "spectrum deleted because peaks list is empty after removing peaks out of mz range choiced by the user" => deletion_report.no_peaks_in_mz_range += count,
            "spectrum deleted because peaks list does not contain minimum number of high peaks required according to the value choiced by the user" => deletion_report.minimum_high_peaks_not_requiered += count,
            _ => continue,
        }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Writing deletion logs...",))?; }

    let final_list = py.allow_threads(|| {
        // 4. Écriture des CSV de suppression
        if !deleted_spectra.is_empty() {
            let del_dir = Path::new(&output_directory).join("DELETED_SPECTRUMS");
            fs::create_dir_all(&del_dir)?;

            let deleted_spectra_vec: Vec<_> = deleted_spectra.into_iter().collect();
            deleted_spectra_vec.into_par_iter().for_each(|(reason, group)| {
                let file_name = match reason.as_str() {
                    "spectrum deleted because peaks list is empty" => "peaks_list_is_empty.csv",
                    "spectrum deleted because precursor mz is less than or equal to zero." => "precursor_mz_less_than_or_equal_zero.csv",
                    "spectrum deleted because it's entropy score is lower than the threshold selected by the user." => "entropy_score_lower_than_threshold.csv",
                    "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number)." => "precursor_mz_invalid_or_empty.csv",
                    "spectrum deleted because it has neither inchi nor smiles nor inchikey" => "no_inchi_smiles_or_inchikey.csv",
                    "spectrum deleted because its adduct field is empty or the value entered is not an adduct" => "adduct_empty_or_invalid.csv",
                    "spectrum deleted because the adduct corresponds to the wrong ionization mode (neg adduct in pos ionmode)." => "wrong_adduct_neg_in_pos.csv",
                    "spectrum deleted because the adduct corresponds to the wrong ionization mode (pos adduct in neg ionmode)." => "wrong_adduct_pos_in_neg.csv",
                    "spectrum deleted because its number of peaks is below the threshold chosen by the user" => "number_of_peaks_below_threshold.csv",
                    "spectrum deleted because peaks list is empty after removing peaks above precursor m/z" => "peaks_empty_after_above_precursor_mz_removal.csv",
                    "spectrum deleted because peaks list is empty after removing peaks out of mz range choiced by the user" => "peaks_empty_after_mz_range_removal.csv",
                    "spectrum deleted because peaks list does not contain minimum number of high peaks required according to the value choiced by the user" => "insufficient_high_peaks.csv",
                    _ => "other_deletions.csv"
                };

                let file_path = del_dir.join(file_name);
                // On preserve l'historique des runs precedents (mode append, sans reset_updates) :
                // si le fichier existe deja, on ajoute les lignes a la suite au lieu de l'ecraser.
                let is_append = update && file_path.exists();
                let file = std::fs::OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(&file_path).unwrap();
                let mut wtr = WriterBuilder::new().delimiter(b'\t').quote(b'"').has_headers(!is_append).from_writer(file);
                if !is_append {
                    let mut header = ordered_columns.clone();
                    header.push("DELETION_REASON".to_string());
                    wtr.write_record(&header).unwrap();
                }

                for meta in group {
                    let mut record = Vec::with_capacity(ordered_columns.len() + 1);
                    for col in &ordered_columns {
                        record.push(meta.get(col).cloned().unwrap_or_default());
                    }
                    record.push(meta.get("DELETION_REASON").cloned().unwrap_or_default());
                    wtr.write_record(&record).unwrap();
                }
                wtr.flush().unwrap();
            });
        }

        // 5. Reconstruction de la liste finale
        let mut final_list = Vec::with_capacity(kept_spectra.len());
        for (meta, peaks) in kept_spectra {
            let mut spec = Spectrum::default();
            spec.metadata = meta;
            spec.peaks = peaks; // CORRECTION : On réinjecte le tableau de pics dans l'objet final !
            final_list.push(spec);
        }

        Ok::<_, pyo3::PyErr>(final_list)
    })?;

    Ok(final_list)
}