// src/spectrum_cleaning.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use csv::WriterBuilder;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

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

    // Les spectres sont DÉPLACÉS paquet par paquet (et non copiés) : chaque paquet est libéré
    // dès qu'il est traité, en parallèle, au lieu de libérer toute la liste d'un coup à la fin.
    let mut remaining = spectrum_list.into_iter();
    loop {
        let chunk: Vec<Spectrum> = remaining.by_ref().take(chunk_size).collect();
        if chunk.is_empty() {
            break;
        }
        let chunk_len = chunk.len();
        let results: Vec<_> = py.allow_threads(|| {
            chunk.into_par_iter().map(|spec| {
                let Spectrum { metadata: meta, peaks: spec_peaks } = spec;

                if spec_peaks.is_empty() {
                    return Err((meta, "spectrum deleted because peaks list is empty".to_string()));
                }

                let mut deletion_reason = None;
                // Appel du Normalizer AVEC le contexte
                let metadata_opt = crate::normalizer::values_normalizer::normalize_values(meta.clone(), &mut deletion_reason, &context);

                if let Some(mut valid_meta) = metadata_opt {
                    let filename = crate::globals_vars::input_filename(&valid_meta).to_string();
                    let instrument = valid_meta.get("INSTRUMENTTYPE").cloned().unwrap_or_default();
                    let is_gc = crate::globals_vars::is_gc_spectrum(&filename, &instrument);

                    let mut float_pmz: Option<f64> = None;

                    if !is_gc {
                        if let Some(pmz) = valid_meta.get("PRECURSORMZ") {
                            if let Some(caps) = crate::globals_vars::FLOAT_CHECK_PATTERN.captures(pmz) {
                                let matched_str = caps.get(1).unwrap().as_str();
                                let replaced = matched_str.replace(',', ".");
                                if let Ok(val) = replaced.parse::<f64>() {
                                    if val <= 0.0 { return Err((valid_meta, "spectrum deleted because precursor mz is less than or equal to zero.".to_string())); }
                                    // On écrit la valeur corrigée (virgule -> point), pas le texte d'origine
                                    valid_meta.insert("PRECURSORMZ".to_string(), replaced);
                                    float_pmz = Some(val);
                                } else { return Err((valid_meta, "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).".to_string())); }
                            } else { return Err((valid_meta, "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).".to_string())); }
                        } else { return Err((valid_meta, "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).".to_string())); }
                    }

                    // Filtrage des pics (les pics du spectre sont repris tels quels, sans copie)
                    // Les pics inutilisables (m/z ou intensité NaN, infinie, <= 0) sont retirés AVANT le tri :
                    // un NaN rend la comparaison incohérente et peut faire planter le tri.
                    let mut peaks = crate::peaks_filters::filters::remove_non_positive_peaks(spec_peaks);
                    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

                    // Charge du précurseur lue dans l'adduit normalisé ("[M+2H]2+" -> 2) ; 1 si l'adduit n'est pas lisible.
                    let charge = valid_meta.get("PRECURSORTYPE")
                        .and_then(|a| crate::adduct_mass_calculator::adduct_nmol_and_charge(a))
                        .map_or(1, |(_, z)| z);

                    let mut peak_del_reason = None;
                    peaks = crate::peaks_filters::filters::apply_filters(peaks, float_pmz, charge, &parameters_dict, &mut peak_del_reason);

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

        processed += chunk_len;
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
            // Mode d'ionisation inconnu : compté avec les problèmes d'adduit (même rubrique du rapport)
            crate::normalizer::check_for_bad_adduct::IONMODE_MISSING_REASON => deletion_report.no_or_bad_adduct += count,
            "spectrum deleted because its number of peaks is below the threshold chosen by the user" => deletion_report.minimum_peaks_not_requiered += count,
            "spectrum deleted because peaks list is empty after removing peaks above precursor m/z" => deletion_report.all_peaks_above_precursor_mz += count,
            "spectrum deleted because peaks list is empty after removing peaks out of mz range choiced by the user" => deletion_report.no_peaks_in_mz_range += count,
            "spectrum deleted because peaks list does not contain minimum number of high peaks required according to the value choiced by the user" => deletion_report.minimum_high_peaks_not_requiered += count,
            _ => continue,
        }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }

    // 4. Écriture des CSV de suppression.
    // Nouvelle barre de progression (en lignes écrites) : avant, la barre restait figée pendant l'écriture.
    let total_deleted: usize = deleted_spectra.values().map(|g| g.len()).sum();
    if total_deleted > 0 {
        if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_deleted,)); }
    }
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Writing deletion logs...",))?; }
    if total_deleted > 0 {
        if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("rows",)); }
    }

    if !deleted_spectra.is_empty() {
        let del_dir = Path::new(&output_directory).join("DELETED_SPECTRUMS");
        fs::create_dir_all(&del_dir)?;

        let deleted_spectra_vec: Vec<_> = deleted_spectra.into_iter().collect();
        let written = AtomicUsize::new(0);
        // L'écriture tourne dans un thread à part ; pendant ce temps, ce thread-ci (qui a le droit
        // d'appeler Python) met à jour la barre toutes les 100 ms.
        std::thread::scope(|scope| {
            let worker = scope.spawn(|| write_deletion_files(deleted_spectra_vec, &del_dir, &ordered_columns, update, &written));
            while !worker.is_finished() {
                py.allow_threads(|| std::thread::sleep(Duration::from_millis(100)));
                // let _ : comme avant, l'écriture n'est pas interrompue par le bouton stop
                if let Some(cb) = &progress_callback { let _ = cb.call1(py, (written.load(Ordering::Relaxed),)); }
            }
            if let Err(panic) = worker.join() {
                std::panic::resume_unwind(panic);
            }
        });
        if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_deleted,)); }
    }

    let final_list = py.allow_threads(|| {
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

/// Nom du fichier CSV de suppression pour un motif donné (liste inchangée).
pub(crate) fn deletion_file_name(reason: &str) -> &'static str {
    match reason {
                    "spectrum deleted because peaks list is empty" => "peaks_list_is_empty.csv",
                    "spectrum deleted because precursor mz is less than or equal to zero." => "precursor_mz_less_than_or_equal_zero.csv",
                    "spectrum deleted because it's entropy score is lower than the threshold selected by the user." => "entropy_score_lower_than_threshold.csv",
                    "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number)." => "precursor_mz_invalid_or_empty.csv",
                    "spectrum deleted because it has neither inchi nor smiles nor inchikey" => "no_inchi_smiles_or_inchikey.csv",
                    "spectrum deleted because its adduct field is empty or the value entered is not an adduct" => "adduct_empty_or_invalid.csv",
                    "spectrum deleted because the adduct corresponds to the wrong ionization mode (neg adduct in pos ionmode)." => "wrong_adduct_neg_in_pos.csv",
                    "spectrum deleted because the adduct corresponds to the wrong ionization mode (pos adduct in neg ionmode)." => "wrong_adduct_pos_in_neg.csv",
                    crate::normalizer::check_for_bad_adduct::IONMODE_MISSING_REASON => "ionmode_missing_or_invalid.csv",
                    // Étape « Mass & Adduct Correction » (mz_correction.rs)
                    crate::mz_correction::MZ_NO_MATCHING_ADDUCT_REASON => "precursor_mz_matches_no_adduct.csv",
                    crate::mz_correction::LOW_RESOLUTION_MS2_REASON => "low_resolution_ms2.csv",
                    crate::mz_correction::MS2_CHEMICAL_CRASH_REASON => "ms2_chemical_crash_test.csv",
                    "spectrum deleted because its number of peaks is below the threshold chosen by the user" => "number_of_peaks_below_threshold.csv",
                    "spectrum deleted because peaks list is empty after removing peaks above precursor m/z" => "peaks_empty_after_above_precursor_mz_removal.csv",
                    "spectrum deleted because peaks list is empty after removing peaks out of mz range choiced by the user" => "peaks_empty_after_mz_range_removal.csv",
                    "spectrum deleted because peaks list does not contain minimum number of high peaks required according to the value choiced by the user" => "insufficient_high_peaks.csv",
                    _ => "other_deletions.csv"
    }
}

/// Nombre de lignes préparées par paquet, et nombre de paquets préparés en parallèle avant
/// d'être écrits (limite la mémoire utilisée par le texte en attente d'écriture).
const ROWS_PER_CHUNK: usize = 2000;
const CHUNKS_PER_BATCH: usize = 64;

/// Écrit les CSV de DELETED_SPECTRUMS : un fichier par motif, comme avant.
///
/// Le texte des lignes est préparé en parallèle par paquets (chaque paquet avec un writer CSV
/// configuré à l'identique), puis les paquets sont écrits dans l'ordre. Un writer CSV produit
/// exactement le même texte pour une ligne quel que soit ce qu'il a écrit avant : le fichier
/// obtenu est donc identique, octet pour octet, à celui de l'ancienne écriture ligne par ligne.
pub(crate) fn write_deletion_files(
    groups: Vec<(String, Vec<HashMap<String, String>>)>,
    del_dir: &Path,
    ordered_columns: &[String],
    update: bool,
    written: &AtomicUsize,
) {
    let new_writer = |is_append: bool| {
        WriterBuilder::new().delimiter(b'\t').quote(b'"').has_headers(!is_append).from_writer(Vec::new())
    };
    groups.into_par_iter().for_each(|(reason, group)| {
        let file_path = del_dir.join(deletion_file_name(&reason));
        // On preserve l'historique des runs precedents (mode append, sans reset_updates) :
        // si le fichier existe deja, on ajoute les lignes a la suite au lieu de l'ecraser.
        let is_append = update && file_path.exists();
        let file = std::fs::OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(&file_path).unwrap();
        let mut out = BufWriter::with_capacity(4 * 1024 * 1024, file);
        if !is_append {
            let mut header = ordered_columns.to_vec();
            header.push("DELETION_REASON".to_string());
            let mut wtr = new_writer(is_append);
            wtr.write_record(&header).unwrap();
            out.write_all(&wtr.into_inner().unwrap()).unwrap();
        }

        for batch in group.chunks(ROWS_PER_CHUNK * CHUNKS_PER_BATCH) {
            let texts: Vec<Vec<u8>> = batch.par_chunks(ROWS_PER_CHUNK).map(|rows| {
                let mut wtr = new_writer(is_append);
                for meta in rows {
                    let mut record = Vec::with_capacity(ordered_columns.len() + 1);
                    for col in ordered_columns {
                        record.push(meta.get(col).map(|s| s.as_str()).unwrap_or(""));
                    }
                    record.push(meta.get("DELETION_REASON").map(|s| s.as_str()).unwrap_or(""));
                    wtr.write_record(&record).unwrap();
                }
                wtr.into_inner().unwrap()
            }).collect();
            for text in &texts {
                out.write_all(text).unwrap();
            }
            written.fetch_add(batch.len(), Ordering::Relaxed);
        }
        out.flush().unwrap();
    });
}


#[cfg(test)]
mod tests {
    use super::*;

    /// Copie exacte de l'ancienne écriture (un writer CSV directement sur le fichier, ligne par ligne).
    fn ancienne_ecriture(groups: &[(String, Vec<HashMap<String, String>>)], del_dir: &Path, ordered_columns: &[String], update: bool) {
        for (reason, group) in groups {
            let file_path = del_dir.join(deletion_file_name(reason));
            let is_append = update && file_path.exists();
            let file = std::fs::OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(&file_path).unwrap();
            let mut wtr = WriterBuilder::new().delimiter(b'\t').quote(b'"').has_headers(!is_append).from_writer(file);
            if !is_append {
                let mut header = ordered_columns.to_vec();
                header.push("DELETION_REASON".to_string());
                wtr.write_record(&header).unwrap();
            }
            for meta in group {
                let mut record = Vec::with_capacity(ordered_columns.len() + 1);
                for col in ordered_columns {
                    record.push(meta.get(col).cloned().unwrap_or_default());
                }
                record.push(meta.get("DELETION_REASON").cloned().unwrap_or_default());
                wtr.write_record(&record).unwrap();
            }
            wtr.flush().unwrap();
        }
    }

    fn jeu_de_test() -> (Vec<String>, Vec<(String, Vec<HashMap<String, String>>)>) {
        let cols: Vec<String> = ["NAME", "SMILES", "COMMENT", "PEAKS_LIST"].iter().map(|s| s.to_string()).collect();
        let pieges = ["simple", "avec\ttabulation", "avec \"guillemets\"", "multi\nligne", "", " espaces ", "é ü 中", "\"", "a;b,c"];
        let mut groups = Vec::new();
        for (g, reason) in ["spectrum deleted because peaks list is empty", "motif inconnu"].iter().enumerate() {
            let mut rows = Vec::new();
            // Plus de lignes qu'un lot complet (ROWS_PER_CHUNK * CHUNKS_PER_BATCH) pour tester les lots
            for i in 0..(ROWS_PER_CHUNK * CHUNKS_PER_BATCH + 4321 + g) {
                let mut m = HashMap::new();
                m.insert("NAME".to_string(), format!("cpd {i} {}", pieges[i % pieges.len()]));
                if i % 3 != 0 { m.insert("SMILES".to_string(), pieges[(i / 3) % pieges.len()].to_string()); }
                m.insert("PEAKS_LIST".to_string(), "10.00000000 20.00000000\n30.00000000 40.00000000".to_string());
                m.insert("DELETION_REASON".to_string(), reason.to_string());
                rows.push(m);
            }
            groups.push((reason.to_string(), rows));
        }
        (cols, groups)
    }

    #[test]
    fn csv_de_suppression_identiques_a_l_ancienne_ecriture() {
        let base = std::env::temp_dir().join(format!("fraghub_deletion_test_{}", std::process::id()));
        let (old_dir, new_dir) = (base.join("old"), base.join("new"));
        let _ = fs::remove_dir_all(&base);
        fs::create_dir_all(&old_dir).unwrap();
        fs::create_dir_all(&new_dir).unwrap();
        let (cols, groups) = jeu_de_test();

        // Premier run, puis un second en mode ajout (update = true, fichiers déjà présents)
        for update in [false, true] {
            ancienne_ecriture(&groups, &old_dir, &cols, update);
            let written = AtomicUsize::new(0);
            write_deletion_files(groups.clone(), &new_dir, &cols, update, &written);
            assert_eq!(written.load(Ordering::Relaxed), groups.iter().map(|(_, g)| g.len()).sum::<usize>());
            for (reason, _) in &groups {
                let name = deletion_file_name(reason);
                let a = fs::read(old_dir.join(name)).unwrap();
                let b = fs::read(new_dir.join(name)).unwrap();
                assert!(a == b, "fichier {name} différent (update = {update})");
            }
        }
        let _ = fs::remove_dir_all(&base);
    }
}
