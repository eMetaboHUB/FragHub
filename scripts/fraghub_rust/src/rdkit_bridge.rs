use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::collections::HashMap;
use rayon::prelude::*;
use crate::globals_vars::{INDIGO_SMILES_CORRECTION_PATTERN, INCHIKEY_PATTERN};
use std::fs;
use std::path::Path;

/// Orchestre le recalcul des formules et des masses à l'aide de RDKit (en C++, voir src/rdkit_native).
///
/// Pour un développeur Python : C'est ici qu'on fait le "pont" (bridge) entre les spectres et
/// RDKit. L'astuce cruciale est l'utilisation d'un `cache`
/// en RAM (HashMap). Comme on croise souvent les mêmes molécules, Rust va mémoriser les résultats
/// de RDKit pour ne jamais recalculer deux fois la même chose. C'est ce qui fait passer le
/// temps de calcul de plusieurs heures à quelques minutes.
pub fn process_mols(
    py: Python,
    mut spectrum_list: Vec<Spectrum>,
    output_directory: &str,
    deletion_report: &mut crate::deletion_report::DeletionReport,
    update: bool,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {

    if let Some(cb) = &prefix_callback { cb.call1(py, ("derivation and calculation (RDKit via Rust):",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("rows",))?; }

    let total = spectrum_list.len();
    if let Some(cb) = &total_items_callback { cb.call1(py, (total, 0))?; }

    let mut valid_list = Vec::new();
    let mut deleted_count = 0;
    let mut deleted_rows = Vec::new();
    // Mêmes colonnes que les autres fichiers de DELETED_SPECTRUMS (DATABASE_NAME manquait)
    let columns = vec!["FILENAME", "DATABASE_NAME", "FILEHASH", "PREDICTED", "SPLASH", "SPECTRUMID", "RESOLUTION", "SYNON", "IONIZATION", "MSLEVEL", "FRAGMENTATIONMETHODE", "NAME", "PRECURSORMZ", "EXACTMASS", "AVERAGEMASS", "PRECURSORTYPE", "INSTRUMENTTYPE", "INSTRUMENT", "SMILES", "INCHI", "INCHIKEY", "COLLISIONENERGY", "FORMULA", "RT", "IONMODE", "COMMENT", "ENTROPY", "CLASSYFIRE_SUPERCLASS", "CLASSYFIRE_CLASS", "CLASSYFIRE_SUBCLASS", "NPCLASS_PATHWAY", "NPCLASS_SUPERCLASS", "NPCLASS_CLASS", "NUM PEAKS", "PEAKS_LIST", "DELETION_REASON"];

    // 1. Extraire les molécules uniques de manière très rapide
    // Nouveau : inclut le PRECURSORTYPE et PRECURSORMZ pour la logique de fragmentation PY1b
    // Calcul en parallèle (sans garder le verrou Python) : chaque spectre donne au plus une clé,
    // calculée exactement comme avant, et l'ensemble final est le même (un ensemble n'a pas d'ordre).
    let unique_mols: std::collections::HashSet<(String, String, String)> = py.allow_threads(|| {
        spectrum_list.par_iter().filter_map(|spec| {
            let inchi = spec.metadata.get("INCHI").map(|s| s.as_str()).unwrap_or("");
            let smiles = spec.metadata.get("SMILES").map(|s| s.as_str()).unwrap_or("");
            let target_mol = if !inchi.is_empty() && inchi != "nan" { inchi } else { smiles };
            if !target_mol.is_empty() && target_mol != "nan" {
                let mut clean_mol = target_mol.to_string();
                if !clean_mol.contains("InChI=") {
                    clean_mol = INDIGO_SMILES_CORRECTION_PATTERN.replace_all(&clean_mol, "").to_string();
                }

                let precursortype = spec.metadata.get("PRECURSORTYPE").cloned().unwrap_or_default();
                let precursormz_raw = spec.metadata.get("PRECURSORMZ").map(|s| s.as_str()).unwrap_or("");
                let precursormz_str = if let Ok(mz) = precursormz_raw.parse::<f64>() {
                    format!("{:.4}", mz)
                } else {
                    "".to_string()
                };

                Some((clean_mol, precursortype, precursormz_str))
            } else {
                None
            }
        }).collect()
    });

    if let Some(cb) = &prefix_callback { cb.call1(py, ("RDKit calculation",))?; }

    // 2. Calcul RDKit : en C++ depuis Rust (src/rdkit_native).
    let total_unique = unique_mols.len();
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_unique, 0))?; }

    let unique_mols_vec: Vec<(String, String, String)> = unique_mols.into_iter().collect();
    let mut cache = run_rdkit(py, unique_mols_vec, &progress_callback)?;

    // ==========================================
    // LOGIQUE DE RETRAITEMENT (Fallback PubChem)
    // ==========================================
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Checking failures & fallback PubChem (Pass 2)...",))?; }

    let state = crate::global_state::STATE.read().unwrap();
    let db = &state.compound_db;
    let smiles_column = db.column_id("SMILES");
    let mut unique_mols_pass2: std::collections::HashSet<(String, String, String)> = std::collections::HashSet::new();

    for spec in &mut spectrum_list {
        let inchi = spec.metadata.get("INCHI").cloned().unwrap_or_default();
        let smiles = spec.metadata.get("SMILES").cloned().unwrap_or_default();
        let mut target_mol = if !inchi.is_empty() && inchi != "nan" { inchi } else { smiles };

        if !target_mol.is_empty() && target_mol != "nan" {
            if !target_mol.contains("InChI=") {
                target_mol = INDIGO_SMILES_CORRECTION_PATTERN.replace_all(&target_mol, "").to_string();
            }
        }

        let precursortype = spec.metadata.get("PRECURSORTYPE").cloned().unwrap_or_default();
        let precursormz_raw = spec.metadata.get("PRECURSORMZ").map(|s| s.as_str()).unwrap_or("");
        let precursormz_str = if let Ok(mz) = precursormz_raw.parse::<f64>() {
            format!("{:.4}", mz)
        } else {
            "".to_string()
        };

        let cache_key = format!("{}|{}|{}", target_mol, precursortype, precursormz_str);

        let mut is_valid = false;
        if let Some(transforms) = cache.get(&cache_key) {
            if transforms.contains_key("EXACTMASS") {
                is_valid = true;
            }
        }

        if !is_valid {
            // RDKit a échoué. On tente de trouver un SMILES de secours via l'INCHIKEY.
            let inchikey = spec.metadata.get("INCHIKEY").cloned().unwrap_or_default();
            if INCHIKEY_PATTERN.is_match(&inchikey) {
                if let Some(loc) = db.locate(&inchikey) {
                    if let Some(fallback_smiles) = smiles_column.and_then(|id| db.value(loc, id)) {
                        if !fallback_smiles.trim().is_empty() && fallback_smiles.to_lowercase() != "nan" {
                            unique_mols_pass2.insert((fallback_smiles.clone(), precursortype, precursormz_str));
                            spec.metadata.insert("FALLBACK_SMILES_PY1".to_string(), fallback_smiles);
                        }
                    }
                }
            }
        }
    }

    if !unique_mols_pass2.is_empty() {
        if let Some(cb) = &prefix_callback { cb.call1(py, ("Processing molecules (Pass 2)...",))?; }
        let total_unique2 = unique_mols_pass2.len();
        if let Some(cb) = &total_items_callback { cb.call1(py, (total_unique2, 0))?; }
        
        let unique_mols_vec2: Vec<(String, String, String)> = unique_mols_pass2.into_iter().collect();
        let cache2 = run_rdkit(py, unique_mols_vec2, &progress_callback)?;
        
        // Fusion du second cache dans le premier
        cache.extend(cache2);
    }
    // ==========================================

    if let Some(cb) = &prefix_callback { cb.call1(py, ("applying RDKit results to all spectra...",))?; }
    if let Some(cb) = &total_items_callback { cb.call1(py, (total, 0))?; }

    // 3. Appliquer les résultats en Rust à vitesse maximale
    for (i, mut spec) in spectrum_list.into_iter().enumerate() {
        let inchi = spec.metadata.get("INCHI").cloned().unwrap_or_default();
        let smiles = spec.metadata.get("SMILES").cloned().unwrap_or_default();
        let target_mol = if !inchi.is_empty() && inchi != "nan" { inchi.clone() } else { smiles.clone() };

        let mut clean_mol = String::new();
        if !target_mol.is_empty() && target_mol != "nan" {
            clean_mol = target_mol.clone();
            if !clean_mol.contains("InChI=") {
                clean_mol = INDIGO_SMILES_CORRECTION_PATTERN.replace_all(&clean_mol, "").to_string();
            }
        }

        // Si on a un fallback défini, on l'utilise à la place du clean_mol original
        if let Some(fallback) = spec.metadata.get("FALLBACK_SMILES_PY1") {
            clean_mol = fallback.clone();
            // On met à jour la colonne SMILES pour que l'export contienne le bon SMILES de PubChem
            spec.metadata.insert("SMILES".to_string(), fallback.clone());
            spec.metadata.remove("FALLBACK_SMILES_PY1");
        }

        let precursortype = spec.metadata.get("PRECURSORTYPE").cloned().unwrap_or_default();
        let precursormz_raw = spec.metadata.get("PRECURSORMZ").map(|s| s.as_str()).unwrap_or("");
        let precursormz_str = if let Ok(mz) = precursormz_raw.parse::<f64>() {
            format!("{:.4}", mz)
        } else {
            "".to_string()
        };
        
        let cache_key = format!("{}|{}|{}", clean_mol, precursortype, precursormz_str);

        if !clean_mol.is_empty() {
            if let Some(transforms) = cache.get(&cache_key) {
                for (k, v) in transforms {
                    spec.metadata.insert(k.clone(), v.clone());
                }
            }
        }

        // Vérification finale
        let final_ik = spec.metadata.get("INCHIKEY").cloned().unwrap_or_default();
        let final_em = spec.metadata.get("EXACTMASS").cloned().unwrap_or_default();
        
        if INCHIKEY_PATTERN.is_match(&final_ik) && !final_em.is_empty() && final_em != "nan" {
            valid_list.push(spec);
        } else {
            deleted_count += 1;
            spec.metadata.insert("DELETION_REASON".to_string(), "spectrum deleted because it has neither inchi nor smiles nor inchikey, even after re calculation (including PubChem fallback)".to_string());
            
            let mut row_vals = Vec::new();
            for col in &columns {
                row_vals.push(spec.metadata.get(*col).cloned().unwrap_or_default());
            }
            deleted_rows.push(row_vals);
        }

        if (i + 1) % 1000 == 0 {
            if let Some(cb) = &progress_callback { cb.call1(py, (i + 1,))?; }
        }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total,))?; }

    deletion_report.no_smiles_no_inchi_no_inchikey += deleted_count;

    // Écriture des suppressions si besoin
    if !deleted_rows.is_empty() {
        let del_dir = Path::new(output_directory).join("DELETED_SPECTRUMS");
        fs::create_dir_all(&del_dir).unwrap_or_default();
        let file_path = del_dir.join("deleted_no_inchi_smiles_inchikey_after_re_calculation.csv");

        // On preserve l'historique des runs precedents (mode append, sans reset_updates) :
        // si le fichier existe deja, on ajoute les lignes a la suite au lieu de l'ecraser.
        let is_append = update && file_path.exists();
        let file = std::fs::OpenOptions::new().write(true).create(true).append(is_append).truncate(!is_append).open(&file_path).unwrap();
        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').has_headers(!is_append).from_writer(file);
        if !is_append {
            wtr.write_record(&columns).unwrap_or_default();
        }
        for row in deleted_rows {
            wtr.write_record(&row).unwrap_or_default();
        }
        wtr.flush().unwrap_or_default();
    }

    Ok(valid_list)
}


/// Calcule les molécules avec RDKit et renvoie le dictionnaire "molécule|adduit|m/z" -> résultats.
/// Calcul en C++ (src/rdkit_native), réparti sur tous les cœurs.
fn run_rdkit(
    py: Python,
    mols: Vec<(String, String, String)>,
    progress_callback: &Option<PyObject>,
) -> PyResult<HashMap<String, HashMap<String, String>>> {
    crate::rdkit_native::run_parallel(py, mols, progress_callback.clone())
}
