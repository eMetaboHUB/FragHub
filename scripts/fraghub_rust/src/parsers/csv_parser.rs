// src/parsers/csv_parser.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

// Import du générateur de hash
use crate::loaders::generate_file_hash;

// 1. Détection du séparateur (comme dans votre Python)
/// Détecte si le fichier CSV utilise des tabulations (\t) ou des points-virgules (;).
///
/// Pour un développeur Python : L'utilisation de `File::open` renvoie un `Result`. 
/// `if let Ok(file)` permet d'ouvrir le fichier seulement si tout s'est bien passé,
/// sinon on passe à la suite sans faire planter le programme avec une exception.
///
/// # Arguments
/// * `file_path` (&str) : Le chemin vers le fichier CSV.
///
/// # Returns
/// * `u8` : L'octet correspondant au séparateur (`b'\t'` ou `b';'`).
fn detect_separator(file_path: &str) -> u8 {
    if let Ok(file) = File::open(file_path) {
        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        if reader.read_line(&mut first_line).is_ok() {
            let tab_count = first_line.matches('\t').count();
            let semi_count = first_line.matches(';').count();
            if tab_count > semi_count {
                return b'\t'; // `b'\t'` est un octet (byte), pas un caractère (char).
            }
        }
    }
    b';'
}

// 2. Parsing exact des pics selon VOTRE logique regex
/// Parse la colonne des pics (sous forme de chaîne JSON) en un Vecteur de tuples (m/z, intensité).
/// 
/// Pour un développeur Python : La boucle parcourt les résultats d'une regex.
/// Le `if let (Ok(mz), Ok(intensity))` vérifie simultanément que la masse et l'intensité 
/// peuvent bien être converties de `str` à `float` (`f64`). Si une valeur échoue (ex: texte malformé),
/// le point de donnée est ignoré silencieusement.
///
/// # Arguments
/// * `peak_list_string` (&str) : La chaîne de caractères brute de la liste des pics.
///
/// # Returns
/// * `Vec<(f64, f64)>` : La liste des tuples `(m/z, intensité)`.
fn parse_peak_list_native(peak_list_string: &str) -> Vec<(f64, f64)> {
    let mut peaks = Vec::new();
    for cap in crate::globals_vars::PEAK_LIST_JSON_PATTERN.captures_iter(peak_list_string) {
        let mz_str = cap[1].replace(",", ".");
        let int_str = cap[2].replace(",", ".");
        if let (Ok(mz), Ok(intensity)) = (mz_str.parse::<f64>(), int_str.parse::<f64>()) {
            peaks.push((mz, intensity)); // Utilisation de tuples pour optimiser la RAM
        }
    }
    peaks
}

// 3. La fonction principale qui remplace Pandas et csv_to_dict_processing
/// La fonction principale qui remplace l'utilisation de Pandas et `csv_to_dict_processing` de Python.
///
/// Pour un développeur Python : Cette fonction est appelée directement depuis Python (grâce à `PyResult`).
/// `py: Python` est un jeton ("token") prouvant que Rust possède le GIL Python. 
/// Il permet d'appeler de façon sécurisée les callbacks Python (pour la barre de progression).
///
/// # Arguments
/// * `py` (Python) : Le token PyO3.
/// * `csv_files` (Vec<String>) : Liste des chemins vers les fichiers CSV.
/// * `keys_dict` (HashMap<String, String>) : Dictionnaire de mapping des clés.
/// * `keys_list` (Vec<String>) : Liste des clés officielles à conserver.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression.
/// * `total_items_callback` (Option<PyObject>) : Callback pour le total.
/// * `prefix_callback` (Option<PyObject>) : Callback pour le préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback pour le type d'élément.
///
/// # Returns
/// * `PyResult<Vec<Spectrum>>` : La liste consolidée de tous les spectres parsés.
pub fn load_and_parse_csv(
    py: Python,
    csv_files: Vec<String>,
    keys_dict: HashMap<String, String>,
    keys_list: Vec<String>,
    input_db_names: std::collections::HashMap<String, String>,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {

    let total_files = csv_files.len();
    if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_files, 0)); }
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, ("Reading CSV files:",)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("csv_files",)); }

    let mut result_list = Vec::new();
    let mut processed_files = 0;

    for file_path in csv_files {
        let db_name = input_db_names.get(&file_path).cloned().unwrap_or_else(|| "Unknown".to_string());
        let file_hash = generate_file_hash(&file_path);
        let filename = std::path::Path::new(&file_path)
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();

        let separator = detect_separator(&file_path);

        // COMPTAGE RAPIDE DES LIGNES (très rapide via BufReader)
        let total_records = py.allow_threads(|| {
            let f = std::fs::File::open(&file_path).map_err(|e| e.to_string())?;
            let mut reader = std::io::BufReader::with_capacity(1024 * 1024, f);
            let mut count = 0;
            let mut buf = Vec::new();
            while let Ok(n) = std::io::BufRead::read_until(&mut reader, b'\n', &mut buf) {
                if n == 0 { break; }
                count += 1;
                buf.clear();
            }
            Ok::<usize, String>(if count > 0 { count - 1 } else { 0 })
        }).map_err(|e| pyo3::exceptions::PyIOError::new_err(e))?;

        if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_records, 0)); }
        if let Some(cb) = &prefix_callback { let _ = cb.call1(py, (format!("Parsing {}", filename),)); }
        if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("rows",)); }

        let mut local_result = py.allow_threads(|| {
            let mut inner_result_list = Vec::new();
            
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(separator)
                .has_headers(true)
                .from_path(&file_path)
                .map_err(|e| e.to_string())?;

            let headers: Vec<String> = rdr.headers()
                .map_err(|e| e.to_string())?
                .iter()
                .map(|h| h.to_lowercase())
                .collect();

            let mut current_row = 0;
            let mut last_update = std::time::Instant::now();

            for result in rdr.records() {
                current_row += 1;
                let record = result.map_err(|e| e.to_string())?;
                let mut spec = Spectrum::default();

                spec.metadata.insert("FILENAME".to_string(), filename.clone());
                // Nom du fichier d'entrée, gardé en interne (voir globals_vars::INPUT_FILENAME_KEY)
                spec.metadata.insert(crate::globals_vars::INPUT_FILENAME_KEY.to_string(), filename.clone());
                let mut original_filename: Option<String> = None;
                spec.metadata.insert("FILEHASH".to_string(), file_hash.clone());
                spec.metadata.insert("DATABASE_NAME".to_string(), db_name.clone());

                for (i, field) in record.iter().enumerate() {
                    if i >= headers.len() { continue; }
                    let header = &headers[i];

                    if header == "peaks" || header == "peaks_list" {
                        spec.peaks = parse_peak_list_native(field);
                        continue;
                    }

                    if let Some(mapped_key) = keys_dict.get(header) {
                        // FILEHASH : celui du fichier CSV. FILENAME : une colonne remplie donne le nom du fichier
                        // de la base d'origine (la première), sinon c'est le nom du fichier CSV
                        if mapped_key == "FILEHASH" { continue; }
                        if mapped_key == "FILENAME" {
                            let v = field.trim();
                            if original_filename.is_none() && !crate::globals_vars::EMPTY_PATTERN.is_match(v) {
                                original_filename = Some(v.to_string());
                            }
                            continue;
                        }
                        if keys_list.contains(mapped_key) {
                            // Unité portée par le nom de colonne (ex. RTINSECONDS -> "70.99 s"), comme pour MSP/MGF/JSON
                            spec.metadata.insert(mapped_key.clone(),
                                crate::parsers::raw_metadata::keep_unit_from_key(header, field.to_string()));
                        }
                    }
                }
                if let Some(original) = original_filename {
                    spec.metadata.insert("FILENAME".to_string(), original);
                }

                for key in &keys_list {
                    if !spec.metadata.contains_key(key) && key != "PEAKS_LIST" {
                        spec.metadata.insert(key.clone(), "".to_string());
                    }
                }

                inner_result_list.push(spec);

                // Update progress every 500ms
                if last_update.elapsed() > std::time::Duration::from_millis(500) {
                    if let Some(cb) = &progress_callback {
                        // cb is a PyObject, but we are inside allow_threads, so we need a Python block
                        // Wait, we can't acquire the GIL here easily unless we use Python::with_gil
                        Python::with_gil(|py| {
                            let _ = cb.call1(py, (current_row,));
                        });
                    }
                    last_update = std::time::Instant::now();
                }
            }
            Ok::<Vec<Spectrum>, String>(inner_result_list)
        }).map_err(|e| pyo3::exceptions::PyIOError::new_err(e))?;
        
        result_list.append(&mut local_result);

        processed_files += 1;
        if let Some(cb) = &progress_callback { let _ = cb.call1(py, (processed_files,)); }

        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(1)); });
    }

    Ok(result_list)
}