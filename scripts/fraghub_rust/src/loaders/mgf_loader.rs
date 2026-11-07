// src/loaders/mgf_loader.rs
//! Lecture des fichiers MGF.
use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::{Duration, Instant};
use regex::Regex;

use super::generate_file_hash;

// ----------------------------------------------------------------------
// MGF FORMAT
// ----------------------------------------------------------------------
/// RUST NATIVE MGF LOADER
///
/// # Arguments
/// * `py` (Python) : Token PyO3.
/// * `mgf_file_path` (&str) : Le chemin vers le fichier MGF.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression en MB.
/// * `total_items_callback` (Option<PyObject>) : Callback pour la taille totale en MB.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément ("MB").
///
/// # Returns
/// * `PyResult<Vec<String>>` : Une liste de spectres bruts sous forme de chaînes de caractères.
pub fn load_spectrum_list_from_mgf(
    py: Python,
    mgf_file_path: &str,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<String>> {
    let file_hash = generate_file_hash(mgf_file_path);
    let path = Path::new(mgf_file_path);
    let filename = path.file_name().unwrap_or_default().to_string_lossy().to_string();

    let total_bytes = std::fs::metadata(mgf_file_path).map(|m| m.len()).unwrap_or(0);
    let total_mb = std::cmp::max(1, total_bytes / (1024 * 1024));

    if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_mb, 0)); }
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, (format!("loading [{}]:", filename),)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("MB",)); }

    let file = File::open(mgf_file_path)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 4, file);
    let mut spectrum_list = Vec::new();
    let mut buffer = vec![format!("FILENAME={}", filename)];
    let mut total_consumed: u64 = 0;
    let mut last_update = Instant::now();

    let re = Regex::new(r"(?i)FILENAME=.*\n").unwrap();
    let replacement = format!("FILENAME={}\nFILEHASH={}\n", filename, file_hash);

    let mut line = String::with_capacity(1024);
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 { break; }
        total_consumed += bytes_read as u64;

        let trimmed = line.trim();
        if trimmed == "END IONS" {
            if !buffer.is_empty() {
                let spectrum = buffer.join("\n") + "\n";
                let replaced = re.replace(&spectrum, replacement.as_str()).to_string();
                spectrum_list.push(replaced.trim_end().to_string());
                buffer.clear();
                buffer.push(format!("FILENAME={}", filename));

                if last_update.elapsed() >= Duration::from_millis(50) {
                    let current_mb = total_consumed / (1024 * 1024);
                    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (current_mb,)); }
                    py.allow_threads(|| std::thread::sleep(Duration::from_millis(1)));
                    last_update = Instant::now();
                }
            }
        } else if !trimmed.is_empty() {
            // Lignes vides ignorées : certains fichiers en mettent entre les pics ; avant, la lecture des pics
            // s'arrêtait à la première ligne vide (les pics suivants étaient perdus).
            buffer.push(trimmed.to_string());
        }
    }

    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_mb,)); }

    Ok(spectrum_list)
}
