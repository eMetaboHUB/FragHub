// src/loaders/json_loader.rs
//! Lecture des fichiers JSON en tableau (format GNPS).
use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::{Duration, Instant};

use super::generate_file_hash;

// ----------------------------------------------------------------------
// RUST NATIVE JSON LOADER (Array Format - GNPS)
// ----------------------------------------------------------------------
/// Lecteur de très gros fichiers JSON.
///
/// Pour un développeur Python : Au lieu de faire un `json.load()` qui chargerait
/// un fichier de 5 Go entièrement en RAM (et ferait probablement planter la machine par un OutOfMemoryError),
/// on lit le fichier "morceau par morceau" (Chunks) grâce à un `BufReader`.
/// On parcourt directement les octets (`bytes`) du disque pour extraire les spectres un à un très rapidement.
///
/// # Arguments
/// * `py` (Python) : Token PyO3.
/// * `json_file_path` (&str) : Le chemin vers le fichier JSON.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression en MB.
/// * `total_items_callback` (Option<PyObject>) : Callback pour la taille totale en MB.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément ("MB").
///
/// # Returns
/// * `PyResult<Vec<String>>` : Une liste de spectres bruts sous forme de chaînes de caractères.
pub fn load_spectrum_list_json(
    py: Python,
    json_file_path: &str,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<String>> {
    let file_hash = generate_file_hash(json_file_path);
    let path = Path::new(json_file_path);
    let filename = path.file_name().unwrap_or_default().to_string_lossy().to_string();

    let total_bytes = std::fs::metadata(json_file_path).map(|m| m.len()).unwrap_or(0);
    let total_mb = std::cmp::max(1, total_bytes / (1024 * 1024));

    if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_mb, 0)); }
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, (format!("loading [{}]:", filename),)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("MB",)); }

    let file = File::open(json_file_path)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 4, file);
    let mut spectrum_list = Vec::new();

    let mut depth = 0;
    let mut in_string = false;
    let mut escaped = false;
    let mut buffer = Vec::with_capacity(65536);
    let mut total_consumed: u64 = 0;
    let mut last_update = Instant::now();

    let injection_prefix = format!(r#"{{"{}":"{}","{}":"{}","#, super::INJECTED_FILENAME, filename, super::INJECTED_FILEHASH, file_hash);

    loop {
        let buf = match reader.fill_buf() {
            Ok(b) => b,
            Err(e) => return Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
        };
        if buf.is_empty() { break; }

        let mut consumed = 0;

        for &b in buf {
            consumed += 1;
            if depth > 0 || b == b'{' { buffer.push(b); }

            if in_string {
                if escaped { escaped = false; }
                else if b == b'\\' { escaped = true; }
                else if b == b'"' { in_string = false; }
            } else {
                match b {
                    b'"' => in_string = true,
                    b'{' => depth += 1,
                    b'}' => {
                        if depth > 0 {
                            depth -= 1;
                            if depth == 0 {
                                let s = std::str::from_utf8(&buffer).unwrap_or("");
                                if s.starts_with('{') {
                                    let mut injected = String::with_capacity(injection_prefix.len() + s.len());
                                    injected.push_str(&injection_prefix);
                                    injected.push_str(&s[1..]);
                                    spectrum_list.push(injected);
                                }
                                buffer.clear();

                                if last_update.elapsed() >= Duration::from_millis(50) {
                                    let current_mb = (total_consumed + consumed as u64) / (1024 * 1024);
                                    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (current_mb,)); }
                                    py.allow_threads(|| std::thread::sleep(Duration::from_millis(1)));
                                    last_update = Instant::now();
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        total_consumed += consumed as u64;
        reader.consume(consumed);
    }
    
    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_mb,)); }

    Ok(spectrum_list)
}
