// src/convertors/loaders.rs

use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use sha2::{Sha256, Digest};
use regex::Regex;
use std::time::{Duration, Instant};


/// Gets the size of a file in bytes, converts it to a string, and returns the SHA-256 hash.
///
/// Pour un développeur Python : La macro `#[pyfunction]` indique à PyO3 que cette fonction
/// Rust doit être compilée et exposée comme une fonction appelable directement depuis Python.
///
/// # Arguments
/// * `file_path` (&str) : Le chemin du fichier.
///
/// # Returns
/// * `String` : Le hash SHA-256 de la taille du fichier.
#[pyfunction]
pub fn generate_file_hash(file_path: &str) -> String {
    if let Ok(metadata) = std::fs::metadata(file_path) {
        let size = metadata.len().to_string();
        let mut hasher = Sha256::new();
        hasher.update(size.as_bytes());
        format!("{:x}", hasher.finalize())
    } else {
        format!("Error: File not found at {}", file_path)
    }
}

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

    let injection_prefix = format!(r#"{{"filename":"{}","filehash":"{}","#, filename, file_hash);

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

// ----------------------------------------------------------------------
// RUST NATIVE JSONL LOADER (Ligne par Ligne)
// ----------------------------------------------------------------------
/// RUST NATIVE JSONL LOADER (Ligne par Ligne)
///
/// # Arguments
/// * `py` (Python) : Token PyO3.
/// * `json_file_path` (&str) : Le chemin vers le fichier JSONL.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression en MB.
/// * `total_items_callback` (Option<PyObject>) : Callback pour la taille totale en MB.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément ("MB").
///
/// # Returns
/// * `PyResult<Vec<String>>` : Une liste de spectres bruts sous forme de chaînes de caractères.
pub fn load_spectrum_list_json_2(
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
    let mut line = String::with_capacity(65536);
    let mut total_consumed: u64 = 0;
    let mut last_update = Instant::now();

    let injection_prefix = format!(r#"{{"filename":"{}","filehash":"{}","#, filename, file_hash);

    loop {
        line.clear();
        match reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(bytes_read) => {
                total_consumed += bytes_read as u64;
                let trimmed = line.trim();
                if !trimmed.is_empty() && trimmed.starts_with('{') {
                    let mut injected = String::with_capacity(injection_prefix.len() + trimmed.len());
                    injected.push_str(&injection_prefix);
                    injected.push_str(&trimmed[1..]);
                    spectrum_list.push(injected);

                    if last_update.elapsed() >= Duration::from_millis(50) {
                        let current_mb = total_consumed / (1024 * 1024);
                        if let Some(cb) = &progress_callback { let _ = cb.call1(py, (current_mb,)); }
                        py.allow_threads(|| std::thread::sleep(Duration::from_millis(1)));
                        last_update = Instant::now();
                    }
                }
            }
            Err(e) => return Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
        }
    }
    
    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_mb,)); }
    
    Ok(spectrum_list)
}

// ----------------------------------------------------------------------
// MSP FORMAT
// ----------------------------------------------------------------------
/// Assemble le texte final d'un spectre MSP : en-tête FILENAME/FILEHASH + lignes du spectre.
///
/// Ancienne méthode : on écrivait "FILENAME: <nom>\n<lignes>\n", puis une regex remplaçait cette
/// première ligne par `replacement`, puis `trim_end()`. Comme la première ligne est toujours celle
/// qu'on vient d'écrire, on peut écrire `replacement` directement : même texte, sans parcourir
/// tout le spectre ni le recopier trois fois.
/// Exception : si le nom contient '$' (interprété par Regex::replace comme un groupe de capture)
/// ou un retour à la ligne, on garde l'ancienne méthode pour obtenir exactement le même texte.
fn assemble_msp_spectrum(filename: &str, replacement: &str, re: &Regex, body: &str) -> String {
    if filename.contains('$') || filename.contains('\n') {
        let spectrum = format!("FILENAME: {}\n{}\n", filename, body);
        return re.replace(&spectrum, replacement).trim_end().to_string();
    }
    let mut s = String::with_capacity(replacement.len() + body.len() + 1);
    s.push_str(replacement);
    s.push_str(body);
    s.push('\n');
    let len = s.trim_end().len();
    s.truncate(len);
    s
}

/// RUST NATIVE MSP LOADER
///
/// # Arguments
/// * `py` (Python) : Token PyO3.
/// * `msp_file_path` (&str) : Le chemin vers le fichier MSP.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression en MB.
/// * `total_items_callback` (Option<PyObject>) : Callback pour la taille totale en MB.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément ("MB").
///
/// # Returns
/// * `PyResult<Vec<String>>` : Une liste de spectres bruts sous forme de chaînes de caractères.
pub fn load_spectrum_list_from_msp(
    py: Python,
    msp_file_path: &str,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<String>> {
    let file_hash = generate_file_hash(msp_file_path);
    let path = Path::new(msp_file_path);
    let filename = path.file_name().unwrap_or_default().to_string_lossy().to_string();

    let total_bytes = std::fs::metadata(msp_file_path).map(|m| m.len()).unwrap_or(0);
    let total_mb = std::cmp::max(1, total_bytes / (1024 * 1024));

    if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_mb, 0)); }
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, (format!("loading [{}]:", filename),)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("MB",)); }

    let file = File::open(msp_file_path)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 4, file);
    let mut spectrum_list = Vec::new();
    // Lignes du spectre en cours, déjà séparées par "\n" (avant : une String par ligne puis join)
    let mut body = String::with_capacity(16 * 1024);
    let mut total_consumed: u64 = 0;
    let mut last_update = Instant::now();

    let re = Regex::new(r"(?i)FILENAME: .*\n").unwrap();
    let replacement = format!("FILENAME: {}\nFILEHASH: {}\n", filename, file_hash);

    let mut line = String::with_capacity(1024);
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 { break; }
        total_consumed += bytes_read as u64;

        let trimmed = line.trim();
        if trimmed.is_empty() {
            if !body.is_empty() {
                spectrum_list.push(assemble_msp_spectrum(&filename, &replacement, &re, &body));
                body.clear();

                if last_update.elapsed() >= Duration::from_millis(50) {
                    let current_mb = total_consumed / (1024 * 1024);
                    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (current_mb,)); }
                    py.allow_threads(|| std::thread::sleep(Duration::from_millis(1)));
                    last_update = Instant::now();
                }
            }
        } else {
            if !body.is_empty() {
                body.push('\n');
            }
            body.push_str(trimmed);
        }
    }

    if !body.is_empty() {
        spectrum_list.push(assemble_msp_spectrum(&filename, &replacement, &re, &body));
    }

    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_mb,)); }

    Ok(spectrum_list)
}

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
        } else {
            buffer.push(trimmed.to_string());
        }
    }

    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_mb,)); }

    Ok(spectrum_list)
}


#[cfg(test)]
mod tests {
    use super::assemble_msp_spectrum;
    use regex::Regex;

    /// Copie exacte de l'ancien assemblage (buffer de lignes + join + regex + trim_end).
    fn ancien(filename: &str, file_hash: &str, lines: &[&str]) -> String {
        let re = Regex::new(r"(?i)FILENAME: .*\n").unwrap();
        let replacement = format!("FILENAME: {}\nFILEHASH: {}\n", filename, file_hash);
        let mut buffer = vec![format!("FILENAME: {}", filename)];
        for l in lines {
            buffer.push(l.to_string());
        }
        let spectrum = buffer.join("\n") + "\n";
        let replaced = re.replace(&spectrum, replacement.as_str()).to_string();
        replaced.trim_end().to_string()
    }

    #[test]
    fn assemblage_identique_a_l_ancien() {
        let re = Regex::new(r"(?i)FILENAME: .*\n").unwrap();
        let noms = ["MoNA-export-LC-MS-MS_Spectra.msp", "", "nom avec espaces.msp", "accentué_é.msp",
                    "prix$1.msp", "a$$b.msp", "${0}.msp", "ligne\nnouvelle.msp", "filename: piège.msp"];
        let spectres: [&[&str]; 4] = [
            &["Name: Caffeine", "Num Peaks: 2", "138.0662 100", "195.0877 45.2"],
            &["seule ligne"],
            &["FILENAME: autre.msp", "Comments: \"filename: x\"", "1 2"],
            &["Name: $1 et \\n", "Num Peaks: 1", "10,5 20,25"],
        ];
        for nom in noms {
            let hash = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
            let replacement = format!("FILENAME: {}\nFILEHASH: {}\n", nom, hash);
            for lignes in spectres {
                let body = lignes.join("\n");
                assert_eq!(assemble_msp_spectrum(nom, &replacement, &re, &body), ancien(nom, hash, lignes),
                           "nom {nom:?}, lignes {lignes:?}");
            }
        }
    }
}
