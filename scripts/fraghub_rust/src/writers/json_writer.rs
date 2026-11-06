// src/writers/json_writer.rs
//! Écriture des fichiers JSON (un tableau de spectres par fichier).
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::fs::{self, OpenOptions};
use std::io::{Write, Seek, SeekFrom, Read};
use std::path::Path;
use once_cell::sync::Lazy;
use regex::Regex;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::{WRITE_CHUNK, WRITE_CHUNKS_PER_ROUND};

static RE_3: Lazy<Regex> = Lazy::new(|| Regex::new(r#"\[\n\s*(-?[\d\.eE\+\-]+),\n\s*(-?[\d\.eE\+\-]+),\n\s*"(.*?)"\n\s*\]"#).unwrap());
static RE_2: Lazy<Regex> = Lazy::new(|| Regex::new(r#"\[\n\s*(-?[\d\.eE\+\-]+),\n\s*(-?[\d\.eE\+\-]+)\n\s*\]"#).unwrap());

/// Regroupe les fonctions d'écriture hautement optimisées pour générer les fichiers JSON.
///
/// # Arguments
/// * `py` (Python) : Token d'accès au GIL.
/// * `update` (bool) : Indicateur pour savoir s'il faut écraser ou ajouter (append) au fichier existant.
/// * `pos_lc_df` (Vec<Spectrum>), `pos_gc_df`... : Les différentes listes de spectres formatés.
/// * `output_directory` (&str) : Dossier cible.
/// * `progress_callback`, `total_items_callback`, `prefix_callback`, `item_type_callback` : Callbacks UI.
///
/// # Returns
/// * `PyResult<()>` : Succès ou erreur IO.
#[allow(clippy::too_many_arguments)]
pub fn writting_json_processing(
    py: Python, update: bool, pos_lc_df: &[Spectrum], pos_gc_df: &[Spectrum], neg_lc_df: &[Spectrum], neg_gc_df: &[Spectrum], pos_lc_df_insilico: &[Spectrum], pos_gc_df_insilico: &[Spectrum], neg_lc_df_insilico: &[Spectrum], neg_gc_df_insilico: &[Spectrum], ordered_columns: Vec<String>, output_directory: &str, progress_callback: Option<PyObject>, total_items_callback: Option<PyObject>, prefix_callback: Option<PyObject>, item_type_callback: Option<PyObject>,
) -> PyResult<()> {

    let tasks = vec![
        (pos_lc_df, "POS_LC.json", "POS"),
        (pos_gc_df, "POS_GC.json", "POS"),
        (neg_lc_df, "NEG_LC.json", "NEG"),
        (neg_gc_df, "NEG_GC.json", "NEG"),
        (pos_lc_df_insilico, "POS_LC_In_Silico.json", "POS"),
        (pos_gc_df_insilico, "POS_GC_In_Silico.json", "POS"),
        (neg_lc_df_insilico, "NEG_LC_In_Silico.json", "NEG"),
        (neg_gc_df_insilico, "NEG_GC_In_Silico.json", "NEG"),
    ];

    let total_items: usize = tasks.iter().map(|(list, _, _)| list.len()).sum();
    if total_items == 0 { return Ok(()); }

    if let Some(cb) = &prefix_callback { cb.call1(py, ("Writing JSON files:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("rows",))?; }
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }

    let global_processed = Arc::new(AtomicUsize::new(0));
    let progress_cb = progress_callback.clone();
    let prefix_cb = prefix_callback.clone();

    let results: Vec<Result<(), String>> = py.allow_threads(|| {
        tasks.par_iter().map(|(data_list, filename, mode)| {
            if data_list.is_empty() { return Ok(()); }

            let path_dir = Path::new(output_directory).join("JSON").join(mode);
            fs::create_dir_all(&path_dir).map_err(|e| e.to_string())?;
            let file_path = path_dir.join(filename);

            write_json_file(data_list, &ordered_columns, &file_path, update, &|n| {
                let current = global_processed.fetch_add(n, Ordering::Relaxed) + n;
                if let Some(ref cb) = progress_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (current,)); });
                }
            }, &|| {
                if let Some(ref cb) = prefix_cb {
                    Python::with_gil(|py| { let _ = cb.call1(py, (format!("Writing {} to disk...", filename),)); });
                }
            })?;
            Ok(())
        }).collect()
    });

    for res in results {
        if let Err(e) = res { return Err(pyo3::exceptions::PyIOError::new_err(e)); }
    }

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }
    Ok(())
}

/// Texte JSON d'un spectre (code inchangé, sorti de writting_json_processing).
fn json_item(spec: &Spectrum, columns: &[String]) -> String {
    let mut map = serde_json::Map::new();

    for col in columns {
        if col == "PEAKS_LIST" || col == "NUM PEAKS" || col.starts_with("TREE_") { continue; }

        let mut val_str = "NaN".to_string();
        if let Some(val) = spec.metadata.get(col) {
            let clean_val = val.replace("\"", "");
            if !clean_val.is_empty() && !clean_val.eq_ignore_ascii_case("nan") { val_str = clean_val; }
        }

        if col == "MSLEVEL" { if let Ok(num) = val_str.parse::<i64>() { map.insert(col.clone(), serde_json::json!(num)); continue; } }
        if ["PRECURSORMZ", "RT", "ENTROPY"].contains(&col.as_str()) {
            if let Ok(num) = val_str.parse::<f64>() { map.insert(col.clone(), serde_json::json!(num)); continue; }
        }
        map.insert(col.clone(), serde_json::Value::String(val_str));
    }

    let num_peaks_str = spec.metadata.get("NUM PEAKS").cloned().unwrap_or_else(|| "0".to_string());
    map.insert("NUM PEAKS".to_string(), serde_json::json!(num_peaks_str.parse::<i64>().unwrap_or(0)));

    let mut peaks_array = Vec::new();

    // --- NOUVEAU : Récupération intelligente et parsing pour le format JSON ---
    if let Some(val) = spec.metadata.get("PEAKS_LIST") {
        if !val.trim().is_empty() && val != "NOT FOUND" {
            let sep = if val.contains(';') { ';' } else { '\n' };
            let lines_count = val.split(sep).filter(|s| !s.trim().is_empty()).count();

            if lines_count == spec.peaks.len() && val.chars().any(|c| c.is_ascii_alphabetic()) {
                for line in val.split(sep).filter(|s| !s.trim().is_empty()) {
                    let parts: Vec<&str> = line.trim().split_whitespace().collect();
                    if parts.len() >= 3 {
                        let mz: f64 = parts[0].parse().unwrap_or(0.0);
                        let int: f64 = parts[1].parse().unwrap_or(0.0);
                        let formula = parts[2..].join(" ");
                        peaks_array.push(serde_json::json!([mz, int, formula]));
                    } else if parts.len() == 2 {
                        let mz: f64 = parts[0].parse().unwrap_or(0.0);
                        let int: f64 = parts[1].parse().unwrap_or(0.0);
                        peaks_array.push(serde_json::json!([mz, int]));
                    }
                }
            }
        }
    }

    // Si aucune donnée de métadonnée valide, on utilise le tableau brute classique
    if peaks_array.is_empty() {
        for &(mz, intensity) in &spec.peaks {
            peaks_array.push(serde_json::json!([mz, intensity]));
        }
    }

    map.insert("PEAKS_LIST".to_string(), serde_json::Value::Array(peaks_array));

    let item_str_pretty = serde_json::to_string_pretty(&map).unwrap();
    // RE_3 compressera parfaitement tes JSON à 3 arguments !
    let compacted_1 = RE_3.replace_all(&item_str_pretty, "[$1, $2, \"$3\"]").to_string();
    let compacted_2 = RE_2.replace_all(&compacted_1, "[$1, $2]").to_string();
    format!("  {}", compacted_2.replace('\n', "\n  "))
}

/// Écrit (ou complète en mode mise à jour) un fichier JSON : tableau de spectres.
/// `progress(n)` après n spectres écrits ; `before_close()` juste avant d'écrire le « ] » final.
pub(crate) fn write_json_file(
    data_list: &[Spectrum],
    ordered_columns: &[String],
    file_path: &Path,
    update: bool,
    progress: &(dyn Fn(usize) + Sync),
    before_close: &dyn Fn(),
) -> Result<(), String> {
    let is_append_mode = update && file_path.exists() && fs::metadata(file_path).map(|m| m.len()).unwrap_or(0) > 2;

    let file_res = if is_append_mode {
        OpenOptions::new().read(true).write(true).open(file_path).and_then(|mut f| {
            let file_len = f.metadata()?.len();
            if file_len > 0 {
                let mut buf = [0u8; 1];
                for offset in 1..=std::cmp::min(10, file_len) {
                    f.seek(SeekFrom::End(-(offset as i64)))?;
                    f.read_exact(&mut buf)?;
                    if buf[0] == b']' {
                        f.set_len(file_len - offset)?;
                        f.seek(SeekFrom::End(0))?;
                        f.write_all(b",\n")?;
                        break;
                    }
                }
            }
            Ok(f)
        })
    } else {
        OpenOptions::new().write(true).create(true).truncate(true).open(file_path).and_then(|mut f| {
            f.write_all(b"[\n")?;
            Ok(f)
        })
    };

    let mut file = std::io::BufWriter::with_capacity(1 << 20, file_res.map_err(|e| e.to_string())?);
    let len = data_list.len();
    let round_size = WRITE_CHUNK * WRITE_CHUNKS_PER_ROUND;

    for (r, round) in data_list.chunks(round_size).enumerate() {
        let base = r * round_size;
        let buffers: Vec<Vec<u8>> = round.par_chunks(WRITE_CHUNK).enumerate().map(|(c, chunk)| {
            let mut out = Vec::with_capacity(chunk.len() * 2048);
            for (k, spec) in chunk.iter().enumerate() {
                let i = base + c * WRITE_CHUNK + k; // position du spectre dans le fichier
                out.extend_from_slice(json_item(spec, ordered_columns).as_bytes());
                if i < len - 1 { out.extend_from_slice(b",\n"); } else { out.extend_from_slice(b"\n"); }
            }
            out
        }).collect();
        for buf in buffers {
            file.write_all(&buf).map_err(|e| e.to_string())?;
        }
        progress(round.len());
    }

    before_close();
    file.write_all(b"]").map_err(|e| e.to_string())?;
    file.flush().map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_data::{jeu, dossier};

    // Ancienne écriture JSON (copie du code d'avant, sans les callbacks).
    fn ancien_json(data_list: &[Spectrum], ordered_columns: &[String], file_path: &Path, update: bool) {
        let is_append_mode = update && file_path.exists() && fs::metadata(file_path).map(|m| m.len()).unwrap_or(0) > 2;
        let mut file = if is_append_mode {
            let mut f = OpenOptions::new().read(true).write(true).open(file_path).unwrap();
            let file_len = f.metadata().unwrap().len();
            let mut buf = [0u8; 1];
            for offset in 1..=std::cmp::min(10, file_len) {
                f.seek(SeekFrom::End(-(offset as i64))).unwrap();
                f.read_exact(&mut buf).unwrap();
                if buf[0] == b']' { f.set_len(file_len - offset).unwrap(); f.seek(SeekFrom::End(0)).unwrap(); f.write_all(b",\n").unwrap(); break; }
            }
            f
        } else {
            let mut f = OpenOptions::new().write(true).create(true).truncate(true).open(file_path).unwrap();
            f.write_all(b"[\n").unwrap();
            f
        };
        let len = data_list.len();
        for (i, spec) in data_list.iter().enumerate() {
            file.write_all(json_item(spec, ordered_columns).as_bytes()).unwrap();
            if i < len - 1 { file.write_all(b",\n").unwrap(); } else { file.write_all(b"\n").unwrap(); }
        }
        file.write_all(b"]").unwrap();
        file.flush().unwrap();
    }

    #[test]
    fn json_identique_a_l_ancienne_ecriture() {
        let d = dossier("json");
        for &n in &[1usize, 2, 1999, 2000, 2001, 2000 * 64 + 3, 140_001] {
            let (cols, specs) = jeu(n, 7 + n as u64);
            let (a, b) = (d.join(format!("a{n}.json")), d.join(format!("b{n}.json")));
            ancien_json(&specs, &cols, &a, false);
            write_json_file(&specs, &cols, &b, false, &|_| {}, &|| {}).unwrap();
            // mode mise à jour : ajout au tableau existant
            ancien_json(&specs[..(n / 2).max(1)], &cols, &a, true);
            write_json_file(&specs[..(n / 2).max(1)], &cols, &b, true, &|_| {}, &|| {}).unwrap();
            assert!(fs::read(&a).unwrap() == fs::read(&b).unwrap(), "JSON différent pour n = {n}");
        }
        let _ = fs::remove_dir_all(&d);
    }
}
