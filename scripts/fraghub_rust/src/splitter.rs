// src/splitter.rs
use pyo3::prelude::*;
use rayon::prelude::*;
use crate::spectrum::Spectrum;

/// Trie les spectres dans 8 catégories (LC/GC, POS/NEG, IN_SILICO/EXP).
///
/// Pour un développeur Python : Voici une merveille de "Programmation Fonctionnelle" en Rust.
/// On utilise `par_iter().fold(...).reduce(...)`. C'est le principe ultra-rapide de "Map-Reduce".
/// Chaque cœur CPU trie son propre petit paquet de spectres (`fold`), puis Rust fusionne (`reduce`)
/// tous ces paquets à la fin. Aucun thread ne bloque un autre, la RAM n'est pas fragmentée !
pub fn master_splitter(
    py: Python,
    spectrum_list: &Vec<Spectrum>,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<(
    Vec<Spectrum>, Vec<Spectrum>, // POS_LC_EXP, POS_LC_INSILICO
    Vec<Spectrum>, Vec<Spectrum>, // POS_GC_EXP, POS_GC_INSILICO
    Vec<Spectrum>, Vec<Spectrum>, // NEG_LC_EXP, NEG_LC_INSILICO
    Vec<Spectrum>, Vec<Spectrum>  // NEG_GC_EXP, NEG_GC_INSILICO
)> {
    if let Some(cb) = &prefix_callback { cb.call1(py, ("Splitting spectra into categories:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    let total_items = spectrum_list.len();
    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }

    let result = py.allow_threads(|| {
        spectrum_list
            .par_iter()
            .fold(
                || (
                    Vec::new(), Vec::new(), Vec::new(), Vec::new(),
                    Vec::new(), Vec::new(), Vec::new(), Vec::new()
                ),
                |mut acc, spec| {
                    let ionmode = spec.metadata.get("IONMODE").map(|s| s.as_str()).unwrap_or("").to_lowercase();

                    // On ignore si ce n'est ni positif ni négatif (comme dans le code initial)
                    if ionmode == "positive" || ionmode == "negative" {
                        let is_pos = ionmode == "positive";
                        let instr = spec.metadata.get("INSTRUMENTTYPE").map(|s| s.as_str()).unwrap_or("");
                        let filename = spec.metadata.get("FILENAME").map(|s| s.as_str()).unwrap_or("");
                        let is_gc = crate::globals_vars::is_gc_spectrum(filename, instr);
                        let pred = spec.metadata.get("PREDICTED").map(|s| s.as_str()).unwrap_or("").to_lowercase();
                        let is_insilico = pred == "true";

                        // On range directement dans la bonne combinaison
                        match (is_pos, is_gc, is_insilico) {
                            (true, false, false)  => acc.0.push(spec.clone()),
                            (true, false, true)   => acc.1.push(spec.clone()),
                            (true, true, false)   => acc.2.push(spec.clone()),
                            (true, true, true)    => acc.3.push(spec.clone()),
                            (false, false, false) => acc.4.push(spec.clone()),
                            (false, false, true)  => acc.5.push(spec.clone()),
                            (false, true, false)  => acc.6.push(spec.clone()),
                            (false, true, true)   => acc.7.push(spec.clone()),
                        }
                    }
                    acc
                }
            )
            .reduce(
                || (
                    Vec::new(), Vec::new(), Vec::new(), Vec::new(),
                    Vec::new(), Vec::new(), Vec::new(), Vec::new()
                ),
                |mut a, mut b| {
                    // On fusionne les résultats de tous les threads très rapidement
                    a.0.append(&mut b.0);
                    a.1.append(&mut b.1);
                    a.2.append(&mut b.2);
                    a.3.append(&mut b.3);
                    a.4.append(&mut b.4);
                    a.5.append(&mut b.5);
                    a.6.append(&mut b.6);
                    a.7.append(&mut b.7);
                    a
                }
            )
    });

    if let Some(cb) = &progress_callback { cb.call1(py, (total_items,))?; }

    Ok(result)
}