// src/de_novo_calculation.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use once_cell::sync::Lazy;
use regex::Regex;

static FORMULA_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"([A-Z][a-z]*)(\d*)").unwrap());

// Dictionnaire des masses exactes (Strictement identique à atoms_of_life)
fn get_atom_mass(symbol: &str) -> Option<f64> {
    match symbol {
        "H" => Some(1.0078250322), "C" => Some(12.000000), "N" => Some(14.003074004),
        "O" => Some(15.994914619), "F" => Some(18.998403162), "Na" => Some(22.98976928),
        "Mg" => Some(23.98504170), "P" => Some(30.973761998), "S" => Some(31.972071174),
        "Cl" => Some(34.9688527), "K" => Some(38.96370649), "Ca" => Some(39.9625909),
        "Mn" => Some(54.938043), "Fe" => Some(55.934936), "Co" => Some(58.933194),
        "Cu" => Some(62.929597), "Zn" => Some(63.929142), "Br" => Some(78.918338),
        "Se" => Some(79.916522), "I" => Some(126.90447),
        _ => None,
    }
}

// Récupération des valences pour le DBE (Cas Complexes)
fn get_valence(symbol: &str) -> i32 {
    match symbol {
        "C" => 4, "N" => 3, "P" => 3, "H" => 1, "F" => 1, "Cl" => 1,
        "Br" => 1, "I" => 1, "Na" => 1, "K" => 1, _ => 2,
    }
}

// Emulation STRICTE du Banker's Rounding de Python (round-half-to-even)
fn py_round(val: f64, decimals: i32) -> f64 {
    let multiplier = 10_f64.powi(decimals);
    let v = val * multiplier;
    let v_abs = v.abs();
    let floor = v_abs.floor();
    let diff = v_abs - floor;

    let rounded_abs = if diff < 0.5 {
        floor
    } else if diff > 0.5 {
        floor + 1.0
    } else {
        if floor as i64 % 2 == 0 { floor } else { floor + 1.0 }
    };

    if val < 0.0 { -rounded_abs / multiplier } else { rounded_abs / multiplier }
}

// Emulation parfaite de str(float) de Python
fn py_format_float(f: f64) -> String {
    let mut s = format!("{}", f);
    if !s.contains('.') && !s.contains('e') && !s.contains('E') {
        s.push_str(".0");
    }
    s
}

// Fonction récursive mathématiquement identique à l'algorithme Numba
fn generate_combinations(
    element_index: usize,
    current_mass: f64,
    composition: &mut [u32],
    element_masses: &[f64],
    max_counts: &[u32],
    min_mass: f64,
    max_mass: f64,
    max_remaining_masses: &[f64],
    results: &mut Vec<Vec<u32>>
) {
    if current_mass + max_remaining_masses[element_index] < min_mass { return; }

    if element_index == element_masses.len() {
        if current_mass >= min_mass {
            results.push(composition.to_vec());
        }
        return;
    }

    let element_mass = element_masses[element_index];
    for count in 0..=max_counts[element_index] {
        let new_mass = current_mass + (count as f64) * element_mass;
        if new_mass > max_mass { break; }

        composition[element_index] = count;
        generate_combinations(
            element_index + 1, new_mass, composition, element_masses,
            max_counts, min_mass, max_mass, max_remaining_masses, results
        );
    }
    composition[element_index] = 0; // Backtracking parfait
}

// Formatage Strict Hill System
fn format_hill_system(composition: &HashMap<&str, u32>) -> String {
    let mut result = String::new();
    let mut comps = composition.clone();

    if let Some(&c) = comps.get("C") {
        if c > 0 { result.push_str(&format!("C{}", if c > 1 { c.to_string() } else { "".to_string() })); }
        comps.remove("C");
    }
    if let Some(&h) = comps.get("H") {
        if h > 0 { result.push_str(&format!("H{}", if h > 1 { h.to_string() } else { "".to_string() })); }
        comps.remove("H");
    }

    let mut keys: Vec<_> = comps.keys().collect();
    keys.sort();
    for k in keys {
        let count = comps[k];
        if count > 0 {
            result.push_str(&format!("{}{}", k, if count > 1 { count.to_string() } else { "".to_string() }));
        }
    }
    result
}

// Moteur de traitement d'un spectre individuel
fn process_spectrum_peaks(formula: &str, peaks_list_str: &str, ppm_tol: f64, ion_mode: &str) -> String {
    let mut max_comp: HashMap<String, u32> = HashMap::new();
    let mut present_elements: HashSet<String> = HashSet::new();

    for caps in FORMULA_PATTERN.captures_iter(formula) {
        let el = caps.get(1).unwrap().as_str().to_string();
        let count = caps.get(2).unwrap().as_str().parse::<u32>().unwrap_or(1);

        if get_atom_mass(&el).is_some() {
            max_comp.insert(el.clone(), count);
        }
        present_elements.insert(el); // Sauvegarde pour la vérification is_simple
    }

    let mut elements_to_test: Vec<String> = max_comp.keys().cloned().collect();
    if elements_to_test.is_empty() { return peaks_list_str.to_string(); }

    elements_to_test.sort_by(|a, b| {
        let mass_a = get_atom_mass(a).unwrap();
        let mass_b = get_atom_mass(b).unwrap();
        mass_a.partial_cmp(&mass_b).unwrap().reverse()
    });

    let simple_atoms = ["C", "H", "N", "O", "P", "S"];
    let is_simple = present_elements.iter().all(|el| simple_atoms.contains(&el.as_str()));

    let c_idx = elements_to_test.iter().position(|x| x == "C");
    let h_idx = elements_to_test.iter().position(|x| x == "H");
    let n_idx = elements_to_test.iter().position(|x| x == "N");
    let o_idx = elements_to_test.iter().position(|x| x == "O");
    let p_idx = elements_to_test.iter().position(|x| x == "P");
    let s_idx = elements_to_test.iter().position(|x| x == "S");

    let element_masses: Vec<f64> = elements_to_test.iter().map(|el| get_atom_mass(el).unwrap()).collect();
    let proton_mass = get_atom_mass("H").unwrap();

    let mut new_peaks_list = String::with_capacity(peaks_list_str.len() + 200);
    let mut has_annotations = false;

    // Split strict (toujours sur saut de ligne pour simuler la fonction Python qui finit par join("\n"))
    let lines: Vec<&str> = peaks_list_str.trim().split('\n').collect();
    let mut is_first_line = true;

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 { continue; }

        let mz_str = parts[0];
        let intensity_str = parts[1];
        let mz: f64 = mz_str.parse().unwrap_or(0.0);
        let is_neg = ion_mode.to_uppercase().contains("N") || ion_mode.to_uppercase().contains("-");
        let target_mass = if is_neg { mz + proton_mass } else { mz - proton_mass };

        let min_mass = target_mass * (1.0 - ppm_tol / 1_000_000.0);
        let max_mass = target_mass * (1.0 + ppm_tol / 1_000_000.0);

        let mut best_formula = String::new();

        if target_mass > 0.0 {
            let max_counts: Vec<u32> = elements_to_test.iter().enumerate().map(|(idx, el)| {
                let max_by_mass = (max_mass / element_masses[idx]) as u32;
                std::cmp::min(*max_comp.get(el).unwrap_or(&0), max_by_mass)
            }).collect();

            let mut max_remaining_masses = vec![0.0; elements_to_test.len() + 1];
            for i in (0..elements_to_test.len()).rev() {
                max_remaining_masses[i] = max_remaining_masses[i + 1] + (max_counts[i] as f64 * element_masses[i]);
            }

            let mut results = Vec::new();
            let mut composition = vec![0; elements_to_test.len()];
            generate_combinations(0, 0.0, &mut composition, &element_masses, &max_counts, min_mass, max_mass, &max_remaining_masses, &mut results);

            for comp in results {
                let c = c_idx.map_or(0, |idx| comp[idx]);
                let h = h_idx.map_or(0, |idx| comp[idx]);
                let n = n_idx.map_or(0, |idx| comp[idx]);
                let p = p_idx.map_or(0, |idx| comp[idx]);

                let current_mass: f64 = comp.iter().zip(&element_masses).map(|(&cnt, &m)| (cnt as f64) * m).sum();

                if is_simple {
                    if c == 0 { continue; }
                    let double_dbe = 2 * (c as i32) - (h as i32) + (n as i32) + (p as i32) + 2;
                    if double_dbe < 0 || double_dbe % 2 != 0 { continue; }
                    if 5 * h < c || h > 3 * c { continue; }
                } else {
                    if c_idx.is_some() && c == 0 { continue; }

                    let mut contains_special = false;
                    for (idx, &cnt) in comp.iter().enumerate() {
                        if cnt > 0 {
                            let el = elements_to_test[idx].as_str();
                            if !["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"].contains(&el) {
                                contains_special = true;
                                break;
                            }
                        }
                    }
                    if !contains_special {
                        let rounded_mass = py_round(current_mass, 0) as i64;
                        let n_count = n as i64 + p as i64;
                        if (rounded_mass % 2).abs() != (n_count % 2) { continue; }
                    }

                    let mut double_dbe_contrib = 0;
                    for (idx, &cnt) in comp.iter().enumerate() {
                        if cnt > 0 {
                            let valence = get_valence(elements_to_test[idx].as_str());
                            double_dbe_contrib += (cnt as i32) * (valence - 2);
                        }
                    }
                    let double_dbe = 2 + double_dbe_contrib;
                    if double_dbe < 0 || double_dbe % 2 != 0 { continue; }

                    if c > 0 {
                        if 5 * h < c || h > 3 * c { continue; }
                        let o = o_idx.map_or(0, |idx| comp[idx]);
                        let s = s_idx.map_or(0, |idx| comp[idx]);
                        if 5 * o > 6 * c || 10 * n > 13 * c || 10 * p > 3 * c || 5 * s > 4 * c { continue; }
                    }
                }

                let mut ion_comp = HashMap::new();
                for (idx, &cnt) in comp.iter().enumerate() {
                    ion_comp.insert(elements_to_test[idx].as_str(), cnt);
                }
                if is_neg {
                    let entry = ion_comp.entry("H").or_insert(0);
                    if *entry > 0 { *entry -= 1; }
                } else {
                    *ion_comp.entry("H").or_insert(0) += 1;
                }

                let error_ppm = ((current_mass - target_mass) / target_mass) * 1_000_000.0;

                // Emulation totale du Python: round(error, 2) + str()
                let rounded_error = py_round(error_ppm, 2);
                let formatted_error = py_format_float(rounded_error);

                best_formula = format!("{}+/{}", format_hill_system(&ion_comp), formatted_error);
                has_annotations = true;
                break;
            }
        }

        if !is_first_line { new_peaks_list.push('\n'); }
        is_first_line = false;

        if best_formula.is_empty() {
            new_peaks_list.push_str(&format!("{} {} ", mz_str, intensity_str));
        } else {
            new_peaks_list.push_str(&format!("{} {} {}", mz_str, intensity_str, best_formula));
        }
    }

    if !has_annotations {
        return peaks_list_str.to_string();
    }

    new_peaks_list
}

/// Calcule et attribue les formules chimiques aux fragments (algorithme de-novo).
///
/// Pour un développeur Python : C'est la fonction la plus complexe mathématiquement (Backtracking/Combinaisons).
/// L'implémentation Rust émule parfaitement l'algorithme "Numba" utilisé côté Python,
/// mais sans le fameux GIL (Global Interpreter Lock). En utilisant `par_iter_mut()`,
/// la charge colossale du calcul "De Novo" est répartie sur tous les cœurs du processeur Mac.
pub fn de_novo_calculation_processing(
    py: Python,
    mut spectrum_list: Vec<Spectrum>,
    parameters_dict: &HashMap<String, f64>,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {

    if let Some(cb) = &prefix_callback { cb.call1(py, ("Calculating de novo formulas:",))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    let total_items = spectrum_list.len();

    if let Some(cb) = &total_items_callback { cb.call1(py, (total_items, 0))?; }
    
    let ppm_tol = parameters_dict.get("de_novo_ppm_tolerance").cloned().unwrap_or(5.0);

    let chunk_size = 500;
    let mut processed = 0;

    for chunk in spectrum_list.chunks_mut(chunk_size) {
        py.allow_threads(|| {
            chunk.par_iter_mut().for_each(|spec| {
                let formula = spec.metadata.get("FORMULA").cloned().unwrap_or_default();
                let ion_mode = spec.metadata.get("IONMODE").cloned().unwrap_or_default();

                // Format the native peaks back to String to reuse the de novo logic (temporarily)
                let mut peaks_list = String::with_capacity(spec.peaks.len() * 20);
                for (i, &(mz, int)) in spec.peaks.iter().enumerate() {
                    if i > 0 { peaks_list.push('\n'); }
                    peaks_list.push_str(&format!("{} {}", mz, int));
                }

                if !formula.is_empty() && !peaks_list.is_empty() && peaks_list != "nan" {
                    let updated_peaks = process_spectrum_peaks(&formula, &peaks_list, ppm_tol, &ion_mode);
                    spec.metadata.insert("PEAKS_LIST".to_string(), updated_peaks);
                }
            });
        });

        processed += chunk.len();
        if let Some(cb) = &progress_callback { cb.call1(py, (processed,))?; }
    }

    Ok(spectrum_list)
}

pub fn is_valid_peak(formula: &str, target_mass: f64, ppm_tol: f64) -> bool {
    let mut max_comp: HashMap<String, u32> = HashMap::new();
    let mut present_elements: HashSet<String> = HashSet::new();

    for caps in FORMULA_PATTERN.captures_iter(formula) {
        let el = caps.get(1).unwrap().as_str().to_string();
        let count = caps.get(2).unwrap().as_str().parse::<u32>().unwrap_or(1);
        if get_atom_mass(&el).is_some() { max_comp.insert(el.clone(), count); }
        present_elements.insert(el);
    }

    let mut elements_to_test: Vec<String> = max_comp.keys().cloned().collect();
    if elements_to_test.is_empty() { return false; }

    elements_to_test.sort_by(|a, b| {
        let mass_a = get_atom_mass(a).unwrap();
        let mass_b = get_atom_mass(b).unwrap();
        mass_a.partial_cmp(&mass_b).unwrap().reverse()
    });

    let element_masses: Vec<f64> = elements_to_test.iter().map(|el| get_atom_mass(el).unwrap()).collect();

    if target_mass <= 0.0 { return false; }

    let min_mass = target_mass * (1.0 - ppm_tol / 1_000_000.0);
    let max_mass = target_mass * (1.0 + ppm_tol / 1_000_000.0);

    let max_counts: Vec<u32> = elements_to_test.iter().enumerate().map(|(idx, el)| {
        let max_by_mass = (max_mass / element_masses[idx]) as u32;
        std::cmp::min(*max_comp.get(el).unwrap_or(&0), max_by_mass)
    }).collect();

    let mut max_remaining_masses = vec![0.0; elements_to_test.len() + 1];
    for i in (0..elements_to_test.len()).rev() {
        max_remaining_masses[i] = max_remaining_masses[i + 1] + (max_counts[i] as f64 * element_masses[i]);
    }

    let mut results = Vec::new();
    let mut composition = vec![0; elements_to_test.len()];
    generate_combinations(0, 0.0, &mut composition, &element_masses, &max_counts, min_mass, max_mass, &max_remaining_masses, &mut results);

    !results.is_empty()
}
