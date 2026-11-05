// src/mz_correction.rs
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::deletion_report::DeletionReport;
use crate::adduct_mass_calculator;
use crate::spectrum::Spectrum;

/// Étape autonome de l'orchestrateur (STEP 6.5) : Mass & Adduct Correction
///
/// Recalcule la masse exacte et corrige l'adduit si la masse expérimentale dévie trop.
pub fn mz_correction_processing(
    py: Python,
    spectrum_list: Vec<Spectrum>,
    deletion_report: &mut DeletionReport,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {
    
    let total_items = spectrum_list.len();
    if let Some(cb) = &total_items_callback {
        cb.call1(py, (total_items,))?;
    }
    if let Some(cb) = &prefix_callback {
        cb.call1(py, ("Mass & Adduct Correction",))?;
    }
    if let Some(cb) = &item_type_callback {
        cb.call1(py, ("spectra",))?;
    }
    if let Some(cb) = &progress_callback {
        cb.call1(py, (0,))?;
    }

    let (adduct_massdiff_pos, adduct_massdiff_neg) = {
        let state = crate::global_state::STATE.read().unwrap();
        (state.adduct_massdiff_dict_pos.clone(), state.adduct_massdiff_dict_neg.clone())
    };

    let batch_size = 1000;
    let mut all_processed = Vec::with_capacity(total_items);
    let mut total_processed = 0;
    
    let mut skipped_empty_formula = 0;
    let mut skipped_theorical_0 = 0;
    let mut skipped_mz_range = 0;
    
    let chunks: Vec<_> = spectrum_list.chunks(batch_size).collect();
    
    for chunk in chunks {
        py.allow_threads(|| {
            let processed_chunk: Vec<Result<Spectrum, &'static str>> = chunk.par_iter().map(|spec| {
                let mut new_spec = spec.clone();
                
                let formula = new_spec.metadata.get("FORMULA").cloned().unwrap_or_default();
                let precursor_type = new_spec.metadata.get("PRECURSORTYPE").cloned().unwrap_or_default();
                let ion_mode = new_spec.metadata.get("IONMODE").cloned().unwrap_or_default();
                
                let precursor_mz_str = new_spec.metadata.get("PRECURSORMZ").cloned().unwrap_or_default();
                let author_mz = precursor_mz_str.parse::<f64>().unwrap_or(0.0);
                
                let mut final_formula = String::new();
                let existing_precursor_formula = new_spec.metadata.get("PRECURSOR_FORMULA").cloned().unwrap_or_default();

                if !formula.is_empty() && !precursor_type.is_empty() && author_mz != 0.0 {
                    let precursor_formula = adduct_mass_calculator::apply_adduct_to_formula(&formula, &precursor_type);
                    let theorical_mz = adduct_mass_calculator::compute_mz(&precursor_formula).unwrap_or(0.0);
                    
                    if theorical_mz != 0.0 {
                        let delta = (author_mz - theorical_mz).abs();
                        
                        if delta <= 2.0 {
                            new_spec.metadata.insert("PRECURSOR_FORMULA".to_string(), precursor_formula.clone());
                            new_spec.metadata.insert("EXACT_MASS_RDKit".to_string(), theorical_mz.to_string());
                            final_formula = precursor_formula;
                        } else {
                            // Cherche un meilleur adduit
                            let adduct_dict = if ion_mode == "positive" { &adduct_massdiff_pos } else { &adduct_massdiff_neg };
                            
                            // On évalue tous les adduits candidats, puis on garde le plus proche (écart < 1.0).
                            // En cas d'égalité (ex. [M+Cl]- et [2M+2Cl]2- ont exactement le même m/z),
                            // on garde l'adduit le plus simple : voir best_candidate. Avant, l'égalité était
                            // tranchée par l'ordre de parcours de la HashMap, tiré au hasard à chaque lancement.
                            let mut candidates: Vec<(f64, &String, String, f64)> = Vec::new();
                            for candidate_adduct in adduct_dict.keys() {
                                let cand_formula = adduct_mass_calculator::apply_adduct_to_formula(&formula, candidate_adduct);
                                if cand_formula.is_empty() { continue; }

                                if let Some(cand_mz) = adduct_mass_calculator::compute_mz(&cand_formula) {
                                    let cand_delta = (author_mz - cand_mz).abs();
                                    candidates.push((cand_delta, candidate_adduct, cand_formula, cand_mz));
                                }
                            }
                            let best = best_candidate(&candidates, 1.0);
                            let best_adduct = best.map(|c| c.1.clone());
                            let (best_precursor_formula, best_mz) = best.map(|c| (c.2.clone(), c.3)).unwrap_or_default();

                            match best_adduct {
                                Some(new_adduct) => {
                                    new_spec.metadata.insert("PRECURSORTYPE".to_string(), new_adduct);
                                    new_spec.metadata.insert("PRECURSOR_FORMULA".to_string(), best_precursor_formula.clone());
                                    new_spec.metadata.insert("EXACT_MASS_RDKit".to_string(), best_mz.to_string());
                                    final_formula = best_precursor_formula;
                                },
                                None => {
                                    return Err("bad_adduct");
                                }
                            }
                        }
                    }
                }

                // Fallback: If we couldn't compute a final formula (because FORMULA is missing or adduct failed)
                // but we already have a PRECURSOR_FORMULA from a previous step, use it for MS2 verification!
                if final_formula.is_empty() && !existing_precursor_formula.is_empty() && existing_precursor_formula != "nan" {
                    final_formula = existing_precursor_formula;
                }
                
                // If STILL empty, it means we have no valid formula to test MS2 against. Skip crash test safely.
                if final_formula.is_empty() {
                    return Ok(new_spec.clone());
                }

                // ==========================================
                // ÉTAPE 5 : VÉRIFICATION MS2 & CRASH TEST
                // ==========================================
                
                let precursor_mz_val = new_spec.metadata.get("PRECURSORMZ").and_then(|v| v.parse::<f64>().ok());
                let mut z_prec = 1;
                if final_formula.contains("]") {
                    let suffix = final_formula.split("]").last().unwrap_or("");
                    let charge_str: String = suffix.chars().filter(|c| c.is_digit(10)).collect();
                    if !charge_str.is_empty() {
                        z_prec = charge_str.parse::<u32>().unwrap_or(1);
                    }
                }

                // Exclure les pics proches du précurseur (DELTA_PRECURSOR = 17.0 Th)
                // et extraire le nombre de décimales significatives
                let mut valid_peaks = Vec::new();
                for &(mz, int) in &new_spec.peaks {
                    if let Some(pmz) = precursor_mz_val {
                        if (mz - pmz).abs() <= 17.0 { continue; }
                    }
                    
                    let mz_str = format!("{:.6}", mz);
                    let mut dec = 0;
                    if let Some(idx) = mz_str.find('.') {
                        let dec_part = mz_str[idx+1..].trim_end_matches('0');
                        dec = dec_part.len();
                    }
                    valid_peaks.push((mz, int, dec));
                }
                
                // Trier par intensité décroissante
                valid_peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
                
                // Prendre les N_TOP = 2 plus intenses
                let top: Vec<&(f64, f64, usize)> = valid_peaks.iter().take(2).collect();
                
                if !top.is_empty() {
                    // Niveau 1 : vérification décimale (MIN_DECIMALS = 3).
                    // Pas pour le GC (règle unique is_gc_spectrum) : un spectre GC-EI est par nature en masse nominale.
                    // Le crash-test ci-dessous reste appliqué : il n'utilise que les pics à 3 décimales ou plus.
                    let is_gc = crate::globals_vars::is_gc_spectrum(
                        new_spec.metadata.get("FILENAME").map(|s| s.as_str()).unwrap_or(""),
                        new_spec.metadata.get("INSTRUMENTTYPE").map(|s| s.as_str()).unwrap_or(""),
                    );
                    if !is_gc && top.iter().all(|p| p.2 < 3) {
                        return Err("low_res");
                    }
                    
                    // Niveau 2 : Crash-Test Chimique (PPM_TOL = 20.0)
                    let precise: Vec<&&(f64, f64, usize)> = top.iter().filter(|p| p.2 >= 3).collect();
                    let mut has_valid_peak = false;
                    
                    if !final_formula.is_empty() {
                        for &&(frag_mz, _, _) in &precise {
                            for z_f in 1..=z_prec {
                                let target_mass = frag_mz * (z_f as f64) + (z_f as f64) * 0.000548579909;
                                if crate::de_novo_calculation::is_valid_peak(&final_formula, target_mass, 20.0) {
                                    has_valid_peak = true;
                                    break;
                                }
                            }
                            if has_valid_peak { break; }
                        }
                    }
                    
                    if !has_valid_peak && !precise.is_empty() && !final_formula.is_empty() {
                        return Err("crash_test");
                    }
                }

                Ok(new_spec)
            }).collect();
            
            for res in processed_chunk {
                match res {
                    Ok(valid_spec) => all_processed.push(valid_spec),
                    Err("bad_adduct") => deletion_report.no_or_bad_adduct += 1,
                    Err("low_res") => deletion_report.low_resolution_ms2 += 1,
                    Err("crash_test") => deletion_report.ms2_chemical_crash += 1,
                    Err("skip_formula") => skipped_empty_formula += 1,
                    Err("skip_mz") => skipped_theorical_0 += 1,
                    Err("skip_mz_range") => skipped_mz_range += 1,
                    _ => {}
                }
            }
        });
        
        total_processed += chunk.len();
        if let Some(cb) = &progress_callback {
            cb.call1(py, (total_processed,))?;
        }
    }
    

    Ok(all_processed)
}


/// Écart de m/z en dessous duquel deux adduits sont considérés à égalité (absorbe les arrondis
/// des calculs en virgule flottante ; très en dessous de toute précision d'instrument).
const ADDUCT_TIE_TOLERANCE: f64 = 1e-9;

/// Ordre de préférence entre adduits à égalité : le moins de molécules M, puis la charge la plus
/// faible, puis l'ordre alphabétique. Un adduit au format non reconnu passe après les autres.
fn adduct_simplicity(adduct: &str) -> (u32, u32, &str) {
    let (n, z) = adduct_mass_calculator::adduct_nmol_and_charge(adduct).unwrap_or((u32::MAX, u32::MAX));
    (n, z, adduct)
}

/// Choisit l'adduit le plus proche parmi `candidates` (écart, adduit, formule, m/z), à condition
/// que l'écart soit < `max_delta`. Parmi les adduits à égalité (à ADDUCT_TIE_TOLERANCE près),
/// le plus simple gagne. Le résultat ne dépend pas de l'ordre des candidats.
fn best_candidate<'a, 'b>(candidates: &'a [(f64, &'b String, String, f64)], max_delta: f64) -> Option<&'a (f64, &'b String, String, f64)> {
    let min_delta = candidates.iter().map(|c| c.0).filter(|d| !d.is_nan()).fold(f64::INFINITY, f64::min);
    if !(min_delta < max_delta) {
        return None;
    }
    candidates.iter()
        .filter(|c| c.0 < max_delta && c.0 - min_delta <= ADDUCT_TIE_TOLERANCE)
        .min_by(|a, b| adduct_simplicity(a.1).cmp(&adduct_simplicity(b.1)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cands(list: &[(f64, &str)]) -> Vec<(f64, String)> {
        list.iter().map(|(d, a)| (*d, a.to_string())).collect()
    }

    fn pick(list: &[(f64, &str)]) -> Option<String> {
        let owned = cands(list);
        let c: Vec<(f64, &String, String, f64)> = owned.iter().map(|(d, a)| (*d, a, String::new(), 0.0)).collect();
        best_candidate(&c, 1.0).map(|x| x.1.clone())
    }

    #[test]
    fn egalite_l_adduit_le_plus_simple_gagne() {
        assert_eq!(pick(&[(0.001, "[2M+2Cl]2-"), (0.001, "[M+Cl]-")]).as_deref(), Some("[M+Cl]-"));
        assert_eq!(pick(&[(0.001, "[M+Cl]-"), (0.001, "[2M+2Cl]2-")]).as_deref(), Some("[M+Cl]-"));
        // Écart minuscule dû aux arrondis : toujours considéré comme une égalité
        assert_eq!(pick(&[(0.001 + 1e-12, "[M+Cl]-"), (0.001, "[2M+2Cl]2-")]).as_deref(), Some("[M+Cl]-"));
    }

    #[test]
    fn le_plus_proche_gagne_toujours_sans_egalite() {
        assert_eq!(pick(&[(0.002, "[M+Cl]-"), (0.001, "[2M+2Cl]2-")]).as_deref(), Some("[2M+2Cl]2-"));
        assert_eq!(pick(&[(0.5, "[M-H]-"), (0.3, "[M+FA-H]-"), (0.9, "[M+Cl]-")]).as_deref(), Some("[M+FA-H]-"));
    }

    #[test]
    fn rien_sous_le_seuil_aucun_adduit() {
        assert_eq!(pick(&[(1.0, "[M-H]-"), (1.5, "[M+Cl]-")]), None);
        assert_eq!(pick(&[]), None);
    }

    #[test]
    fn resultat_independant_de_l_ordre() {
        let base = [(0.004, "[M-H]-"), (0.001, "[2M+2Cl]2-"), (0.001, "[M+Cl]-"), (0.001, "[3M+3Cl]3-"), (0.2, "[M+FA-H]-")];
        let attendu = pick(&base);
        let mut v = base.to_vec();
        for i in 0..v.len() {
            v.rotate_left(1);
            assert_eq!(pick(&v), attendu, "ordre {i}");
            v.reverse();
            assert_eq!(pick(&v), attendu, "ordre inversé {i}");
        }
        assert_eq!(attendu.as_deref(), Some("[M+Cl]-"));
    }
}
