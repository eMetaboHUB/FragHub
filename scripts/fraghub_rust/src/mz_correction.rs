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
                            
                            let mut best_adduct: Option<String> = None;
                            let mut best_delta = 1.0;
                            let mut best_precursor_formula = String::new();
                            let mut best_mz = 0.0;
                            
                            for candidate_adduct in adduct_dict.keys() {
                                let cand_formula = adduct_mass_calculator::apply_adduct_to_formula(&formula, candidate_adduct);
                                if cand_formula.is_empty() { continue; }
                                
                                if let Some(cand_mz) = adduct_mass_calculator::compute_mz(&cand_formula) {
                                    let cand_delta = (author_mz - cand_mz).abs();
                                    if cand_delta < best_delta {
                                        best_delta = cand_delta;
                                        best_adduct = Some(candidate_adduct.clone());
                                        best_precursor_formula = cand_formula;
                                        best_mz = cand_mz;
                                    }
                                }
                            }
                            
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
                    // Niveau 1 : vérification décimale (MIN_DECIMALS = 3)
                    if top.iter().all(|p| p.2 < 3) {
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
