// src/report.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use std::fs::File;
use std::io::Write;
use std::path::Path;

// Fonction utilitaire pour extraire la taille et les INCHIKEYs uniques d'un Vec<Spectrum>

fn get_df_stats(list: &Vec<Spectrum>) -> PyResult<(usize, usize)> {
    let length = list.len();
    let mut unique = 0;
    if length > 0 {
        let mut unique_set = std::collections::HashSet::new();
        for spec in list {
            if let Some(val) = spec.metadata.get("INCHIKEY") {
                let s = val.to_string();
                if !s.is_empty() && s.to_lowercase() != "nan" {
                    unique_set.insert(s);
                }
            }
        }
        unique = unique_set.len();
    }
    Ok((length, unique))
}


/// Génère le rapport final (report.txt) récapitulant les paramètres et le résultat du nettoyage.
///
/// Pour un développeur Python : Ce fichier contient le formatage d'une très longue chaîne.
/// En Rust, on utilise `format!(r#"..."#)` pour manipuler du texte brut multi-lignes ("raw strings") 
/// sans avoir besoin d'échapper les guillemets ou caractères spéciaux.
fn capitalize_words(s: &str) -> String {
    if s.to_uppercase() == "NOT FOUND" || s.to_uppercase() == "UNKNOWN" {
        return "NOT FOUND".to_string();
    }
    s.split_whitespace()
     .map(|word| {
         let mut c = word.chars();
         match c.next() {
             None => String::new(),
             Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
         }
     })
     .collect::<Vec<String>>()
     .join(" ")
}

pub fn generate_report_processing(
    _py: Python,
    output_directory: String,
    date_str: String,
    total_time_str: String,
    parameters_dict: &std::collections::HashMap<String, f64>,
    input_paths: &Vec<String>,
    deletion_report: &crate::deletion_report::DeletionReport,
    pos_lc: &Vec<Spectrum>,
    pos_lc_insilico: &Vec<Spectrum>,
    pos_gc: &Vec<Spectrum>,
    pos_gc_insilico: &Vec<Spectrum>,
    neg_lc: &Vec<Spectrum>,
    neg_lc_insilico: &Vec<Spectrum>,
    neg_gc: &Vec<Spectrum>,
    neg_gc_insilico: &Vec<Spectrum>,
) -> PyResult<()> {

    // 1. Extraction des statistiques des 8 DataFrames
    let (pos_lc_len, pos_lc_uniq) = get_df_stats(pos_lc)?;
    let (pos_lc_in_len, pos_lc_in_uniq) = get_df_stats(pos_lc_insilico)?;
    let (pos_gc_len, pos_gc_uniq) = get_df_stats(pos_gc)?;
    let (pos_gc_in_len, pos_gc_in_uniq) = get_df_stats(pos_gc_insilico)?;
    let (neg_lc_len, neg_lc_uniq) = get_df_stats(neg_lc)?;
    let (neg_lc_in_len, neg_lc_in_uniq) = get_df_stats(neg_lc_insilico)?;
    let (neg_gc_len, neg_gc_uniq) = get_df_stats(neg_gc)?;
    let (neg_gc_in_len, neg_gc_in_uniq) = get_df_stats(neg_gc_insilico)?;

    let total_spectra = pos_lc_len + pos_lc_in_len + pos_gc_len + pos_gc_in_len + neg_lc_len + neg_lc_in_len + neg_gc_len + neg_gc_in_len;
    let total_unique = pos_lc_uniq + pos_lc_in_uniq + pos_gc_uniq + pos_gc_in_uniq + neg_lc_uniq + neg_lc_in_uniq + neg_gc_uniq + neg_gc_in_uniq;

    // 2. Formatage des chaînes (Equivalent exact de votre Python)

    // Helpers natifs
    let get_bool = |key: &str| -> bool {
        parameters_dict.get(key).cloned().unwrap_or(0.0) == 1.0
    };
    let get_f64 = |key: &str| -> f64 {
        parameters_dict.get(key).cloned().unwrap_or(0.0)
    };

    
    
    let mut msp_files = String::new();
    let mut json_files = String::new();
    let mut csv_files = String::new();
    let mut mgf_files = String::new();

    for file_str in input_paths {
        let formatted = format!("<li>{}</li>", file_str);
        if file_str.ends_with(".msp") { msp_files.push_str(&formatted); }
        else if file_str.ends_with(".json") { json_files.push_str(&formatted); }
        else if file_str.ends_with(".csv") { csv_files.push_str(&formatted); }
        else if file_str.ends_with(".mgf") { mgf_files.push_str(&formatted); }
    }

    if msp_files.is_empty() { msp_files = "<li>-- no file --</li>".to_string(); }
    if json_files.is_empty() { json_files = "<li>-- no file --</li>".to_string(); }
    if csv_files.is_empty() { csv_files = "<li>-- no file --</li>".to_string(); }
    if mgf_files.is_empty() { mgf_files = "<li>-- no file --</li>".to_string(); }

    let get_bool_str = |key: &str| -> &str { if get_bool(key) { "ON" } else { "OFF" } };
    let get_class_str = |key: &str| -> &str { if get_bool(key) { "on" } else { "off" } };

    let total_deletions = deletion_report.duplicatas_removed +
        deletion_report.previously_cleaned +
        deletion_report.no_peaks_list +
        deletion_report.no_smiles_no_inchi_no_inchikey +
        deletion_report.no_precursor_mz +
        deletion_report.no_or_bad_adduct +
        deletion_report.low_entropy_score +
        deletion_report.minimum_peaks_not_requiered +
        deletion_report.all_peaks_above_precursor_mz +
        deletion_report.no_peaks_in_mz_range +
        deletion_report.minimum_high_peaks_not_requiered +
        deletion_report.low_resolution_ms2 +
        deletion_report.ms2_chemical_crash;

    let total_input = total_spectra + total_deletions;
    let max_h = if total_input > 0 { total_input as f64 } else { 1.0 };
    let pct_del = (total_deletions as f64 / max_h) * 100.0;
    let pct_out = (total_spectra as f64 / max_h) * 100.0;

    let calc_h = |val: usize| {
        if total_deletions > 0 { (val as f64 / total_deletions as f64) * 100.0 } else { 0.0 }
    };

    let parts: Vec<&str> = date_str.split("__").collect();
    let display_date = if parts.len() == 2 {
        let d = parts[0].replace("_", "-");
        let t = parts[1].replace("_", ":");
        format!("{} at {}", d, t)
    } else {
        date_str.replace("_", " ")
    };

    let mut html = String::from(include_str!("report_template.html"));
    html = html.replace("{TOTAL_TIME}", &total_time_str);

    html = html.replace("{DATE}", &display_date);
    html = html.replace("{MSP}", &msp_files);
    html = html.replace("{JSON}", &json_files);
    html = html.replace("{CSV}", &csv_files);
    html = html.replace("{MGF}", &mgf_files);
    html = html.replace("{OUT_DIR}", &output_directory);
    html = html.replace("{OUT_CSV}", if get_bool("csv") { "YES" } else { "NO" });
    html = html.replace("{OUT_MSP}", if get_bool("msp") { "YES" } else { "NO" });
    html = html.replace("{OUT_JSON}", if get_bool("json") { "YES" } else { "NO" });
    html = html.replace("{OUT_MZSPECLIB_JSON}", if get_bool("mzspeclib_json") { "YES" } else { "NO" });

    html = html.replace("{V1}", get_bool_str("normalize_intensity")).replace("{C1}", get_class_str("normalize_intensity"));
    html = html.replace("{V2}", get_bool_str("remove_peak_above_precursormz")).replace("{C2}", get_class_str("remove_peak_above_precursormz"));
    html = html.replace("{V3}", get_bool_str("check_minimum_peak_requiered")).replace("{C3}", get_class_str("check_minimum_peak_requiered"));
    html = html.replace("{V3_1}", &get_f64("check_minimum_peak_requiered_n_peaks").to_string());
    html = html.replace("{V4}", get_bool_str("reduce_peak_list")).replace("{C4}", get_class_str("reduce_peak_list"));
    html = html.replace("{V4_1}", &get_f64("reduce_peak_list_max_peaks").to_string());
    html = html.replace("{V5}", get_bool_str("remove_spectrum_under_entropy_score")).replace("{C5}", get_class_str("remove_spectrum_under_entropy_score"));
    html = html.replace("{V5_1}", &get_f64("remove_spectrum_under_entropy_score_value").to_string());
    html = html.replace("{V6}", get_bool_str("keep_mz_in_range")).replace("{C6}", get_class_str("keep_mz_in_range"));
    html = html.replace("{V6_1}", &get_f64("keep_mz_in_range_from_mz").to_string()).replace("{V6_2}", &get_f64("keep_mz_in_range_to_mz").to_string());
    html = html.replace("{V7}", get_bool_str("check_minimum_of_high_peaks_requiered")).replace("{C7}", get_class_str("check_minimum_of_high_peaks_requiered"));
    html = html.replace("{V7_1}", &get_f64("check_minimum_of_high_peaks_requiered_intensity_percent").to_string());
    html = html.replace("{V7_2}", &get_f64("check_minimum_of_high_peaks_requiered_no_peaks").to_string());
    html = html.replace("{V8}", get_bool_str("reset_updates")).replace("{C8}", get_class_str("reset_updates"));

    html = html.replace("{TOTAL_IN}", &total_input.to_string());
    html = html.replace("{TOTAL_DEL}", &total_deletions.to_string());
    html = html.replace("{TOTAL_OUT}", &total_spectra.to_string());
    html = html.replace("{TOTAL_UNIQ}", &total_unique.to_string());
    
    html = html.replace("{PCT_DEL}", &pct_del.to_string());
    html = html.replace("{PCT_OUT}", &pct_out.to_string());

    html = html.replace("{H_PEAKS}", &calc_h(deletion_report.no_peaks_list).to_string());
    html = html.replace("{V_PEAKS}", &deletion_report.no_peaks_list.to_string());
    html = html.replace("{H_SMILES}", &calc_h(deletion_report.no_smiles_no_inchi_no_inchikey).to_string());
    html = html.replace("{V_SMILES}", &deletion_report.no_smiles_no_inchi_no_inchikey.to_string());
    html = html.replace("{H_PRECURSOR}", &calc_h(deletion_report.no_precursor_mz).to_string());
    html = html.replace("{V_PRECURSOR}", &deletion_report.no_precursor_mz.to_string());
    html = html.replace("{H_ADDUCT}", &calc_h(deletion_report.no_or_bad_adduct).to_string());
    html = html.replace("{V_ADDUCT}", &deletion_report.no_or_bad_adduct.to_string());
    html = html.replace("{H_ENTROPY}", &calc_h(deletion_report.low_entropy_score).to_string());
    html = html.replace("{V_ENTROPY}", &deletion_report.low_entropy_score.to_string());
    html = html.replace("{H_MIN_PEAKS}", &calc_h(deletion_report.minimum_peaks_not_requiered).to_string());
    html = html.replace("{V_MIN_PEAKS}", &deletion_report.minimum_peaks_not_requiered.to_string());
    html = html.replace("{H_ABOVE_PREC}", &calc_h(deletion_report.all_peaks_above_precursor_mz).to_string());
    html = html.replace("{V_ABOVE_PREC}", &deletion_report.all_peaks_above_precursor_mz.to_string());
    html = html.replace("{H_RANGE}", &calc_h(deletion_report.no_peaks_in_mz_range).to_string());
    html = html.replace("{V_RANGE}", &deletion_report.no_peaks_in_mz_range.to_string());
    html = html.replace("{H_HIGH_PEAKS}", &calc_h(deletion_report.minimum_high_peaks_not_requiered).to_string());
    html = html.replace("{V_HIGH_PEAKS}", &deletion_report.minimum_high_peaks_not_requiered.to_string());
    html = html.replace("{H_DUP}", &calc_h(deletion_report.duplicatas_removed).to_string());
    html = html.replace("{V_DUP}", &deletion_report.duplicatas_removed.to_string());
    html = html.replace("{H_PREV}", &calc_h(deletion_report.previously_cleaned).to_string());
    html = html.replace("{V_PREV}", &deletion_report.previously_cleaned.to_string());
    html = html.replace("{H_LOW_RES}", &calc_h(deletion_report.low_resolution_ms2).to_string());
    html = html.replace("{V_LOW_RES}", &deletion_report.low_resolution_ms2.to_string());
    html = html.replace("{H_CRASH}", &calc_h(deletion_report.ms2_chemical_crash).to_string());
    html = html.replace("{V_CRASH}", &deletion_report.ms2_chemical_crash.to_string());


    let lc_pos_tot = pos_lc_len + pos_lc_in_len;
    let lc_neg_tot = neg_lc_len + neg_lc_in_len;
    let gc_pos_tot = pos_gc_len + pos_gc_in_len;
    let gc_neg_tot = neg_gc_len + neg_gc_in_len;
    let lc_tot = lc_pos_tot + lc_neg_tot;
    let gc_tot = gc_pos_tot + gc_neg_tot;

    let lc_pos_uniq = pos_lc_uniq + pos_lc_in_uniq;
    let lc_neg_uniq = neg_lc_uniq + neg_lc_in_uniq;
    let gc_pos_uniq = pos_gc_uniq + pos_gc_in_uniq;
    let gc_neg_uniq = neg_gc_uniq + neg_gc_in_uniq;
    let lc_uniq = lc_pos_uniq + lc_neg_uniq;
    let gc_uniq = gc_pos_uniq + gc_neg_uniq;

    html = html.replace("{LC_TOT}", &lc_tot.to_string());
    html = html.replace("{GC_TOT}", &gc_tot.to_string());
    html = html.replace("{LC_POS_TOT}", &lc_pos_tot.to_string());
    html = html.replace("{LC_NEG_TOT}", &lc_neg_tot.to_string());
    html = html.replace("{GC_POS_TOT}", &gc_pos_tot.to_string());
    html = html.replace("{GC_NEG_TOT}", &gc_neg_tot.to_string());

    html = html.replace("{LC_UNIQ}", &lc_uniq.to_string());
    html = html.replace("{GC_UNIQ}", &gc_uniq.to_string());
    html = html.replace("{LC_POS_UNIQ}", &lc_pos_uniq.to_string());
    html = html.replace("{LC_NEG_UNIQ}", &lc_neg_uniq.to_string());
    html = html.replace("{GC_POS_UNIQ}", &gc_pos_uniq.to_string());
    html = html.replace("{GC_NEG_UNIQ}", &gc_neg_uniq.to_string());

    html = html.replace("{T1}", &pos_lc_len.to_string()); html = html.replace("{U1}", &pos_lc_uniq.to_string());
    html = html.replace("{T2}", &neg_lc_len.to_string()); html = html.replace("{U2}", &neg_lc_uniq.to_string());
    html = html.replace("{T3}", &pos_lc_in_len.to_string()); html = html.replace("{U3}", &pos_lc_in_uniq.to_string());
    html = html.replace("{T4}", &neg_lc_in_len.to_string()); html = html.replace("{U4}", &neg_lc_in_uniq.to_string());
    html = html.replace("{T5}", &pos_gc_len.to_string()); html = html.replace("{U5}", &pos_gc_uniq.to_string());
    html = html.replace("{T6}", &neg_gc_len.to_string()); html = html.replace("{U6}", &neg_gc_uniq.to_string());
    html = html.replace("{T7}", &pos_gc_in_len.to_string()); html = html.replace("{U7}", &pos_gc_in_uniq.to_string());
    html = html.replace("{T8}", &neg_gc_in_len.to_string()); html = html.replace("{U8}", &neg_gc_in_uniq.to_string());

    
    let all_spectra = pos_lc.iter()
        .chain(pos_lc_insilico.iter())
        .chain(pos_gc.iter())
        .chain(pos_gc_insilico.iter())
        .chain(neg_lc.iter())
        .chain(neg_lc_insilico.iter())
        .chain(neg_gc.iter())
        .chain(neg_gc_insilico.iter());

    let mut cf_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut np_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut inst_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut db_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    let root_cf = "ClassyFire";
    let root_np = "NPClassifier";
    let root_inst = "Instrument";
    let root_db = "Databases";

    cf_counts.insert(root_cf.to_string(), 0);
    np_counts.insert(root_np.to_string(), 0);
    db_counts.insert(root_db.to_string(), 0);

    let mut cf_parents: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut cf_labels: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    let mut np_parents: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut np_labels: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    let mut db_parents: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut db_labels: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    cf_parents.insert(root_cf.to_string(), "".to_string());
    cf_labels.insert(root_cf.to_string(), root_cf.to_string());

    np_parents.insert(root_np.to_string(), "".to_string());
    np_labels.insert(root_np.to_string(), root_np.to_string());

    db_parents.insert(root_db.to_string(), "".to_string());
    db_labels.insert(root_db.to_string(), root_db.to_string());

    let clean_str = |s: &str| -> String {
        let t = s.trim().trim_matches('"').trim().to_string();
        if t.is_empty() || t.to_lowercase() == "nan" || t.to_lowercase() == "unknown" { "NOT FOUND".to_string() } else { t }
    };

    for spec in all_spectra {
        let cf_sup = clean_str(spec.metadata.get("CLASSYFIRE_SUPERCLASS").map(|s| s.as_str()).unwrap_or(""));
        let cf_cla = clean_str(spec.metadata.get("CLASSYFIRE_CLASS").map(|s| s.as_str()).unwrap_or(""));
        let cf_sub = clean_str(spec.metadata.get("CLASSYFIRE_SUBCLASS").map(|s| s.as_str()).unwrap_or(""));

        let p1 = format!("{}|{}", root_cf, cf_sup);
        let p2 = format!("{}|{}", p1, cf_cla);
        let p3 = format!("{}|{}", p2, cf_sub);

        cf_labels.insert(p1.clone(), cf_sup);
        cf_parents.insert(p1.clone(), root_cf.to_string());
        
        cf_labels.insert(p2.clone(), cf_cla);
        cf_parents.insert(p2.clone(), p1.clone());

        cf_labels.insert(p3.clone(), cf_sub);
        cf_parents.insert(p3.clone(), p2.clone());

        *cf_counts.entry(root_cf.to_string()).or_insert(0) += 1;
        *cf_counts.entry(p1).or_insert(0) += 1;
        *cf_counts.entry(p2).or_insert(0) += 1;
        *cf_counts.entry(p3).or_insert(0) += 1;

        let np_path = clean_str(spec.metadata.get("NPCLASS_PATHWAY").map(|s| s.as_str()).unwrap_or(""));
        let np_sup = clean_str(spec.metadata.get("NPCLASS_SUPERCLASS").map(|s| s.as_str()).unwrap_or(""));
        let np_cla = clean_str(spec.metadata.get("NPCLASS_CLASS").map(|s| s.as_str()).unwrap_or(""));

        let n1 = format!("{}|{}", root_np, np_path);
        let n2 = format!("{}|{}", n1, np_sup);
        let n3 = format!("{}|{}", n2, np_cla);

        np_labels.insert(n1.clone(), np_path);
        np_parents.insert(n1.clone(), root_np.to_string());

        np_labels.insert(n2.clone(), np_sup);
        np_parents.insert(n2.clone(), n1.clone());

        np_labels.insert(n3.clone(), np_cla);
        np_parents.insert(n3.clone(), n2.clone());

        *np_counts.entry(root_np.to_string()).or_insert(0) += 1;
        *np_counts.entry(n1).or_insert(0) += 1;
        *np_counts.entry(n2).or_insert(0) += 1;
        *np_counts.entry(n3).or_insert(0) += 1;

        let db_name = clean_str(spec.metadata.get("DATABASE_NAME").map(|s| s.as_str()).unwrap_or(""));
        let p1_db = format!("{}|{}", root_db, db_name);

        db_labels.insert(p1_db.clone(), db_name);
        db_parents.insert(p1_db.clone(), root_db.to_string());

        *db_counts.entry(root_db.to_string()).or_insert(0) += 1;
        *db_counts.entry(p1_db).or_insert(0) += 1;

        let t_marque = capitalize_words(spec.metadata.get("TREE_MARQUE").map(|s| s.as_str()).unwrap_or("NOT FOUND"));
        let t_modele = capitalize_words(spec.metadata.get("TREE_MODELE").map(|s| s.as_str()).unwrap_or("NOT FOUND"));
        let t_spec = capitalize_words(spec.metadata.get("TREE_SPECTYPE").map(|s| s.as_str()).unwrap_or("NOT FOUND"));
        let t_inst = capitalize_words(spec.metadata.get("TREE_INSTYPE").map(|s| s.as_str()).unwrap_or("NOT FOUND"));
        let t_ioni = capitalize_words(spec.metadata.get("TREE_IONI").map(|s| s.as_str()).unwrap_or("NOT FOUND"));

        let p1_inst = format!("{}|{}", root_inst, t_marque);
        let p2_inst = format!("{}|{}", p1_inst, t_modele);
        let p3_inst = format!("{}|{}", p2_inst, t_spec);
        let p4_inst = format!("{}|{}", p3_inst, t_inst);
        let p5_inst = format!("{}|{}", p4_inst, t_ioni);

        *inst_counts.entry(root_inst.to_string()).or_insert(0) += 1;
        *inst_counts.entry(p1_inst).or_insert(0) += 1;
        *inst_counts.entry(p2_inst).or_insert(0) += 1;
        *inst_counts.entry(p3_inst).or_insert(0) += 1;
        *inst_counts.entry(p4_inst).or_insert(0) += 1;
        *inst_counts.entry(p5_inst).or_insert(0) += 1;
    }

    let mut cf_ids_list = Vec::new();
    let mut cf_labels_list = Vec::new();
    let mut cf_parents_list = Vec::new();
    let mut cf_values_list = Vec::new();

    for (id, val) in &cf_counts {
        if *val == 0 { continue; }
        cf_ids_list.push(id.clone());
        cf_labels_list.push(cf_labels.get(id).cloned().unwrap_or_else(|| "".to_string()));
        cf_parents_list.push(cf_parents.get(id).cloned().unwrap_or_else(|| "".to_string()));
        cf_values_list.push(*val);
    }

    let mut np_ids_list = Vec::new();
    let mut np_labels_list = Vec::new();
    let mut np_parents_list = Vec::new();
    let mut np_values_list = Vec::new();

    for (id, val) in &np_counts {
        if *val == 0 { continue; }
        np_ids_list.push(id.clone());
        np_labels_list.push(np_labels.get(id).cloned().unwrap_or_else(|| "".to_string()));
        np_parents_list.push(np_parents.get(id).cloned().unwrap_or_else(|| "".to_string()));
        np_values_list.push(*val);
    }

    let mut inst_ids = Vec::new();
    let mut inst_labels = Vec::new();
    let mut inst_parents = Vec::new();
    let mut inst_values = Vec::new();

    for (id, val) in &inst_counts {
        if *val == 0 { continue; }
        inst_ids.push(id.clone());
        inst_values.push(*val);
        
        if id == root_inst {
            inst_labels.push(root_inst.to_string());
            inst_parents.push("".to_string());
        } else {
            let parts: Vec<&str> = id.split('|').collect();
            inst_labels.push(parts.last().unwrap().to_string());
            let parent_id = parts[..parts.len() - 1].join("|");
            inst_parents.push(parent_id);
        }
    }

    let mut db_ids_list = Vec::new();
    let mut db_labels_list = Vec::new();
    let mut db_parents_list = Vec::new();
    let mut db_values_list = Vec::new();

    for (id, val) in &db_counts {
        if *val == 0 { continue; }
        db_ids_list.push(id.clone());
        db_labels_list.push(db_labels.get(id).cloned().unwrap_or_else(|| "".to_string()));
        db_parents_list.push(db_parents.get(id).cloned().unwrap_or_else(|| "".to_string()));
        db_values_list.push(*val);
    }

    let cf_labels_json = serde_json::to_string(&cf_labels_list).unwrap_or_else(|_| "[]".to_string());
    let cf_parents_json = serde_json::to_string(&cf_parents_list).unwrap_or_else(|_| "[]".to_string());
    let cf_values_json = serde_json::to_string(&cf_values_list).unwrap_or_else(|_| "[]".to_string());
    let cf_ids_json = serde_json::to_string(&cf_ids_list).unwrap_or_else(|_| "[]".to_string());

    let np_labels_json = serde_json::to_string(&np_labels_list).unwrap_or_else(|_| "[]".to_string());
    let np_parents_json = serde_json::to_string(&np_parents_list).unwrap_or_else(|_| "[]".to_string());
    let np_values_json = serde_json::to_string(&np_values_list).unwrap_or_else(|_| "[]".to_string());
    let np_ids_json = serde_json::to_string(&np_ids_list).unwrap_or_else(|_| "[]".to_string());

    let db_labels_json = serde_json::to_string(&db_labels_list).unwrap_or_else(|_| "[]".to_string());
    let db_parents_json = serde_json::to_string(&db_parents_list).unwrap_or_else(|_| "[]".to_string());
    let db_values_json = serde_json::to_string(&db_values_list).unwrap_or_else(|_| "[]".to_string());
    let db_ids_json = serde_json::to_string(&db_ids_list).unwrap_or_else(|_| "[]".to_string());

    html = html.replace("{CF_LABELS}", &cf_labels_json);
    html = html.replace("{CF_PARENTS}", &cf_parents_json);
    html = html.replace("{CF_VALUES}", &cf_values_json);
    html = html.replace("{CF_IDS}", &cf_ids_json);

    html = html.replace("{NP_LABELS}", &np_labels_json);
    html = html.replace("{NP_PARENTS}", &np_parents_json);
    html = html.replace("{NP_VALUES}", &np_values_json);
    html = html.replace("{NP_IDS}", &np_ids_json);

    html = html.replace("{DB_LABELS}", &db_labels_json);
    html = html.replace("{DB_PARENTS}", &db_parents_json);
    html = html.replace("{DB_VALUES}", &db_values_json);
    html = html.replace("{DB_IDS}", &db_ids_json);

    let inst_labels_json = serde_json::to_string(&inst_labels).unwrap_or_else(|_| "[]".to_string());
    let inst_parents_json = serde_json::to_string(&inst_parents).unwrap_or_else(|_| "[]".to_string());
    let inst_values_json = serde_json::to_string(&inst_values).unwrap_or_else(|_| "[]".to_string());
    let inst_ids_json = serde_json::to_string(&inst_ids).unwrap_or_else(|_| "[]".to_string());

    html = html.replace("{INST_LABELS}", &inst_labels_json);
    html = html.replace("{INST_PARENTS}", &inst_parents_json);
    html = html.replace("{INST_VALUES}", &inst_values_json);
    html = html.replace("{INST_IDS}", &inst_ids_json);
    
    
    let mut combinations_map: std::collections::HashMap<Vec<String>, usize> = std::collections::HashMap::new();
    let mut spec_to_dbs: std::collections::HashMap<String, std::collections::HashSet<String>> = std::collections::HashMap::new();
    let mut db_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    let all_spectra_2 = pos_lc.iter()
        .chain(pos_lc_insilico.iter())
        .chain(pos_gc.iter())
        .chain(pos_gc_insilico.iter())
        .chain(neg_lc.iter())
        .chain(neg_lc_insilico.iter())
        .chain(neg_gc.iter())
        .chain(neg_gc_insilico.iter());

    for spec in all_spectra_2 {
        let splash = spec.metadata.get("SPLASH").map(|s| s.as_str()).unwrap_or("");
        let inchikey = spec.metadata.get("INCHIKEY").map(|s| s.as_str()).unwrap_or("");
        let dbname = spec.metadata.get("DATABASE_NAME").map(|s| s.as_str()).unwrap_or("Unknown").to_string();
        
        let id = format!("{}_{}", splash, inchikey);
        spec_to_dbs.entry(id).or_default().insert(dbname.clone());
        *db_counts.entry(dbname).or_insert(0) += 1;
    }

    let mut dbs_vec: Vec<(String, usize)> = db_counts.into_iter().collect();
    dbs_vec.sort_by(|a, b| b.1.cmp(&a.1));
    let top_dbs: std::collections::HashSet<String> = dbs_vec.into_iter().take(30).map(|(f, _)| f).collect();

    for (_, dbs) in spec_to_dbs {
        let mut sorted_dbs: Vec<String> = dbs.into_iter().filter(|f| top_dbs.contains(f)).collect();
        if sorted_dbs.is_empty() { continue; }
        sorted_dbs.sort();
        *combinations_map.entry(sorted_dbs).or_insert(0) += 1;
    }

    let mut upset_data = Vec::new();
    for (dbs, count) in combinations_map {
        upset_data.push(serde_json::json!({"sets": dbs, "count": count}));
    }
    let upset_json = serde_json::to_string(&upset_data).unwrap_or_else(|_| "[]".to_string());
    
    html = html.replace("{UPSET_DATA}", &upset_json);

    let final_cols = vec!["FILENAME", "FILEHASH", "PREDICTED", "SPLASH", "SPECTRUMID", "RESOLUTION", "SYNON", "IONIZATION", "MSLEVEL", "FRAGMENTATIONMETHODE", "NAME", "PRECURSORMZ", "EXACTMASS", "AVERAGEMASS", "PRECURSORTYPE", "INSTRUMENTTYPE", "INSTRUMENT", "SMILES", "INCHI", "INCHIKEY", "COLLISIONENERGY", "FORMULA", "RT", "IONMODE", "COMMENT", "ENTROPY", "CLASSYFIRE_SUPERCLASS", "CLASSYFIRE_CLASS", "CLASSYFIRE_SUBCLASS", "NPCLASS_PATHWAY", "NPCLASS_SUPERCLASS", "NPCLASS_CLASS", "NUM PEAKS", "PEAKS_LIST"];
    let mut col_counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    let mut total_spectra_count = 0;

    let all_spectra_3 = pos_lc.iter()
        .chain(pos_lc_insilico.iter())
        .chain(pos_gc.iter())
        .chain(pos_gc_insilico.iter())
        .chain(neg_lc.iter())
        .chain(neg_lc_insilico.iter())
        .chain(neg_gc.iter())
        .chain(neg_gc_insilico.iter());

    for spec in all_spectra_3 {
        total_spectra_count += 1;
        for &col in &final_cols {
            if col == "NUM PEAKS" || col == "PEAKS_LIST" {
                if !spec.peaks.is_empty() {
                    *col_counts.entry(col).or_insert(0) += 1;
                }
            } else if let Some(val) = spec.metadata.get(col) {
                let s = val.trim();
                let lower = s.to_lowercase();
                if !s.is_empty() && lower != "nan" && lower != "not found" && lower != "unknown" {
                    *col_counts.entry(col).or_insert(0) += 1;
                }
            }
        }
    }

    let mut col_stats_json = Vec::new();
    for &col in &final_cols {
        let count = *col_counts.get(col).unwrap_or(&0);
        let pct = if total_spectra_count > 0 { (count as f64 / total_spectra_count as f64) * 100.0 } else { 0.0 };
        col_stats_json.push(serde_json::json!({
            "column": col,
            "count": count,
            "percentage": pct
        }));
    }
    let col_stats_str = serde_json::to_string(&col_stats_json).unwrap_or_else(|_| "[]".to_string());
    html = html.replace("{COLUMNS_DATA}", &col_stats_str);

    let file_name = format!("report_{}.html", date_str);
    let report_path = Path::new(&output_directory).join(file_name);
    let mut file = File::create(report_path)?;
    file.write_all(html.as_bytes())?;
    
    Ok(())
}
