// src/convertors/msp_to_dict.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use crate::convertors::raw_metadata::RawMetadata;
use once_cell::sync::Lazy;
use regex::Regex;
use std::borrow::Cow;


/// Structure interne pour le multithreading MSP.
///
/// Pour un développeur Python : Remarquez le `Vec<(f64, f64)>` au lieu d'une liste de listes `[[f64, f64]]`.
/// En Rust, les tuples `(f64, f64)` sont stockés de manière contiguë (les uns à côté des autres)
/// en mémoire RAM, ce qui évite la fragmentation mémoire et accélère drastiquement les accès.
struct ParsedSpectrum {
    metadata: RawMetadata, // champs bruts, dans l'ordre du fichier
    peaks: Vec<(f64, f64)>, // <-- TUPLES : Fini la fragmentation mémoire !
}

/// Lit un nombre décimal qui peut utiliser la virgule (ex. "12,5").
///
/// Même résultat que `s.replace(",", ".").parse::<f64>().unwrap_or(0.0)`, mais sans créer
/// de nouvelle chaîne quand il n'y a pas de virgule (le cas de presque tous les pics).
fn parse_decimal(s: &str) -> f64 {
    if s.contains(',') {
        s.replace(",", ".").parse::<f64>().unwrap_or(0.0)
    } else {
        s.parse::<f64>().unwrap_or(0.0)
    }
}

/// Ligne de pics au format NIST : plusieurs pics séparés par « ; » (`85 999; 86 23; 87 5;`), chaque pic
/// pouvant porter une annotation entre guillemets (`85 999 "C5H9+";`).
static NIST_PEAK_LINE: Lazy<Regex> = Lazy::new(|| {
    let num = r"-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?";
    let pair = format!(r#"{num}(?:[ \t]+|:){num}(?:[ \t]+"[^"]*")?"#);
    Regex::new(&format!(r"^\s*(?:{pair}\s*;\s*)+(?:{pair}\s*;?\s*)?$")).unwrap()
});

/// Met un pic par ligne dans les lignes au format NIST (`85 999; 86 23;` -> `85 999` / `86 23`).
/// Avant, seul le premier pic de ces lignes était lu. Les autres lignes ne sont pas touchées : un fichier
/// sans « ; » dans ses pics est lu exactement comme avant.
fn split_nist_peak_lines(peaks_str: &str) -> Cow<'_, str> {
    if !peaks_str.contains(';') {
        return Cow::Borrowed(peaks_str);
    }
    let mut out = String::with_capacity(peaks_str.len() + 16);
    for (i, line) in peaks_str.split('\n').enumerate() {
        if i > 0 { out.push('\n'); }
        if NIST_PEAK_LINE.is_match(line) {
            let parts: Vec<&str> = line.split(';').map(|p| p.trim()).filter(|p| !p.is_empty()).collect();
            out.push_str(&parts.join("\n"));
        } else {
            out.push_str(line);
        }
    }
    Cow::Owned(out)
}

/// Fonction principale de parsing MSP via Rayon.
///
/// Pour un développeur Python : Le `filter_map` combiné au `par_iter` agit comme une compréhension 
/// de liste multi-threadée qui écarte automatiquement les valeurs `None`. Cela évite d'avoir à faire 
/// un filtrage manuel des "mauvais" spectres à la fin.
///
/// # Arguments
/// * `py` (Python) : Le token PyO3.
/// * `rust_strings` (Vec<String>) : Les spectres MSP bruts.
/// * `keys_dict` (HashMap<String, String>) : Mapping des clés.
/// * `keys_list` (Vec<String>) : Liste des clés à conserver.
/// * `progress_callback` (Option<PyObject>) : Callback de progression.
/// * `total_items_callback` (Option<PyObject>) : Callback de total.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type.
///
/// # Returns
/// * `PyResult<Vec<Spectrum>>` : La liste de spectres parsés.
pub fn msp_to_dict_processing(
    py: Python,
    rust_strings: Vec<String>,
    keys_dict: HashMap<String, String>,
    keys_list: Vec<String>,
    db_name: String,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<Spectrum>> {



    let total = rust_strings.len();
    if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total, 0)); }
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, ("Parsing MSP spectra:",)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("spectra",)); }

    let mut result_list = Vec::new();
    let mut processed = 0;
    let chunk_size = 2000;
    // Ensemble des clés à garder : même réponse que keys_list.contains(), sans parcourir la liste à chaque champ
    let keys_set: HashSet<&str> = keys_list.iter().map(|s| s.as_str()).collect();

    for chunk in rust_strings.chunks(chunk_size) {

        let parsed_chunk: Vec<ParsedSpectrum> = py.allow_threads(|| {
            chunk.par_iter().filter_map(|spectrum| {
                if let Some(caps) = crate::globals_vars::METADATA_PEAK_LIST_SPLIT_PATTERN_MSP.captures(spectrum) {
                    let metadata_str = caps.get(1).map_or("", |m| m.as_str());
                    let peaks_str = caps.get(2).map_or("", |m| m.as_str());

                    if !metadata_str.is_empty() && !peaks_str.is_empty() {
                        let mut metadata_matches = Vec::new();
                        for m in crate::globals_vars::METADATA_PATTERN_MSP.captures_iter(metadata_str) {
                            metadata_matches.push((m.get(1).unwrap().as_str().to_string(), m.get(2).unwrap().as_str().to_string()));
                        }

                        let mut new_matches = Vec::new();
                        let mut valid_comments = true;
                        for (k, v) in metadata_matches {
                            if crate::globals_vars::COMMENT_PATTERN.is_match(&k) {
                                new_matches.push((k.clone(), v.clone()));
                                if v.contains("=") {
                                    let mut found_sub = false;
                                    for sub_m in crate::globals_vars::SUB_FIELDS_PATTERN.captures_iter(&v) {
                                        let mut g_k = String::new();
                                        let mut g_v = String::new();
                                        let mut count = 0;
                                        for i in 1..sub_m.len() {
                                            if let Some(grp) = sub_m.get(i) {
                                                if !grp.as_str().is_empty() {
                                                    if count == 0 { g_k = grp.as_str().to_string(); count += 1; }
                                                    else if count == 1 { g_v = grp.as_str().to_string(); count += 1; break; }
                                                }
                                            }
                                        }
                                        if count == 2 {
                                            new_matches.push((g_k, g_v));
                                            found_sub = true;
                                        }
                                    }
                                    if !found_sub { valid_comments = false; break; }
                                }
                            } else { new_matches.push((k, v)); }
                        }

                        if valid_comments && !new_matches.is_empty() {
                            let mut metadata_dict = RawMetadata::default();
                            for (k_raw, v_raw) in new_matches.into_iter().filter(|(k, _)| !crate::globals_vars::COMPUTED_PATTERN.is_match(k)) {
                                let k = crate::globals_vars::METADATA_FIELDS_NAME_PATTERN.replace_all(&k_raw, "").to_lowercase().trim().to_string();
                                let v = crate::globals_vars::METADATA_STRIP_VALUE_PATTERN.replace_all(&v_raw, "").to_string();
                                metadata_dict.insert(k, v);
                            }

                            let mut peaks_array = Vec::new();
                            let peaks_str = split_nist_peak_lines(peaks_str);
                            for p in crate::globals_vars::PEAK_LIST_SPLIT_PATTERN.captures_iter(&peaks_str) {
                                let mz = parse_decimal(p.get(1).unwrap().as_str());
                                let int = parse_decimal(p.get(2).unwrap().as_str());
                                peaks_array.push((mz, int)); // ✅ TUPLE AU LIEU D'UN VEC
                            }

                            if !peaks_array.is_empty() {
                                return Some(ParsedSpectrum { metadata: metadata_dict, peaks: peaks_array });
                            }
                        }
                    }
                }
                None
            }).collect()
        });

        for parsed in parsed_chunk {
            let mut spec = Spectrum::default();
            // Renommage vers les clés FragHub : si deux champs donnent la même clé, le premier du fichier gagne
            parsed.metadata.map_into(&keys_dict, |mapped| keys_set.contains(mapped), &mut spec.metadata);
            spec.metadata.insert("DATABASE_NAME".to_string(), db_name.clone());

            spec.peaks = parsed.peaks;

            for key in &keys_list {
                if !spec.metadata.contains_key(key) && key != "PEAKS_LIST" { 
                    spec.metadata.insert(key.clone(), "".to_string()); 
                }
            }

            result_list.push(spec);
        }

        processed += chunk.len();
        if let Some(cb) = &progress_callback { let _ = cb.call1(py, (processed,)); }
        py.allow_threads(|| { std::thread::sleep(std::time::Duration::from_millis(1)); });
    }

    Ok(result_list)
}

#[cfg(test)]
mod tests {
    use super::parse_decimal;
    use super::split_nist_peak_lines;

    fn peaks(s: &str) -> Vec<(String, String)> {
        let s = split_nist_peak_lines(s);
        crate::globals_vars::PEAK_LIST_SPLIT_PATTERN.captures_iter(&s)
            .map(|c| (c[1].to_string(), c[2].to_string())).collect()
    }

    #[test]
    fn format_nist_plusieurs_pics_par_ligne() {
        let p = peaks("85 999; 86 23; 87 5;\n88 1; 89 2\n");
        assert_eq!(p.len(), 5);
        assert_eq!(p[2], ("87".to_string(), "5".to_string()));
        assert_eq!(p[4], ("89".to_string(), "2".to_string()));
        // annotations entre guillemets et séparateur « : »
        let p = peaks("41 47 \"C3H5+\"; 43:100 \"C2H3O+\";\n");
        assert_eq!(p, vec![("41".into(), "47".into()), ("43".into(), "100".into())]);
    }

    #[test]
    fn spectre_nist_complet() {
        // Bloc MSP NIST complet : la séparation métadonnées / pics garde bien toutes les lignes de pics
        let spec = "FILENAME: nist.msp\nName: Caffeine\nPrecursor_type: [M+H]+\nPrecursorMZ: 195.0877\nNum Peaks: 5\n42 12; 69 30; 110 45;\n138 999; 195 120;\n";
        let caps = crate::globals_vars::METADATA_PEAK_LIST_SPLIT_PATTERN_MSP.captures(spec).unwrap();
        assert!(caps.get(1).unwrap().as_str().contains("Num Peaks: 5"));
        assert_eq!(peaks(caps.get(2).unwrap().as_str()).len(), 5);
    }

    #[test]
    fn autres_formats_lus_comme_avant() {
        for s in ["85 999\n86 23\n", "85\t999\n86\t23", "85:999\n86:23\n", "85 999 \"a;b\"\n86 23\n", "85 999;\n86 23\n", "1.5e2 3,5\n"] {
            let avant: Vec<(String, String)> = crate::globals_vars::PEAK_LIST_SPLIT_PATTERN.captures_iter(s)
                .map(|c| (c[1].to_string(), c[2].to_string())).collect();
            assert_eq!(peaks(s), avant, "{s:?}");
        }
        // Sans « ; », la chaîne n'est même pas recopiée
        assert!(matches!(split_nist_peak_lines("85 999\n"), std::borrow::Cow::Borrowed(_)));
    }

    #[test]
    fn parse_decimal_identique_a_l_ancienne_ecriture() {
        let cas = ["123.4567", "12,5", "1,2,3", "-0.5", "1e-3", "1E+5", "7", "", "abc", "1.", ",5", "  1.0", "NaN", "inf"];
        for s in cas {
            let ancien = s.replace(",", ".").parse::<f64>().unwrap_or(0.0);
            let nouveau = parse_decimal(s);
            // to_bits : compare aussi les NaN et le signe de zéro
            assert_eq!(ancien.to_bits(), nouveau.to_bits(), "valeur {s:?}");
        }
    }
}
