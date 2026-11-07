// src/convertors/json_to_dict.rs
use pyo3::prelude::*;
use crate::spectrum::Spectrum;
use rayon::prelude::*;
use std::collections::HashMap;
use crate::convertors::raw_metadata::RawMetadata;
use std::time::Duration;


struct ParsedJsonSpectrum {
    metadata: RawMetadata, // champs bruts, dans l'ordre de l'objet JSON
    peaks: Vec<(f64, f64)>,
}

/// Parseur natif et ultra-rapide de listes de pics JSON.
///
/// Pour un développeur Python : Observez le mot-clé `unsafe`. En Rust, lire directement
/// la RAM en ignorant la vérification UTF-8 peut faire crasher le programme si les données 
/// sont invalides. Ici, on utilise `unsafe` car on a déjà vérifié à la main (via `is_ascii_digit()`) 
/// que les octets étaient valides. C'est ce type d'astuce qui rend Rust aussi rapide que le C ou le C++.
///
/// # Arguments
/// * `peak_list_string` (&str) : La liste de pics sous forme de texte brut.
///
/// # Returns
/// * `Vec<(f64, f64)>` : La liste des tuples `(m/z, intensité)`.
fn parse_json_peak_list_rust(peak_list_string: &str) -> Vec<(f64, f64)> {
    let bytes = peak_list_string.as_bytes();
    let mut numbers = Vec::with_capacity(1024);
    let mut start = None;

    for (i, &b) in bytes.iter().enumerate() {
        let is_num_char = b.is_ascii_digit() || b == b'.' || b == b'-' || b == b'e' || b == b'E' || b == b'+';
        if is_num_char {
            if start.is_none() { start = Some(i); }
        } else if let Some(s) = start {
            // SAFETY: Les octets sont garantis être des chiffres/symboles ASCII valides
            let s_val = unsafe { std::str::from_utf8_unchecked(&bytes[s..i]) };
            if let Ok(num) = s_val.parse::<f64>() { numbers.push(num); }
            start = None;
        }
    }
    if let Some(s) = start {
        let s_val = unsafe { std::str::from_utf8_unchecked(&bytes[s..]) };
        if let Ok(num) = s_val.parse::<f64>() { numbers.push(num); }
    }

    let mut peaks = Vec::with_capacity(numbers.len() / 2);
    for chunk in numbers.chunks_exact(2) {
        peaks.push((chunk[0], chunk[1]));
    }
    peaks
}

use crate::loaders::{INJECTED_FILENAME, INJECTED_FILEHASH, INJECTED_PREFIX};

/// Tableau JSON non vide de nombres ou de textes -> "v1, v2, ..." ; `None` pour tout autre tableau ou valeur.
fn scalar_list(val: &serde_json::Value) -> Option<String> {
    let arr = val.as_array()?;
    if arr.is_empty() {
        return None;
    }
    let mut parts = Vec::with_capacity(arr.len());
    for x in arr {
        if let Some(s) = x.as_str() {
            parts.push(s.to_string());
        } else if x.is_number() {
            parts.push(x.to_string());
        } else {
            return None;
        }
    }
    Some(parts.join(", "))
}

/// Lit une liste MoNA [{name, value, computed}] : valeur texte ou nombre. Les entrées « computed: true »
/// sont mises dans `computed` (valeurs de repli), les autres dans `metadata`.
fn read_mona_metadata(items: &[serde_json::Value], metadata: &mut RawMetadata, computed: &mut Vec<(String, String)>) {
    for item in items {
        if let Some(name) = item["name"].as_str() {
            let name_lower = name.to_lowercase();
            let val = if let Some(val) = item["value"].as_str() {
                val.to_string()
            } else if let Some(val) = item["value"].as_f64() {
                val.to_string()
            } else {
                continue;
            };
            if item["computed"].as_bool() == Some(true) {
                computed.push((name_lower, val));
            } else {
                metadata.insert(name_lower, val);
            }
        }
    }
}

fn parse_json_peak_array(val: &serde_json::Value) -> Vec<(f64, f64)> {
    let mut peaks = Vec::new();
    if let Some(arr) = val.as_array() {
        peaks.reserve(arr.len());
        for item in arr {
            if let Some(pair) = item.as_array() {
                if pair.len() >= 2 {
                    let mz = pair[0].as_f64().unwrap_or(0.0);
                    let int = pair[1].as_f64().unwrap_or(0.0);
                    peaks.push((mz, int));
                }
            }
        }
    }
    peaks
}

/// Analyse UN spectre JSON (texte venant du lecteur, avec le nom et l'empreinte du fichier en tête) :
/// champs bruts (avant renommage) et pics. `None` si le texte n'est pas un spectre lisible.
/// Sortie de la boucle parallèle telle quelle, pour pouvoir la tester.
fn parse_json_spectrum(json_str: &str) -> Option<ParsedJsonSpectrum> {
    let v: serde_json::Value = serde_json::from_str(&json_str).ok()?;
    let mut metadata = RawMetadata::with_capacity(32);
    let mut peaks = Vec::new();

    // Nom et empreinte du fichier : clés ajoutées par le lecteur (loaders::INJECTED_*). Avant, le
    // lecteur ajoutait "filename" / "filehash" : une clé "filename" du fichier les remplaçait.
    let is_mona = v.get("compound").is_some() && v.get("id").is_some() && v.get("metaData").is_some() && v.get("spectrum").is_some() && v.get(INJECTED_FILENAME).is_some();

    if let Some(fh) = v.get(INJECTED_FILEHASH).and_then(|h| h.as_str()) { metadata.insert("filehash".to_string(), fh.to_string()); }
    if let Some(fnm) = v.get(INJECTED_FILENAME).and_then(|f| f.as_str()) { metadata.insert("filename".to_string(), fnm.to_string()); }

    if is_mona {
        if let Some(comp_name) = v["compound"].get(0).and_then(|c| c["names"].get(0)).and_then(|n| n["name"].as_str()) {
            metadata.insert("compound_name".to_string(), comp_name.to_string());
        } else {
            metadata.insert("compound_name".to_string(), "".to_string());
        }

        // Valeurs « computed » de MoNA (recalculées par MoNA : SMILES, InChI, formule, masses...) :
        // gardées de côté et utilisées seulement à défaut de la valeur d'origine (voir plus bas).
        // Avant, un SMILES « computed » (sans stéréochimie) remplaçait celui du déposant.
        let mut computed: Vec<(String, String)> = Vec::new();

        if let Some(meta) = v["compound"].get(0).and_then(|c| c["metaData"].as_array()) {
            read_mona_metadata(meta, &mut metadata, &mut computed);
        }

        // Extraction de la classification (MoNA utilise ClassyFire par défaut)
        if let Some(classifications) = v["compound"].get(0).and_then(|c| c["classification"].as_array()) {
            for item in classifications {
                if let Some(name) = item["name"].as_str() {
                    let name_lower = format!("classyfire_{}", name.to_lowercase().replace(" ", "_"));
                    if let Some(val) = item["value"].as_str() {
                        metadata.insert(name_lower, val.to_string());
                    }
                }
            }
        }

        if let Some(id) = v["id"].as_str() { metadata.insert("spectrum_id".to_string(), id.to_string()); }

        if let Some(meta) = v["metaData"].as_array() {
            read_mona_metadata(meta, &mut metadata, &mut computed);
        }

        metadata.insert("predicted".to_string(), "false".to_string());
        if let Some(tags) = v["tags"].as_array() {
            for tag in tags {
                if tag["text"].as_str() == Some("In-Silico") {
                    metadata.insert("predicted".to_string(), "true".to_string());
                    break;
                }
            }
        }

        if let Some(splash) = v.get("splash").and_then(|s| s.get("splash")).and_then(|s| s.as_str()) {
            metadata.insert("splash".to_string(), splash.to_string());
        }

        if let Some(spec) = v["spectrum"].as_str() {
            peaks = parse_json_peak_list_rust(spec);
        }

        // Valeurs « computed » en dernier : seulement là où il n'y a pas de valeur d'origine
        // (ni même champ, ni synonyme rempli avant : elles passent après tous les autres champs)
        for (name, val) in computed {
            metadata.insert_if_absent(name, val);
        }
    } else {
        let mut matched = false;
        for pk in ["spectrum", "peaks_json", "peaks"] {
            if let Some(pval) = v.get(pk) {
                if let Some(pstr) = pval.as_str() {
                    peaks = parse_json_peak_list_rust(pstr);
                } else {
                    peaks = parse_json_peak_array(pval);
                }
                matched = true;
                break;
            }
        }
        if !matched { return None; }

        if let Some(obj) = v.as_object() {
            for (k, val) in obj {
                let k_str = k.as_str();
                if k_str != "spectrum" && k_str != "peaks_json" && k_str != "peaks" && !k_str.starts_with(INJECTED_PREFIX) {
                    if let Some(s) = val.as_str() {
                        metadata.insert(k.to_lowercase(), s.to_string());
                    } else if val.is_number() || val.is_boolean() {
                        metadata.insert(k.to_lowercase(), val.to_string());
                    } else if let Some(text) = scalar_list(val) {
                        // Liste de valeurs simples (MSnLib : "collision_energy": [20.0, 30.0, 60.0, 40.0] pour
                        // un spectre fusionné) : gardée en texte "20.0, 30.0, 60.0, 40.0". Avant, ignorée.
                        metadata.insert(k.to_lowercase(), text);
                    }
                }
            }
        }
        
        // Extraction flexible de CollisionEnergy du nom pour GNPS via Regex compilée
        if let Some(name) = metadata.get("compound_name").cloned() {
            if let Some(caps) = crate::globals_vars::GNPS_COLLISION_ENERGY_EXTRACT_PATTERN.captures(&name) {
                if let Some(energy_match) = caps.get(1) {
                    let energy = energy_match.as_str().to_string();
                    let real_name = crate::globals_vars::GNPS_COLLISION_ENERGY_EXTRACT_PATTERN.replace(&name, " ").to_string();
                    metadata.insert("compound_name".to_string(), real_name.trim().to_string());
                    metadata.insert("collision_energy".to_string(), energy.trim().to_string());
                }
            }
        }
    }

    Some(ParsedJsonSpectrum { metadata, peaks })
}

/// Fonction principale de parsing JSON via Rayon.
///
/// Pour un développeur Python : C'est ici que l'on remplace le `multiprocessing` de Python.
/// `.into_par_iter()` (fourni par la librairie Rayon) transforme instantanément notre boucle for
/// en une boucle multithreadée. Pas besoin de gérer les Process, les Queues ou de `pickle` les données :
/// Rayon s'occupe de répartir la charge sur tous les cœurs du processeur sans aucune fuite de mémoire.
///
/// # Arguments
/// * `py` (Python) : Le token PyO3 pour appeler les callbacks.
/// * `rust_strings` (Vec<String>) : Les spectres JSON sous forme de chaîne de caractères.
/// * `keys_dict` (HashMap<String, String>) : Le mapping des clés.
/// * `keys_list` (Vec<String>) : La liste des colonnes requises.
/// * `progress_callback` (Option<PyObject>) : Callback de progression.
/// * `total_items_callback` (Option<PyObject>) : Callback du nombre total.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe texte.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément.
///
/// # Returns
/// * `PyResult<Vec<Spectrum>>` : La liste de spectres parsée en parallèle.
pub fn json_to_dict_processing(
    py: Python,
    mut rust_strings: Vec<String>,
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
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, ("Parsing JSON spectrums:",)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("spectra",)); }

    let mut result_list = Vec::new();

    let mut processed = 0;
    let chunk_size = 3000;

    // --- BOUCLIER RAM 2 : Extraction rapide par la fin (O(1) déplacement) ---
    while !rust_strings.is_empty() {
        let end = rust_strings.len();
        let start = if end > chunk_size { end - chunk_size } else { 0 };
        let chunk: Vec<String> = rust_strings.drain(start..end).collect();
        let current_chunk_len = chunk.len();

        let parsed_chunk: Vec<ParsedJsonSpectrum> = py.allow_threads(|| {
            chunk.into_par_iter().filter_map(|json_str| parse_json_spectrum(&json_str)).collect()
        });

        for parsed in parsed_chunk {
            let mut spec = Spectrum::default();

            // Renommage vers les clés FragHub : si deux champs donnent la même clé, le premier du fichier gagne
            parsed.metadata.map_into(&keys_dict, |mapped| keys_list.iter().any(|k| k == mapped), &mut spec.metadata);
            spec.metadata.insert("DATABASE_NAME".to_string(), db_name.clone());

            spec.peaks = parsed.peaks;

            for key in &keys_list {
                if !spec.metadata.contains_key(key) && key != "PEAKS_LIST" {
                    spec.metadata.insert(key.clone(), "".to_string());
                }
            }

            result_list.push(spec);
        }

        processed += current_chunk_len;
        if let Some(cb) = &progress_callback { let _ = cb.call1(py, (processed,)); }

        py.allow_threads(|| { std::thread::sleep(Duration::from_millis(1)); });
    }

    Ok(result_list)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Texte tel que le donne le lecteur : nom et empreinte du fichier ajoutés en tête de l'objet
    fn as_loaded(obj: &str, filename: &str) -> String {
        format!(r#"{{"{}":"{}","{}":"h","#, INJECTED_FILENAME, filename, INJECTED_FILEHASH) + &obj.trim()[1..]
    }

    fn mapped(obj: &str, filename: &str, keys: &[(&str, &str)]) -> HashMap<String, String> {
        let k: HashMap<String, String> = keys.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect();
        let parsed = parse_json_spectrum(&as_loaded(obj, filename)).expect("spectre lu");
        let mut out = HashMap::new();
        parsed.metadata.map_into(&k, |_| true, &mut out);
        out
    }

    // Spectre MoNA réduit (structure réelle de l'export MoNA, valeurs du spectre WA002994)
    const MONA: &str = r#"{"compound":[{"names":[{"name":"Cefuroxime"}],
        "metaData":[
          {"name":"SMILES","value":"CO/N=C(/C1=CC=CO1)\\C(=O)N","computed":false},
          {"name":"molecular formula","value":"C16H16N4O8S","computed":false},
          {"name":"molecular formula","value":"C16H16N4O8S","computed":true},
          {"name":"total exact mass","value":424.068884472,"computed":true},
          {"name":"SMILES","value":"O=C(O)C1=C(COC(=N)O)","computed":true},
          {"name":"InChI","value":"InChI=1S/C16H16N4O8S/calcule","computed":true}]}],
        "id":"WA002994",
        "metaData":[{"name":"exact mass","value":"424.06888","computed":false},
                    {"name":"spectral entropy","value":4.12,"computed":true}],
        "tags":[{"text":"MassBank"}],
        "spectrum":"100:3.1 101:1.2"}"#;

    #[test]
    fn mona_valeur_d_origine_gardee_computed_en_repli() {
        let out = mapped(MONA, "mona.json", &[("smiles", "SMILES"), ("inchi", "INCHI"), ("molecular formula", "FORMULA"),
            ("total exact mass", "EXACTMASS"), ("exact mass", "EXACTMASS"), ("spectral entropy", "ENTROPY"),
            ("spectrum_id", "SPECTRUMID"), ("filename", "FILENAME")]);
        // SMILES du déposant, pas celui recalculé par MoNA
        assert_eq!(out["SMILES"], "CO/N=C(/C1=CC=CO1)\\C(=O)N");
        // Pas d'InChI d'origine : celui de MoNA sert de repli
        assert_eq!(out["INCHI"], "InChI=1S/C16H16N4O8S/calcule");
        // Masse exacte : celle d'origine (« exact mass ») passe avant la valeur calculée (« total exact mass »)
        assert_eq!(out["EXACTMASS"], "424.06888");
        assert_eq!(out["FORMULA"], "C16H16N4O8S");
        assert_eq!(out["ENTROPY"], "4.12");
        assert_eq!(out["SPECTRUMID"], "WA002994");
        assert_eq!(out["FILENAME"], "mona.json");
        assert_eq!(parse_json_spectrum(&as_loaded(MONA, "mona.json")).unwrap().peaks, vec![(100.0, 3.1), (101.0, 1.2)]);
    }

    #[test]
    fn json_generique_nom_de_fichier_d_origine() {
        let obj = r#"{"filename":"acquisition.mzML","compound_name":"X","peaks":[[100.0,1.0]]}"#;
        let out = mapped(obj, "entree.json", &[("filename", "FILENAME"), ("filehash", "FILEHASH"), ("compound_name", "NAME")]);
        // Nom du fichier de la base d'origine ; fichier d'entrée gardé en interne ; empreinte du fichier d'entrée
        assert_eq!(out["FILENAME"], "acquisition.mzML");
        assert_eq!(out[crate::globals_vars::INPUT_FILENAME_KEY], "entree.json");
        assert_eq!(out["FILEHASH"], "h");
        assert_eq!(out["NAME"], "X");
        // Sans nom d'origine : nom du fichier d'entrée
        let out = mapped(r#"{"compound_name":"X","peaks":[[100.0,1.0]]}"#, "entree.json", &[("filename", "FILENAME")]);
        assert_eq!(out["FILENAME"], "entree.json");
        // Les clés ajoutées par le lecteur ne deviennent pas des champs bruts
        let parsed = parse_json_spectrum(&as_loaded(obj, "entree.json")).unwrap();
        assert!(parsed.metadata.get(INJECTED_FILENAME).is_none());
    }

    #[test]
    fn json_liste_d_energies_lue_en_texte() {
        let obj = r#"{"compound_name":"X","collision_energy":[20.0, 30.0, 60.0, 40.0],"scan_number":[176, 177],"peaks":[[100.0,1.0]],"annotation_history":[{"a":1}]}"#;
        let out = mapped(obj, "msnlib.json", &[("collision_energy", "COLLISIONENERGY"), ("compound_name", "NAME")]);
        assert_eq!(out["COLLISIONENERGY"], "20.0, 30.0, 60.0, 40.0");
        // Énergie seule : inchangée
        let out = mapped(r#"{"collision_energy":60.0,"peaks":[[100.0,1.0]]}"#, "a.json", &[("collision_energy", "COLLISIONENERGY")]);
        assert_eq!(out["COLLISIONENERGY"], "60.0");
        // Liste d'objets : toujours ignorée
        let parsed = parse_json_spectrum(&as_loaded(obj, "msnlib.json")).unwrap();
        assert!(parsed.metadata.get("annotation_history").is_none());
    }

    #[test]
    fn gnps_json_inchi_aux_vide_ne_bloque_plus_l_inchikey() {
        let obj = r#"{"spectrum_id":"CCMSLIB00000001547","Compound_Name":"3-Des-Microcystein_LR","peaks_json":"[[289.28,8068.0]]",
            "INCHI_AUX":"N/A","InChIKey_smiles":"IYDKWWDUBYWQGF-NNAZGLEUSA-N","InChIKey_inchi":"","Formula_smiles":"C48H72N10O12","Formula_inchi":""}"#;
        let keys = [("inchi_aux", "INCHIKEY"), ("inchikey_smiles", "INCHIKEY"), ("inchikey_inchi", "INCHIKEY"),
            ("formula_smiles", "FORMULA"), ("formula_inchi", "FORMULA"), ("compound_name", "NAME")];
        let out = mapped(obj, "ALL_GNPS.json", &keys);
        assert_eq!(out["INCHIKEY"], "IYDKWWDUBYWQGF-NNAZGLEUSA-N");
        assert_eq!(out["FORMULA"], "C48H72N10O12");
        // INCHI_AUX rempli : il garde la priorité, comme avant
        let obj2 = obj.replace(r#""INCHI_AUX":"N/A""#, r#""INCHI_AUX":"CGZKSPLDUIRCIO-RPCRKUJJSA-N""#);
        assert_eq!(mapped(&obj2, "ALL_GNPS.json", &keys)["INCHIKEY"], "CGZKSPLDUIRCIO-RPCRKUJJSA-N");
    }
}
