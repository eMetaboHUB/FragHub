// src/normalizer/normalize_adduct.rs
use std::collections::HashMap;


/// Standardise l'ordre des éléments d'un adduit après la partie M
/// Par exemple: 2M-2H+Ca-H devient 2M+Ca-2H-H (trié de manière déterministe)
pub fn standardize_adduct_key(adduct: &str) -> String {
    let first_sign_idx = adduct.find(|c| c == '+' || c == '-');
    match first_sign_idx {
        Some(idx) => {
            let (base, rest) = adduct.split_at(idx);
            let mut parts: Vec<String> = Vec::new();
            let mut current_part = String::new();
            for c in rest.chars() {
                if c == '+' || c == '-' {
                    if !current_part.is_empty() {
                        parts.push(current_part.clone());
                    }
                    current_part.clear();
                }
                current_part.push(c);
            }
            if !current_part.is_empty() {
                parts.push(current_part);
            }
            parts.sort_unstable();
            format!("{}{}", base, parts.join(""))
        }
        None => adduct.to_string(),
    }
}

/// Normalise la forme de l'adduct (ex: "[M+H]+" au lieu de "M+H").
///
/// Pour un développeur Python : Remarquez comment le `NormalizerContext` est passé par référence (`&super::NormalizerContext`).
/// Cela permet à cette fonction de consulter les grands dictionnaires d'adducts en RAM sans les copier.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Dictionnaire contenant les métadonnées.
/// * `context` (&super::NormalizerContext) : Contexte global pour vérifier la validité de l'adduct.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec l'adduct standardisé.
pub fn normalize_adduct(mut metadata_dict: HashMap<String, String>, context: &super::NormalizerContext) -> HashMap<String, String> {
    let filename = metadata_dict.get("FILENAME").map(|s| s.as_str()).unwrap_or("");
    let inst_type = metadata_dict.get("INSTRUMENTTYPE").map(|s| s.as_str()).unwrap_or("");
    if crate::globals_vars::is_gc_spectrum(filename, inst_type) {
        return metadata_dict; // On ignore les GC
    }

    // `.cloned()` permet d'extraire une copie de la valeur pour la modifier librement.
    if let Some(adduct) = metadata_dict.get("PRECURSORTYPE").cloned() {
        // `.replace_all` renvoie un `Cow` (Copy on Write), `.into_owned()` en fait une `String` mutable.
        let mut cleaned = crate::globals_vars::SUB_ADDUCT_PATTERN.replace_all(&adduct, "").into_owned();

        // Remplacement de sub_signe_end_adduct_pattern: r"(?<!M)(\-|\+)$"
        if (cleaned.ends_with('+') || cleaned.ends_with('-')) && !cleaned.ends_with("M+") && !cleaned.ends_with("M-") {
            cleaned.pop(); // Retire très efficacement le dernier caractère de la chaîne.
        }

        let standard_key = standardize_adduct_key(&cleaned);

        // On consulte le contexte avec la clé standardisée
        if let Some(canonical) = context.adduct_pos.get(&standard_key) {
            metadata_dict.insert("PRECURSORTYPE".to_string(), canonical.clone());
        } else if let Some(canonical) = context.adduct_neg.get(&standard_key) {
            metadata_dict.insert("PRECURSORTYPE".to_string(), canonical.clone());
        } else {
            metadata_dict.insert("PRECURSORTYPE".to_string(), cleaned);
        }
    }

    metadata_dict
}
