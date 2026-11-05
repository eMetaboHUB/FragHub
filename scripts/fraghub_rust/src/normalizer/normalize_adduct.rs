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
            // Adduit absent des tables : on le garde sous la forme standard "[...]charge".
            // Avant, le texte nettoyé était gardé tel quel : "[M+HCOOH-H]-" devenait "M+HCOOH-H".
            let ion_mode = metadata_dict.get("IONMODE").cloned().unwrap_or_default();
            metadata_dict.insert("PRECURSORTYPE".to_string(), unknown_adduct_form(&adduct, &cleaned, &ion_mode));
        }
    }

    metadata_dict
}

/// Forme écrite d'un adduit absent des tables : "[" + adduit nettoyé + "]" + charge.
///
/// La charge est celle du texte d'origine (après le dernier "]", ou signe final sans crochets :
/// "M+CH3COO-"), sinon celle du mode d'ionisation ("+" / "-"), sinon elle est omise.
/// Un signe collé au M ("M+", "2M-") est la charge de l'ion moléculaire : "[M]+", "[2M]-".
/// Un adduit vide reste vide.
fn unknown_adduct_form(original: &str, cleaned: &str, ion_mode: &str) -> String {
    if cleaned.is_empty() {
        return String::new();
    }
    if cleaned.ends_with("M+") || cleaned.ends_with("M-") {
        let (body, sign) = cleaned.split_at(cleaned.len() - 1);
        return format!("[{}]{}", body, sign);
    }
    let written = match original.rfind(']') {
        Some(close) => original[close + 1..].trim().to_string(),
        None => {
            // Sans crochets : le nettoyage a pu retirer un signe final, c'est la charge écrite
            let stripped = crate::globals_vars::SUB_ADDUCT_PATTERN.replace_all(original, "");
            stripped.get(cleaned.len()..).unwrap_or("").trim().to_string()
        }
    };
    let charge = if crate::globals_vars::ADDUCT_CHARGE_PATTERN.is_match(&written) {
        written
    } else {
        match ion_mode {
            "positive" => "+".to_string(),
            "negative" => "-".to_string(),
            _ => String::new(),
        }
    };
    format!("[{}]{}", cleaned, charge)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn context() -> super::super::NormalizerContext {
        let mut pos = HashMap::new();
        pos.insert("M+H".to_string(), "[M+H]+".to_string());
        super::super::NormalizerContext {
            adduct_pos: pos, adduct_neg: HashMap::new(),
            adduct_massdiff_pos: HashMap::new(), adduct_massdiff_neg: HashMap::new(),
            instrument_tree: serde_json::Value::Null,
        }
    }

    fn adduct(raw: &str, ion_mode: &str) -> String {
        let mut m = HashMap::new();
        m.insert("PRECURSORTYPE".to_string(), raw.to_string());
        m.insert("IONMODE".to_string(), ion_mode.to_string());
        m.insert("FILENAME".to_string(), "a.msp".to_string());
        m.insert("INSTRUMENTTYPE".to_string(), "LC-ESI-QTOF".to_string());
        normalize_adduct(m, &context())["PRECURSORTYPE"].clone()
    }

    #[test]
    fn adduit_connu_inchange() {
        assert_eq!(adduct("[M+H]+", "positive"), "[M+H]+");
        assert_eq!(adduct("M+H", ""), "[M+H]+");
    }

    #[test]
    fn adduit_inconnu_garde_crochets_et_charge() {
        // Cas réels de l'échantillon (MoNA, MassBank, GNPS)
        assert_eq!(adduct("[M+HCOOH-H]-", "negative"), "[M+HCOOH-H]-");
        assert_eq!(adduct("[M+Hac-H]-", "negative"), "[M+Hac-H]-");
        assert_eq!(adduct("[4M+Ca]2+", "positive"), "[4M+Ca]2+");
        assert_eq!(adduct("[M-C3H6NO2]-", "negative"), "[M-C3H6NO2]-");
        // Sans crochets ni charge : charge prise du mode d'ionisation
        assert_eq!(adduct("M-H1", "negative"), "[M-H1]-");
        assert_eq!(adduct("M+NH3", "positive"), "[M+NH3]+");
        assert_eq!(adduct("M+NH3", ""), "[M+NH3]");
        // Sans crochets, signe final écrit
        assert_eq!(adduct("M+CH3COO-", "negative"), "[M+CH3COO]-");
        // Ion moléculaire
        assert_eq!(adduct("M", "positive"), "[M]+");
        assert_eq!(adduct("M+", ""), "[M]+");
        // Vide : reste vide
        assert_eq!(adduct("", "positive"), "");
    }

    #[test]
    fn spectre_gc_non_touche() {
        let mut m = HashMap::new();
        m.insert("PRECURSORTYPE".to_string(), "M-H1".to_string());
        m.insert("FILENAME".to_string(), "a.msp".to_string());
        m.insert("INSTRUMENTTYPE".to_string(), "GC-EI-TOF".to_string());
        assert_eq!(normalize_adduct(m, &context())["PRECURSORTYPE"], "M-H1");
    }
}
