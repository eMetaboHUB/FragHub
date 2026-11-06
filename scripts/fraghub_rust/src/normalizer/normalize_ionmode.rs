// src/normalizer/normalize_ionmode.rs
use std::collections::HashMap;

/// Mode d'ionisation lu dans le texte : "positive", "negative", ou `None` si le texte ne dit rien.
///
/// Règle d'origine : le texte commence par p / pos / + (positif) ou n / neg / - (négatif).
/// En plus : les espaces autour sont retirés (GNPS : " Positive"), et si le début ne dit rien, on cherche
/// le mot positive / negative ailleurs dans le texte (commentaires MoNA : "ionmode=ESI Positive").
/// (Les valeurs vides comme "N/A" ou "unknown" ont déjà été vidées par normalize_empties.)
pub fn ion_mode_from_text(raw: &str) -> Option<&'static str> {
    let t = raw.trim();
    if t.is_empty() {
        return None;
    }
    if crate::globals_vars::IONMODE_POS_PATTERN.is_match(t) {
        return Some("positive");
    }
    if crate::globals_vars::IONMODE_NEG_PATTERN.is_match(t) {
        return Some("negative");
    }
    let pos = crate::globals_vars::IONMODE_POS_WORD_PATTERN.is_match(t);
    let neg = crate::globals_vars::IONMODE_NEG_WORD_PATTERN.is_match(t);
    match (pos, neg) {
        (true, false) => Some("positive"),
        (false, true) => Some("negative"),
        _ => None, // aucun des deux, ou les deux à la fois : illisible
    }
}

/// Mode d'ionisation déduit de la charge écrite après le dernier "]" de l'adduit :
/// "[M+H]+" ou "[M+2H]2+" -> positive, "[M-H]-" -> negative. `None` sans crochet, sans signe, ou signes mêlés.
pub fn ion_mode_from_adduct(adduct: &str) -> Option<&'static str> {
    let t = adduct.trim();
    let close = t.rfind(']')?;
    let signs: String = t[close + 1..].chars().filter(|c| *c == '+' || *c == '-').collect();
    if signs.is_empty() {
        None
    } else if signs.chars().all(|c| c == '+') {
        Some("positive")
    } else if signs.chars().all(|c| c == '-') {
        Some("negative")
    } else {
        None
    }
}

/// Standardise le champ IONMODE en "positive" ou "negative".
/// Un texte illisible est laissé tel quel (le mode pourra être déduit de l'adduit, voir plus bas).
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec IONMODE standardisé.
pub fn normalize_ion_mode(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    if let Some(mode) = metadata_dict.get("IONMODE").and_then(|m| ion_mode_from_text(m)) {
        metadata_dict.insert("IONMODE".to_string(), mode.to_string());
    }
    metadata_dict
}

/// Mode d'ionisation absent ou illisible : on le déduit du signe de l'adduit (déjà normalisé).
/// Avant, ces spectres arrivaient au splitter sans mode et y disparaissaient sans laisser de trace.
/// Spectre GC sans adduit exploitable : positif (l'ionisation électronique donne des ions positifs).
/// Sinon IONMODE reste tel quel et check_for_bad_adduct supprime le spectre (avec un motif).
pub fn ion_mode_from_adduct_if_missing(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    let known = matches!(metadata_dict.get("IONMODE").map(|s| s.as_str()), Some("positive") | Some("negative"));
    if !known {
        let from_adduct = metadata_dict.get("PRECURSORTYPE").and_then(|a| ion_mode_from_adduct(a));
        let mode = from_adduct.or_else(|| {
            let filename = metadata_dict.get("FILENAME").map(|s| s.as_str()).unwrap_or("");
            let instrument_type = metadata_dict.get("INSTRUMENTTYPE").map(|s| s.as_str()).unwrap_or("");
            crate::globals_vars::is_gc_spectrum(filename, instrument_type).then_some("positive")
        });
        if let Some(mode) = mode {
            metadata_dict.insert("IONMODE".to_string(), mode.to_string());
        }
    }
    metadata_dict
}

#[cfg(test)]
mod tests {
    use super::*;

    fn with(fields: &[(&str, &str)]) -> HashMap<String, String> {
        fields.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect()
    }

    #[test]
    fn mode_lu_dans_le_texte() {
        // Règle d'origine, inchangée
        assert_eq!(ion_mode_from_text("Positive"), Some("positive"));
        assert_eq!(ion_mode_from_text("P"), Some("positive"));
        assert_eq!(ion_mode_from_text("+"), Some("positive"));
        assert_eq!(ion_mode_from_text("Positive-40eV"), Some("positive"));
        assert_eq!(ion_mode_from_text("negative"), Some("negative"));
        assert_eq!(ion_mode_from_text("N"), Some("negative"));
        assert_eq!(ion_mode_from_text("-"), Some("negative"));
        // Nouveaux cas : espace devant (GNPS), mot ailleurs dans le texte (MoNA)
        assert_eq!(ion_mode_from_text(" Positive"), Some("positive"));
        assert_eq!(ion_mode_from_text(" Negative"), Some("negative"));
        assert_eq!(ion_mode_from_text("ESI Positive"), Some("positive"));
        assert_eq!(ion_mode_from_text("ESI Negative"), Some("negative"));
        // Illisible
        assert_eq!(ion_mode_from_text(""), None);
        assert_eq!(ion_mode_from_text("E"), None);
        assert_eq!(ion_mode_from_text("ESI positive/negative"), None);
    }

    #[test]
    fn mode_deduit_de_l_adduit() {
        assert_eq!(ion_mode_from_adduct("[M+H]+"), Some("positive"));
        assert_eq!(ion_mode_from_adduct("[M+2H]2+"), Some("positive"));
        assert_eq!(ion_mode_from_adduct("[M]+*"), Some("positive"));
        assert_eq!(ion_mode_from_adduct("[M-H]-"), Some("negative"));
        assert_eq!(ion_mode_from_adduct("[M-2H]2-"), Some("negative"));
        assert_eq!(ion_mode_from_adduct("[(M+CH3COOH)-H]-"), Some("negative"));
        assert_eq!(ion_mode_from_adduct("M+H"), None); // sans crochets : on ne devine pas
        assert_eq!(ion_mode_from_adduct("[M+H]"), None);
        assert_eq!(ion_mode_from_adduct(""), None);
    }

    #[test]
    fn deduction_seulement_si_le_mode_manque() {
        let m = ion_mode_from_adduct_if_missing(with(&[("IONMODE", ""), ("PRECURSORTYPE", "[M+H]+")]));
        assert_eq!(m["IONMODE"], "positive");
        let m = ion_mode_from_adduct_if_missing(with(&[("PRECURSORTYPE", "[M-H]-")]));
        assert_eq!(m["IONMODE"], "negative");
        let m = ion_mode_from_adduct_if_missing(with(&[("IONMODE", "E"), ("PRECURSORTYPE", "[M+H]+")]));
        assert_eq!(m["IONMODE"], "positive");
        // Mode déjà connu : jamais remplacé (l'incohérence éventuelle est traitée par check_for_bad_adduct)
        let m = ion_mode_from_adduct_if_missing(with(&[("IONMODE", "negative"), ("PRECURSORTYPE", "[M+H]+")]));
        assert_eq!(m["IONMODE"], "negative");
        // Pas d'adduit exploitable : inchangé
        let m = ion_mode_from_adduct_if_missing(with(&[("IONMODE", ""), ("PRECURSORTYPE", "")]));
        assert_eq!(m["IONMODE"], "");
        // Spectre GC sans mode ni adduit : positif d'office
        let m = ion_mode_from_adduct_if_missing(with(&[("IONMODE", ""), ("PRECURSORTYPE", ""), ("INSTRUMENTTYPE", "GC-EI-TOF")]));
        assert_eq!(m["IONMODE"], "positive");
        // Spectre GC dont l'adduit dit négatif : négatif
        let m = ion_mode_from_adduct_if_missing(with(&[("IONMODE", ""), ("PRECURSORTYPE", "[M-H]-"), ("INSTRUMENTTYPE", "GC-APCI-QTOF")]));
        assert_eq!(m["IONMODE"], "negative");
    }
}
