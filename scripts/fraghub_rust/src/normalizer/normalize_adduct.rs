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
    let filename = crate::globals_vars::input_filename(&metadata_dict);
    let inst_type = metadata_dict.get("INSTRUMENTTYPE").map(|s| s.as_str()).unwrap_or("");
    if crate::globals_vars::is_gc_spectrum(filename, inst_type) {
        return metadata_dict; // On ignore les GC
    }

    // `.cloned()` permet d'extraire une copie de la valeur pour la modifier librement.
    if let Some(adduct) = metadata_dict.get("PRECURSORTYPE").cloned() {
        let cleaned = clean_adduct(&adduct);
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

/// Adduit sans crochets, parenthèses ni charge finale (ex. "[M+H]+" -> "M+H", "(M+Na)" -> "M+Na").
/// Un signe collé au M ("M+", "2M-") est gardé : c'est l'ion moléculaire.
fn clean_adduct(adduct: &str) -> String {
    // `.replace_all` renvoie un `Cow` (Copy on Write), `.into_owned()` en fait une `String` mutable.
    let mut cleaned = crate::globals_vars::SUB_ADDUCT_PATTERN.replace_all(adduct, "").into_owned();
    // Remplacement de sub_signe_end_adduct_pattern: r"(?<!M)(\-|\+)$"
    if (cleaned.ends_with('+') || cleaned.ends_with('-')) && !cleaned.ends_with("M+") && !cleaned.ends_with("M-") {
        cleaned.pop(); // Retire très efficacement le dernier caractère de la chaîne.
    }
    cleaned
}

/// Vrai si le texte est un adduit des tables (adduct_to_convert.csv), quelle que soit son écriture.
fn is_known_adduct(text: &str, context: &super::NormalizerContext) -> bool {
    let key = standardize_adduct_key(&clean_adduct(text));
    context.adduct_pos.contains_key(&key) || context.adduct_neg.contains_key(&key)
}

/// Adduit écrit à la fin du nom quand le champ adduit est vide.
///
/// Cas du MGF de GNPS : pas de champ adduit, le nom se termine par l'adduit ("Hectochlorin M+Na",
/// "Cyclomarin A M+H-H2O"). Vaut pour tout fichier qui fait pareil. Le dernier mot du nom est pris
/// comme adduit seulement s'il en a la forme (un M suivi d'un signe) ET s'il est dans les tables
/// d'adduits : un nom ordinaire n'est jamais coupé. L'adduit est alors retiré du nom ("Hectochlorin"),
/// comme dans les exports JSON et MSP de GNPS. Un adduit déjà renseigné n'est jamais remplacé, et les
/// spectres GC (sans adduit) ne sont pas touchés. L'adduit trouvé passe ensuite par les mêmes contrôles
/// que les autres (normalisation, mode d'ionisation, cohérence avec le m/z du précurseur).
pub fn adduct_from_name_if_missing(mut metadata_dict: HashMap<String, String>, context: &super::NormalizerContext) -> HashMap<String, String> {
    if !metadata_dict.get("PRECURSORTYPE").map_or(true, |a| a.trim().is_empty()) {
        return metadata_dict;
    }
    let filename = crate::globals_vars::input_filename(&metadata_dict);
    let inst_type = metadata_dict.get("INSTRUMENTTYPE").map(|s| s.as_str()).unwrap_or("");
    if crate::globals_vars::is_gc_spectrum(filename, inst_type) {
        return metadata_dict;
    }
    let name = match metadata_dict.get("NAME") {
        Some(n) => n.trim_end().to_string(),
        None => return metadata_dict,
    };
    let Some(cut) = name.rfind(char::is_whitespace) else { return metadata_dict };
    let (head, last) = (name[..cut].trim_end(), &name[cut..].trim());
    if head.is_empty() || !crate::globals_vars::IS_ADDUCT_PATTERN.is_match(last) || !is_known_adduct(last, context) {
        return metadata_dict;
    }
    metadata_dict.insert("PRECURSORTYPE".to_string(), last.to_string());
    metadata_dict.insert("NAME".to_string(), head.to_string());
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

    fn context_reel_reduit() -> super::super::NormalizerContext {
        // Quelques adduits des tables, sous leur clé standardisée (comme loading_db.rs)
        let mut pos = HashMap::new();
        for (k, v) in [("M+H", "[M+H]+"), ("M+Na", "[M+Na]+"), ("M+H-H2O", "[M+H-H2O]+"), ("M+", "[M]+"), ("2M+Na", "[2M+Na]+")] {
            pos.insert(standardize_adduct_key(k), v.to_string());
        }
        let mut neg = HashMap::new();
        neg.insert(standardize_adduct_key("M-H"), "[M-H]-".to_string());
        super::super::NormalizerContext {
            adduct_pos: pos, adduct_neg: neg,
            adduct_massdiff_pos: HashMap::new(), adduct_massdiff_neg: HashMap::new(),
            instrument_tree: serde_json::Value::Null,
        }
    }

    fn from_name(name: &str, adduct: &str, inst: &str) -> (String, String) {
        let mut m = HashMap::new();
        m.insert("NAME".to_string(), name.to_string());
        m.insert("PRECURSORTYPE".to_string(), adduct.to_string());
        m.insert("FILENAME".to_string(), "ALL_GNPS.mgf".to_string());
        m.insert("INSTRUMENTTYPE".to_string(), inst.to_string());
        let m = adduct_from_name_if_missing(m, &context_reel_reduit());
        (m["NAME"].clone(), m["PRECURSORTYPE"].clone())
    }

    #[test]
    fn adduit_lu_a_la_fin_du_nom() {
        // Noms réels du MGF GNPS
        assert_eq!(from_name("3-Des-Microcystein_LR M+H", "", "LC-ESI-qTof"), ("3-Des-Microcystein_LR".into(), "M+H".into()));
        assert_eq!(from_name("Cyclomarin A M+H-H2O", "", ""), ("Cyclomarin A".into(), "M+H-H2O".into()));
        assert_eq!(from_name("Hectochlorin M+Na", "", ""), ("Hectochlorin".into(), "M+Na".into()));
        assert_eq!(from_name("unpublished peptide from CNL 643 M+H", "", ""), ("unpublished peptide from CNL 643".into(), "M+H".into()));
        // Autres écritures
        assert_eq!(from_name("Glucose [M-H]-", "", ""), ("Glucose".into(), "[M-H]-".into()));
        assert_eq!(from_name("Glucose (M+Na)", "", ""), ("Glucose".into(), "(M+Na)".into()));
        assert_eq!(from_name("Truc 2M+Na  ", "", ""), ("Truc".into(), "2M+Na".into()));
        // Normalisé ensuite comme tout adduit
        let mut m = HashMap::new();
        m.insert("PRECURSORTYPE".to_string(), "M+H-H2O".to_string());
        m.insert("FILENAME".to_string(), "a.mgf".to_string());
        m.insert("INSTRUMENTTYPE".to_string(), "".to_string());
        assert_eq!(normalize_adduct(m, &context_reel_reduit())["PRECURSORTYPE"], "[M+H-H2O]+");
    }

    #[test]
    fn nom_sans_adduit_inchange() {
        // Dernier mot qui n'est pas un adduit des tables, ou nom d'un seul mot
        assert_eq!(from_name("Vitamin M", "", ""), ("Vitamin M".into(), "".into()));
        assert_eq!(from_name("Compound X-M+Zz", "", ""), ("Compound X-M+Zz".into(), "".into()));
        assert_eq!(from_name("Tryptophan", "", ""), ("Tryptophan".into(), "".into()));
        assert_eq!(from_name("M+H", "", ""), ("M+H".into(), "".into()));
        assert_eq!(from_name("Acide 2-amino-3-methyl", "", ""), ("Acide 2-amino-3-methyl".into(), "".into()));
        // Adduit déjà renseigné : ni l'adduit ni le nom ne changent
        assert_eq!(from_name("Hectochlorin M+Na", "[M+H]+", ""), ("Hectochlorin M+Na".into(), "[M+H]+".into()));
        // Spectre GC : pas touché
        assert_eq!(from_name("Hectochlorin M+Na", "", "GC-EI-TOF"), ("Hectochlorin M+Na".into(), "".into()));
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
