// src/normalizer/check_for_bad_adduct.rs
use std::collections::HashMap;

/// Motif de suppression d'un spectre dont le mode d'ionisation reste inconnu (ni écrit, ni déductible de l'adduit).
pub const IONMODE_MISSING_REASON: &str =
    "spectrum deleted because its ionization mode is missing or invalid and cannot be deduced from the adduct";

/// Complète l'adduit quand il manque ou est écrit en abrégé. Passe AVANT le recalcul du m/z manquant,
/// qui a besoin de l'adduit :
/// - spectre in silico sans adduit : [M+H]+ en mode positif, [M-H]- en mode négatif ;
/// - forme courte "M" : [M]+ ou [M]- selon le mode.
/// Avant, ces compléments étaient faits dans check_for_bad_adduct, donc après le recalcul : le m/z
/// manquant de ces spectres ne pouvait jamais être recalculé.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Dictionnaire contenant les métadonnées.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire, avec l'adduit complété si besoin.
pub fn complete_adduct(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    let adduct = metadata_dict.get("PRECURSORTYPE").cloned().unwrap_or_default();
    let predicted = metadata_dict.get("PREDICTED").cloned().unwrap_or_default();
    let positive = match metadata_dict.get("IONMODE").map(|s| s.as_str()) {
        Some("positive") => true,
        Some("negative") => false,
        _ => return metadata_dict, // mode inconnu : rien à compléter
    };

    if predicted == "true" && adduct.is_empty() {
        let default = if positive { "[M+H]+" } else { "[M-H]-" };
        metadata_dict.insert("PRECURSORTYPE".to_string(), default.to_string());
    } else if adduct == "M" {
        let molecular_ion = if positive { "[M]+" } else { "[M]-" };
        metadata_dict.insert("PRECURSORTYPE".to_string(), molecular_ion.to_string());
    }
    metadata_dict
}

/// Vérifie que l'adduct est cohérent avec le mode d'ionisation, et que ce mode est connu.
///
/// Pour un développeur Python : Retourne `Option<HashMap>` pour propager une suppression silencieusement.
/// Le `deletion_reason` est passé par `&mut` pour permettre à la fonction de modifier la raison
/// visible par le rapporteur final en cas de suppression (renvoi de `None`).
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Dictionnaire contenant les métadonnées.
/// * `deletion_reason` (&mut Option<String>) : Pointeur mutable pour stocker la raison de suppression.
/// * `context` (&NormalizerContext) : Dictionnaires de référence pour les adducts.
///
/// # Returns
/// * `Option<HashMap<String, String>>` : Dictionnaire valide ou `None` si incohérent.
pub fn check_for_bad_adduct(
    metadata_dict: HashMap<String, String>,
    deletion_reason: &mut Option<String>,
    context: &super::NormalizerContext
) -> Option<HashMap<String, String>> {

    let adduct = metadata_dict.get("PRECURSORTYPE").cloned().unwrap_or_default();
    let ion_mode = metadata_dict.get("IONMODE").cloned().unwrap_or_default();
    let instrument_type = metadata_dict.get("INSTRUMENTTYPE").cloned().unwrap_or_default();
    let filename = crate::globals_vars::input_filename(&metadata_dict).to_string();

    // 1. Spectres GC (ionisation électronique) : pas d'adduit, c'est normal
    if crate::globals_vars::is_gc_spectrum(&filename, &instrument_type) && adduct.is_empty() {
        return keep_if_ion_mode_known(metadata_dict, &ion_mode, deletion_reason);
    }

    // 2. Valider le format de l'adduct
    if !crate::globals_vars::IS_ADDUCT_PATTERN.is_match(&adduct) {
        *deletion_reason = Some("spectrum deleted because its adduct field is empty or the value entered is not an adduct".to_string());
        return None; // Suppression du spectre
    }

    // 3. Valider la cohérence Adduct / Ion Mode grâce aux dictionnaires en RAM (`context.adduct_massdiff_neg`)
    if ion_mode == "positive" {
        if context.adduct_massdiff_neg.contains_key(&adduct) {
            *deletion_reason = Some("spectrum deleted because the adduct corresponds to the wrong ionization mode (neg adduct in pos ionmode).".to_string());
            return None;
        }
    } else if ion_mode == "negative" {
        if context.adduct_massdiff_pos.contains_key(&adduct) {
            *deletion_reason = Some("spectrum deleted because the adduct corresponds to the wrong ionization mode (pos adduct in neg ionmode).".to_string());
            return None;
        }
    }

    // 4. Le mode d'ionisation doit être connu (sinon le spectre ne peut être rangé ni en POS ni en NEG :
    //    avant, il disparaissait au splitter sans apparaître dans DELETED_SPECTRUMS)
    keep_if_ion_mode_known(metadata_dict, &ion_mode, deletion_reason)
}

/// Garde le spectre si son mode d'ionisation est "positive" ou "negative", sinon le supprime avec un motif.
fn keep_if_ion_mode_known(
    metadata_dict: HashMap<String, String>,
    ion_mode: &str,
    deletion_reason: &mut Option<String>,
) -> Option<HashMap<String, String>> {
    if ion_mode == "positive" || ion_mode == "negative" {
        Some(metadata_dict)
    } else {
        *deletion_reason = Some(IONMODE_MISSING_REASON.to_string());
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn context() -> super::super::NormalizerContext {
        super::super::NormalizerContext {
            adduct_pos: HashMap::new(), adduct_neg: HashMap::new(),
            adduct_massdiff_pos: HashMap::new(), adduct_massdiff_neg: HashMap::new(),
            instrument_tree: serde_json::Value::Null,
        }
    }

    fn fields(f: &[(&str, &str)]) -> HashMap<String, String> {
        f.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect()
    }

    fn check(f: &[(&str, &str)]) -> bool {
        let mut reason = None;
        check_for_bad_adduct(fields(f), &mut reason, &context()).is_some()
    }

    fn reason(f: &[(&str, &str)]) -> Option<String> {
        let mut reason = None;
        check_for_bad_adduct(fields(f), &mut reason, &context());
        reason
    }

    #[test]
    fn spectres_gc_et_ei_sans_adduit_gardes() {
        // Cas réels supprimés à tort avant (MoNA, MassBank) : ionisation électronique, pas d'adduit
        assert!(check(&[("INSTRUMENTTYPE", "EI-B"), ("IONMODE", "positive"), ("FILENAME", "MoNA.msp")]));
        assert!(check(&[("INSTRUMENTTYPE", "EI-EBEB"), ("IONMODE", "positive"), ("FILENAME", "MassBank.msp")]));
        assert!(check(&[("INSTRUMENTTYPE", "GC-EI-TOF"), ("IONMODE", "positive"), ("FILENAME", "a.msp")]));
        // Fichier marqué _GC (README), type d'instrument inconnu
        assert!(check(&[("INSTRUMENTTYPE", ""), ("IONMODE", "positive"), ("FILENAME", "lib_GC.msp")]));
    }

    #[test]
    fn spectre_lc_sans_adduit_toujours_supprime() {
        assert!(!check(&[("INSTRUMENTTYPE", "LC-ESI-QTOF"), ("IONMODE", "positive"), ("FILENAME", "a.msp")]));
        assert!(!check(&[("INSTRUMENTTYPE", ""), ("IONMODE", "positive"), ("FILENAME", "a.msp")]));
    }

    #[test]
    fn mode_inconnu_supprime_avec_un_motif() {
        // Avant : gardé ici, puis perdu sans trace au splitter
        let r = reason(&[("PRECURSORTYPE", "[M+H]"), ("IONMODE", ""), ("FILENAME", "a.msp")]);
        assert_eq!(r.as_deref(), Some(IONMODE_MISSING_REASON));
        let r = reason(&[("INSTRUMENTTYPE", "GC-EI-TOF"), ("IONMODE", ""), ("FILENAME", "a.msp")]);
        assert_eq!(r.as_deref(), Some(IONMODE_MISSING_REASON));
        // Adduit invalide : le motif « adduit » reste prioritaire, comme avant
        let r = reason(&[("PRECURSORTYPE", ""), ("IONMODE", ""), ("FILENAME", "a.msp")]);
        assert!(r.unwrap().contains("adduct field is empty"));
        // Mode connu : gardé
        assert!(check(&[("PRECURSORTYPE", "[M+H]+"), ("IONMODE", "positive"), ("FILENAME", "a.msp")]));
    }

    #[test]
    fn adduit_complete_avant_le_recalcul_du_mz() {
        let m = complete_adduct(fields(&[("PREDICTED", "true"), ("PRECURSORTYPE", ""), ("IONMODE", "positive")]));
        assert_eq!(m["PRECURSORTYPE"], "[M+H]+");
        let m = complete_adduct(fields(&[("PREDICTED", "true"), ("PRECURSORTYPE", ""), ("IONMODE", "negative")]));
        assert_eq!(m["PRECURSORTYPE"], "[M-H]-");
        let m = complete_adduct(fields(&[("PREDICTED", "false"), ("PRECURSORTYPE", "M"), ("IONMODE", "negative")]));
        assert_eq!(m["PRECURSORTYPE"], "[M]-");
        // Expérimental sans adduit, ou mode inconnu : inchangé
        let m = complete_adduct(fields(&[("PREDICTED", "false"), ("PRECURSORTYPE", ""), ("IONMODE", "positive")]));
        assert_eq!(m["PRECURSORTYPE"], "");
        let m = complete_adduct(fields(&[("PREDICTED", "true"), ("PRECURSORTYPE", ""), ("IONMODE", "")]));
        assert_eq!(m["PRECURSORTYPE"], "");
        // Adduit déjà présent : inchangé
        let m = complete_adduct(fields(&[("PREDICTED", "true"), ("PRECURSORTYPE", "[M+Na]+"), ("IONMODE", "positive")]));
        assert_eq!(m["PRECURSORTYPE"], "[M+Na]+");
    }
}
