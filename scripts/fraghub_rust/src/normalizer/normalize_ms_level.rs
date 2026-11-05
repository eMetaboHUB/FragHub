use std::collections::HashMap;

/// Standardise le niveau MS (MS1, MS2, MS3, etc.).
///
/// Seul le chiffre est gardé : "MS2", "ms2" et "2" donnent tous "2", "MS2/MS3" donne "2-3".
///
/// Pour un développeur Python : `.captures_iter()` sur une regex en Rust renvoie un itérateur de captures.
/// On prend le groupe 1 (le chiffre) de chacune, puis `.collect()` rassemble le tout dans un vecteur
/// (`Vec<&str>`). C'est l'équivalent de `[m.group(1) for m in pattern.finditer(text)]`.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec le niveau MS standardisé.
pub fn normalize_ms_level(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    let ms_level = metadata_dict.get("MSLEVEL").cloned().unwrap_or_default();

    if !ms_level.is_empty() {
        // Extraction de tous les niveaux MS mentionnés (ex: "MS2/MS3" -> ["2", "3"]).
        // On garde le groupe 1 (le chiffre) et non tout le texte trouvé, sinon "MS2" restait "MS2"
        // alors que "2" restait "2" : les deux écritures se mélangeaient dans les sorties.
        let matched_levels: Vec<&str> = crate::globals_vars::MS_LEVEL_PATTERN
            .captures_iter(&ms_level)
            .filter_map(|c| c.get(1).map(|m| m.as_str()))
            .collect();

        // Si on a trouvé exactement 1 niveau, on l'utilise.
        if matched_levels.len() == 1 {
            metadata_dict.insert("MSLEVEL".to_string(), matched_levels[0].to_string());
        } 
        // Si on en a trouvé plusieurs, on fait un format combiné (ex: "2-3")
        else if matched_levels.len() >= 2 {
            // `format!` est l'équivalent de f"{matched_levels[0]}-{matched_levels[1]}"
            metadata_dict.insert("MSLEVEL".to_string(), format!("{}-{}", matched_levels[0], matched_levels[1]));
        } 
        // Par défaut, si rien ne match (ou MS inconnu), on suppose MS2.
        else {
            metadata_dict.insert("MSLEVEL".to_string(), "2".to_string());
        }
    } else {
        // Par défaut, on suppose MS2 si le champ est vide.
        metadata_dict.insert("MSLEVEL".to_string(), "2".to_string());
    }

    metadata_dict
}

#[cfg(test)]
mod tests {
    use super::*;

    fn level(raw: &str) -> String {
        let mut m = HashMap::new();
        m.insert("MSLEVEL".to_string(), raw.to_string());
        normalize_ms_level(m)["MSLEVEL"].clone()
    }

    #[test]
    fn niveau_ms_reduit_au_chiffre() {
        assert_eq!(level("MS2"), "2");
        assert_eq!(level("ms2"), "2");
        assert_eq!(level("2"), "2");
        assert_eq!(level("MS1"), "1");
        assert_eq!(level("MS3"), "3");
        assert_eq!(level("MS2/MS3"), "2-3");
    }

    #[test]
    fn niveau_ms_absent_ou_illisible_vaut_2() {
        assert_eq!(level(""), "2");
        assert_eq!(level("MSn"), "2");
        assert_eq!(normalize_ms_level(HashMap::new())["MSLEVEL"], "2");
    }
}
