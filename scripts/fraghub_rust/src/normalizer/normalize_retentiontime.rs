use std::collections::HashMap;

/// Convertit le temps de rétention en minutes (float) peu importe l'unité d'origine.
///
/// Pour un développeur Python : Observez le mot-clé `match`. C'est le cousin très puissant
/// du `switch`/`case` ou du récent `match` en Python 3.10. Le compilateur Rust nous obligera
/// TOUJOURS à traiter tous les cas possibles (exhaustivité), d'où le `_ =>` (le cas "else" par défaut).
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec RT normalisé en minutes.
pub fn normalize_retention_time(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    let rt = metadata_dict.get("RT").cloned().unwrap_or_default();

    if let Some(caps) = crate::globals_vars::RETENTION_TIME_PATTERN.captures(&rt) {
        if let Some(time_match) = caps.get(1) {
            // `.parse::<f64>()` remplace le `float(valeur)` de Python. 
            // Si la conversion réussit (`Ok`), on rentre dans le bloc.
            if let Ok(time_val) = time_match.as_str().replace(',', ".").parse::<f64>() { // virgule décimale acceptée
                // On récupère l'unité (groupe 2 de la regex) ou "" par défaut
                let unit = caps.get(2).map(|m| m.as_str().to_lowercase()).unwrap_or_default();
                
                // Le `match` (pattern matching) pour la conversion d'unités
                let final_rt = match unit.as_str() {
                    "m" | "min" | "minute" | "minutes" | "" => time_val, // Déjà en minutes
                    "s" | "sec" | "second" | "seconds" => time_val / 60.0, // Secondes -> Minutes
                    "ms" | "millisecond" | "milliseconds" => time_val / 60000.0, // Millisecondes -> Minutes
                    _ => time_val, // Par sécurité (le compilateur l'exige)
                };
                
                metadata_dict.insert("RT".to_string(), final_rt.to_string());
            }
        }
    }
    metadata_dict
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rt(v: &str) -> String {
        let mut m = HashMap::new();
        m.insert("RT".to_string(), v.to_string());
        normalize_retention_time(m)["RT"].clone()
    }

    #[test]
    fn unites_converties_en_minutes() {
        assert_eq!(rt("70.99 s"), (70.99_f64 / 60.0).to_string()); // RTINSECONDS marqué à la lecture
        assert_eq!(rt("90 sec"), "1.5");
        assert_eq!(rt("12 ms"), (12.0_f64 / 60000.0).to_string()); // était lu comme 12 minutes
        assert_eq!(rt("120000 milliseconds"), "2");
        assert_eq!(rt("5.2 min"), "5.2");
        assert_eq!(rt("5 minutes"), "5");
        assert_eq!(rt("3 m"), "3");
        assert_eq!(rt("12,5 s"), (12.5_f64 / 60.0).to_string()); // virgule décimale
        assert_eq!(rt("5,2 min"), "5.2");
    }

    #[test]
    fn nombre_seul_reste_en_minutes() {
        assert_eq!(rt("5.2"), "5.2");
        assert_eq!(rt("61.971333"), "61.971333");
        assert_eq!(rt(""), "");
    }
}
