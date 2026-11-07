// src/normalizer/normalize_name.rs
use std::collections::HashMap;

/// Retire de la fin du nom les morceaux qui ne sont pas le nom de la molécule.
///
/// MS-DIAL (MSMS-Public) écrit parfois l'énergie de collision et l'InChIKey dans le nom, séparés par « ; » :
/// "Piperine; CE0; MXXWOMGUGJBKIW-YPCIICBESA-N" -> "Piperine", "Betaine; CE30; CE10; KWIU…-N" -> "Betaine".
/// Seuls les derniers morceaux qui sont exactement une InChIKey ou "CE" suivi d'un nombre sont retirés :
/// deux noms séparés par « ; » ("Isosaxalin; Heraclenol") restent tels quels, et le premier morceau est
/// toujours gardé.
pub fn clean_name_suffixes(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    if let Some(name) = metadata_dict.get("NAME") {
        let cleaned = strip_suffixes(name);
        if cleaned.len() != name.len() {
            metadata_dict.insert("NAME".to_string(), cleaned);
        }
    }
    metadata_dict
}

fn strip_suffixes(name: &str) -> String {
    let mut end = name.trim_end();
    while let Some(cut) = end.rfind(';') {
        let last = end[cut + 1..].trim();
        let head = end[..cut].trim_end();
        if head.is_empty() || !(is_inchikey(last) || is_collision_energy(last)) {
            break;
        }
        end = head;
    }
    if end.len() == name.trim_end().len() { name.to_string() } else { end.to_string() }
}

fn is_inchikey(s: &str) -> bool {
    crate::globals_vars::INCHIKEY_FULL_PATTERN.is_match(s)
}

fn is_collision_energy(s: &str) -> bool {
    let rest = match s.strip_prefix("CE") { Some(r) => r, None => return false };
    !rest.is_empty() && rest.replace(',', ".").parse::<f64>().map_or(false, |x| x.is_finite())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn name(n: &str) -> String {
        let mut m = HashMap::new();
        m.insert("NAME".to_string(), n.to_string());
        clean_name_suffixes(m)["NAME"].clone()
    }

    #[test]
    fn energie_et_inchikey_retirees_du_nom() {
        // Noms réels de MSMS-Public (MS-DIAL)
        assert_eq!(name("Piperine; CE0; MXXWOMGUGJBKIW-YPCIICBESA-N"), "Piperine");
        assert_eq!(name("Betaine; CE30; CE10; KWIUHFFTVRNATP-UHFFFAOYSA-N"), "Betaine");
        assert_eq!(name("Truc; MXXWOMGUGJBKIW-YPCIICBESA-N"), "Truc");
        assert_eq!(name("Truc; CE12.5"), "Truc");
    }

    #[test]
    fn autres_noms_inchanges() {
        assert_eq!(name("Isosaxalin; Heraclenol"), "Isosaxalin; Heraclenol");
        assert_eq!(name("Chrysin"), "Chrysin");
        assert_eq!(name("PC 34:1; CE"), "PC 34:1; CE");
        assert_eq!(name("Acide; CEtone"), "Acide; CEtone");
        // Le nom n'est jamais vidé
        assert_eq!(name("MXXWOMGUGJBKIW-YPCIICBESA-N"), "MXXWOMGUGJBKIW-YPCIICBESA-N");
        assert_eq!(name("CE10; MXXWOMGUGJBKIW-YPCIICBESA-N"), "CE10");
    }
}
