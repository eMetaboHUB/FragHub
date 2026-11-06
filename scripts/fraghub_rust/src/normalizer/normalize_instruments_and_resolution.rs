// src/normalizer/normalize_instruments_and_resolution.rs
use std::collections::HashMap;
use regex::Regex;
use once_cell::sync::Lazy;

/// Nettoie les informations de base des instruments.
///
/// Pour un développeur Python : L'utilisation de `Lazy::new()` permet de compiler les regex
/// une seule fois au premier appel (équivalent à compiler ses regex globalement en Python),
/// évitant de tuer les performances en les recompilant à chaque passage.
fn apply_common_cleaning(mut s: String) -> String {
    // "q-tof", "it-tof", "qit-tof" -> un seul mot. (Avant : tout "-tof" était collé, et dans un type
    // "GC-EI-TOF" cela donnait "eitof" : ni l'ionisation EI ni l'analyseur TOF n'étaient plus reconnus.)
    static ANALYZER_TOF: Lazy<Regex> = Lazy::new(|| Regex::new(r"\b(q|it|qit)[- ]tof").unwrap());
    s = ANALYZER_TOF.replace_all(&s, "${1}tof").into_owned();
    s = s.replace("q-", "q");
    s = s.replace("q exactive", " qexactive ");
    s = s.replace("applied biosystems", " sciex ");
    s = s.replace(" ab ", " sciex ");
    s = s.replace("sciex", " sciex ");

    // TripleTOF (Sciex) est un Q-TOF, pas un triple quadripôle (avant : classé "qqq")
    static TRIPLE_TOF: Lazy<Regex> = Lazy::new(|| Regex::new(r"triple[- ]?tof").unwrap());
    s = TRIPLE_TOF.replace_all(&s, " qtof ").into_owned();

    static TRIPLE_QUAD: Lazy<Regex> = Lazy::new(|| Regex::new(r"triple[- ]?quad").unwrap());
    s = TRIPLE_QUAD.replace_all(&s, " qqq ").into_owned();

    static UPLC: Lazy<Regex> = Lazy::new(|| Regex::new(r".{3} uplc .{3}").unwrap());
    s = UPLC.replace_all(&s, " ").into_owned();

    s
}

/// Nettoie un champ : minuscules, nettoyage commun PUIS tirets remplacés par des espaces.
/// (Avant, les tirets étaient remplacés en premier : "Q-TOF" devenait "q tof" et n'était jamais reconnu.)
fn clean_field(value: &str) -> String {
    apply_common_cleaning(value.to_lowercase()).replace('-', " ")
}

/// Texte de recherche : ". mots ." (caractères spéciaux retirés, espaces simples).
fn search_text(parts: &[String]) -> String {
    static CLEAN_CHARS: Lazy<Regex> = Lazy::new(|| Regex::new(r"[^-\w\s]").unwrap());
    let cleaned = CLEAN_CHARS.replace_all(&parts.join(" "), " ").into_owned();
    format!(". {} .", cleaned.split_whitespace().collect::<Vec<&str>>().join(" "))
}

/// Texte dans lequel on cherche l'instrument : INSTRUMENT et INSTRUMENTTYPE.
/// COMMENT n'est utilisé que si ces deux champs sont vides : sinon, une date ("2020") ou le nom
/// de la chromatographie ("Waters Acquity") pouvaient être pris pour le modèle ou la marque.
fn clean_spectrum_instrument_info(metadata_dict: &HashMap<String, String>) -> String {
    let get = |k: &str| metadata_dict.get(k).map(|s| s.as_str()).unwrap_or("");
    let instrument = clean_field(get("INSTRUMENT"));
    let instrument_type = clean_field(get("INSTRUMENTTYPE"));
    if instrument.trim().is_empty() && instrument_type.trim().is_empty() {
        search_text(&[clean_field(get("COMMENT"))])
    } else {
        search_text(&[instrument, instrument_type])
    }
}

// Remplace la lenteur du regex rf"(\b|^|$){key}(\b|^|$)" de Python par une vérification native ultra-rapide
// Python utilise un moteur Regex basé sur le langage C, mais les appels répétés en Python coûtent cher.
// En Rust, analyser directement les octets (`bytes`) est immensément plus rapide.
fn is_word_char(b: u8) -> bool { b.is_ascii_alphanumeric() || b == b'_' }

fn contains_word(text: &str, word: &str) -> bool {
    if word.is_empty() { return false; }
    let text_bytes = text.as_bytes(); // Travaille directement sur les octets ASCII
    let word_bytes = word.as_bytes();
    let word_len = word_bytes.len();

    let mut start = 0;
    while let Some(pos) = text[start..].find(word) {
        let actual_pos = start + pos;
        let end_pos = actual_pos + word_len;

        // Vérifie les limites du mot (boundaries)
        let is_left_boundary = actual_pos == 0 || !is_word_char(text_bytes[actual_pos - 1]);
        let is_right_boundary = end_pos == text_bytes.len() || !is_word_char(text_bytes[end_pos]);

        if is_left_boundary && is_right_boundary { return true; }
        start = actual_pos + 1;
    }
    false
}

/// Parcourt l'arborescence (arbre JSON de configuration) pour résoudre l'instrument.
/// La syntaxe `<'a>` décrit une "Lifetime" (durée de vie). Elle indique au compilateur que les 
/// références retournées par la fonction vivront aussi longtemps que l'arbre `tree` passé en paramètre.
/// Si aucune clé du niveau n'est dans le texte, on prend la branche "not found" (information absente),
/// comme le prévoit l'arbre. Avant, la recherche s'arrêtait : 96,6 % des spectres ne trouvaient rien.
fn search_level<'a>(tree: &'a serde_json::Value, infos: &str) -> Option<(&'a String, &'a serde_json::Value)> {
    let map = tree.as_object()?;
    for (key, next_level) in map {
        if contains_word(infos, key) { return Some((key, next_level)); }
    }
    map.get_key_value(NOT_FOUND)
}

const NOT_FOUND: &str = "not found";

/// Éléments d'un type d'instrument que l'arbre ne doit pas faire disparaître.
const TYPE_TOKENS: &[&str] = &[
    "lc", "gc", "ce",
    "esi", "apci", "appi", "ei", "ci", "fab", "fd", "maldi", "ptr",
    "qtof", "qft", "itft", "tof", "qqq", "qq", "q", "it", "ittof", "qittof", "qit", "b", "eb", "ebeb", "ft",
];

/// Le type proposé par l'arbre contient-il tous les éléments connus du type écrit dans le fichier ?
fn keeps_type_tokens(raw_type: &str, new_type: &str) -> bool {
    let tokens = |t: &str| -> Vec<String> {
        clean_field(t).split_whitespace().map(|w| w.to_string()).filter(|w| TYPE_TOKENS.contains(&w.as_str())).collect()
    };
    let new_tokens = tokens(new_type);
    tokens(raw_type).iter().all(|t| new_tokens.contains(t))
}

/// Une valeur de l'arbre qui ne dit rien ("NOT FOUND") ne remplace jamais la valeur du fichier.
fn is_known(value: &str) -> bool {
    !value.is_empty() && !value.eq_ignore_ascii_case(NOT_FOUND)
}

/// Identifie avec précision le modèle d'instrument, le type d'ionisation et la résolution.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
/// * `context` (&super::NormalizerContext) : Contexte contenant l'arbre de décision JSON des instruments.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire mis à jour avec INSTRUMENT, INSTRUMENTTYPE, RESOLUTION.
pub fn normalize_instruments_and_resolution(mut metadata_dict: HashMap<String, String>, context: &super::NormalizerContext) -> HashMap<String, String> {
    let infos = clean_spectrum_instrument_info(&metadata_dict);

    let mut current_level = &context.instrument_tree;
    let mut tree_levels = Vec::new();

    // On descend les 5 niveaux (Marque, Modèle, Type Spectre, Type Instrument, Ionisation)
    for _ in 0..5 {
        if let Some((key, next)) = search_level(current_level, &infos) {
            tree_levels.push(key.clone());
            current_level = next;
        } else {
            return metadata_dict; // Niveau sans clé correspondante ni branche "not found" : dico intact
        }
    }

    if tree_levels.len() == 5 {
        metadata_dict.insert("TREE_MARQUE".to_string(), tree_levels[0].clone());
        metadata_dict.insert("TREE_MODELE".to_string(), tree_levels[1].clone());
        metadata_dict.insert("TREE_SPECTYPE".to_string(), tree_levels[2].clone());
        metadata_dict.insert("TREE_INSTYPE".to_string(), tree_levels[3].clone());
        metadata_dict.insert("TREE_IONI".to_string(), tree_levels[4].clone());
    }

    // Niveau 6 : Résolution
    // L'arbre complète les champs ; il ne remplace jamais une valeur du fichier par une valeur moins précise.
    if let Some(res_level) = current_level.as_object() {
        let resolution = if res_level.contains_key("high") { "high" }
        else if res_level.contains_key("low") { "low" }
        else { NOT_FOUND };

        // Chaînage de méthodes sécurisé via `and_then` pour naviguer dans l'objet JSON.
        if let Some(solution_str) = res_level.get(resolution).and_then(|r| r.get("SOLUTION")).and_then(|s| s.as_str()) {
            let parts: Vec<&str> = solution_str.split(',').map(|p| p.trim()).collect();
            if parts.len() >= 3 {
                let raw_instrument = metadata_dict.get("INSTRUMENT").cloned().unwrap_or_default();
                let raw_type = metadata_dict.get("INSTRUMENTTYPE").cloned().unwrap_or_default();
                let raw_res = metadata_dict.get("RESOLUTION").cloned().unwrap_or_default().to_lowercase();

                // INSTRUMENT : une solution générique ("Agilent instrument") ne remplace pas un nom écrit
                // (ex. "Pegasus III TOF-MS system, Leco; GC 6890, Agilent" ou "Bruker timsTOF Pro 2").
                let generic = parts[0].to_lowercase().ends_with(" instrument");
                if is_known(parts[0]) && (raw_instrument.trim().is_empty() || !generic) {
                    metadata_dict.insert("INSTRUMENT".to_string(), parts[0].to_string());
                }

                // INSTRUMENTTYPE : remplacé seulement si le type de l'arbre garde tous les éléments connus
                // du type écrit (LC/GC, ionisation, analyseur). Sinon "LC-ESI-QQ" devenait "LC-ESI",
                // "GC-EI-TOF" devenait "GC", "EI-B" perdait son EI (utilisé pour reconnaître le GC).
                if is_known(parts[1]) && keeps_type_tokens(&raw_type, parts[1]) {
                    metadata_dict.insert("INSTRUMENTTYPE".to_string(), parts[1].to_string());
                }

                // RESOLUTION : une résolution écrite comme nominale / basse n'est pas remplacée
                let says_low = raw_res.contains("nominal") || raw_res.contains("low") || raw_res.contains("unit");
                if is_known(parts[2]) && !says_low {
                    metadata_dict.insert("RESOLUTION".to_string(), parts[2].to_string());
                }
            }
        }
    }

    // IONIZATION : celle du fichier si elle est connue ; sinon celle du niveau « ionisation » de l'arbre ;
    // sinon celle écrite dans le type final (ex. "LC-ESI-QTOF" -> ESI).
    // (Avant : 2e morceau du type de l'arbre, faux pour les types en deux parties : "ESI-QFT" donnait "QFT".)
    let current_ionization = metadata_dict.get("IONIZATION").cloned().unwrap_or_default();
    if current_ionization.trim().is_empty() {
        let from_type = metadata_dict.get("INSTRUMENTTYPE")
            .and_then(|t| crate::globals_vars::IONIZATION_MODE_PATTERN.captures(t))
            .and_then(|c| c.get(1)).map(|m| m.as_str().to_uppercase());
        let ionization = if is_known(&tree_levels[4]) { Some(tree_levels[4].to_uppercase()) } else { from_type };
        if let Some(ion) = ionization {
            metadata_dict.insert("IONIZATION".to_string(), if ion == "ACPI" { "APCI".to_string() } else { ion });
        }
    }

    metadata_dict
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Extrait du vrai arbre (datas/instruments_tree.json), mêmes clés et mêmes solutions.
    fn context() -> super::super::NormalizerContext {
        let tree: serde_json::Value = serde_json::from_str(r#"{
          "not found": {
            "orbitrap": {"not found": {"not found": {"not found": {"high": {"SOLUTION": "Thermo Fisher Scientific Orbitrap, LC-FT, HIGH"}}}}},
            "not found": {
              "not found": {"qtof": {"not found": {"high": {"SOLUTION": "NOT FOUND, QTOF, HIGH"}}}},
              "lc": {"qtof": {"esi": {"high": {"SOLUTION": "NOT FOUND, LC-ESI-QTOF, HIGH"}}}}
            }
          },
          "hitachi": {"not found": {"not found": {"b": {"not found": {"low": {"SOLUTION": "Hitachi instrument, B, LOW"}}}}}},
          "thermo": {"qexactive": {"not found": {"qft": {"esi": {"high": {"SOLUTION": "Thermo Fisher Scientific Q-Exactive, ESI-QFT, HIGH"}}}}}},
          "bruker": {"not found": {"not found": {"qtof": {"not found": {"high": {"SOLUTION": "Bruker instrument, QTOF, HIGH"}}}}}}
        }"#).unwrap();
        super::super::NormalizerContext {
            adduct_pos: HashMap::new(), adduct_neg: HashMap::new(),
            adduct_massdiff_pos: HashMap::new(), adduct_massdiff_neg: HashMap::new(),
            instrument_tree: tree,
        }
    }

    fn run(fields: &[(&str, &str)]) -> HashMap<String, String> {
        let m: HashMap<String, String> = fields.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect();
        normalize_instruments_and_resolution(m, &context())
    }

    #[test]
    fn branche_not_found_quand_l_information_manque() {
        // GNPS : "Orbitrap" seul, ni marque ni type (avant : aucune correspondance)
        let m = run(&[("INSTRUMENT", "Orbitrap"), ("INSTRUMENTTYPE", "")]);
        assert_eq!(m["INSTRUMENT"], "Thermo Fisher Scientific Orbitrap");
        assert_eq!(m["INSTRUMENTTYPE"], "LC-FT");
        assert_eq!(m["RESOLUTION"], "HIGH");
        // La solution "NOT FOUND" ne remplace pas la valeur du fichier
        let m = run(&[("INSTRUMENT", ""), ("INSTRUMENTTYPE", "qtof")]);
        assert_eq!(m["INSTRUMENT"], "");
        assert_eq!(m["INSTRUMENTTYPE"], "QTOF");
    }

    #[test]
    fn q_tof_avec_tiret_reconnu() {
        // MoNA : "Bruker timsTOF Pro 2" / "Q-TOF" (avant : "q tof", jamais reconnu)
        let m = run(&[("INSTRUMENT", "Bruker timsTOF Pro 2"), ("INSTRUMENTTYPE", "Q-TOF")]);
        assert_eq!(m["INSTRUMENTTYPE"], "QTOF");
        let m = run(&[("INSTRUMENTTYPE", "LC-ESI-Q-TOF")]);
        assert_eq!(m["INSTRUMENTTYPE"], "LC-ESI-QTOF");
        assert_eq!(m["IONIZATION"], "ESI");
    }

    #[test]
    fn mention_ei_jamais_effacee() {
        // Hitachi "EI-B" : la solution de l'arbre dit "B" ; on garde "EI-B" pour reconnaître le GC
        let m = run(&[("INSTRUMENT", "HITACHI M-80B"), ("INSTRUMENTTYPE", "EI-B")]);
        assert_eq!(m["INSTRUMENTTYPE"], "EI-B");
        assert_eq!(m["INSTRUMENT"], "HITACHI M-80B");
        assert_eq!(m["RESOLUTION"], "LOW");
    }

    #[test]
    fn jamais_moins_precis_que_le_fichier() {
        // Type écrit plus complet que la solution de l'arbre : gardé
        assert!(!keeps_type_tokens("LC-ESI-QQ", "LC-ESI"));
        assert!(!keeps_type_tokens("GC-EI-TOF", "GC"));
        assert!(!keeps_type_tokens("LC-ESI-QFT", "LC-ESI-ITFT"));
        assert!(!keeps_type_tokens("EI-B", "GC-B"));
        // Type écrit sans élément connu, ou contenu dans la solution : remplacé
        assert!(keeps_type_tokens("Orbitrap", "LC-FT"));
        assert!(keeps_type_tokens("", "QTOF"));
        assert!(keeps_type_tokens("Q-TOF", "QTOF"));
        assert!(keeps_type_tokens("LC-ESI-MS", "LC-ESI-FT"));
        // Nom d'instrument écrit : la solution générique ne le remplace pas
        let m = run(&[("INSTRUMENT", "Bruker timsTOF Pro 2"), ("INSTRUMENTTYPE", "Q-TOF")]);
        assert_eq!(m["INSTRUMENT"], "Bruker timsTOF Pro 2");
        // Résolution nominale : gardée
        let m = run(&[("INSTRUMENT", "Orbitrap"), ("RESOLUTION", "nominal")]);
        assert_eq!(m["RESOLUTION"], "nominal");
    }

    #[test]
    fn ionisation_prise_au_niveau_ionisation() {
        // Avant : 2e morceau du type, "ESI-QFT" donnait IONIZATION = "QFT"
        let m = run(&[("INSTRUMENT", "Thermo Q Exactive HF"), ("INSTRUMENTTYPE", "ESI-QFT"), ("IONIZATION", "ESI")]);
        assert_eq!(m["INSTRUMENTTYPE"], "ESI-QFT");
        assert_eq!(m["IONIZATION"], "ESI");
        // IONIZATION du fichier connue : jamais remplacée
        let m = run(&[("INSTRUMENT", "Orbitrap"), ("IONIZATION", "APCI")]);
        assert_eq!(m["IONIZATION"], "APCI");
        // IONIZATION vide : prise dans le type final
        let m = run(&[("INSTRUMENTTYPE", "LC-ESI-Q-TOF"), ("IONIZATION", "")]);
        assert_eq!(m["IONIZATION"], "ESI");
    }

    #[test]
    fn commentaire_seulement_si_instrument_et_type_vides() {
        let text = clean_spectrum_instrument_info(&[("INSTRUMENT", "Orbitrap"), ("COMMENT", "acquired 2020 on Waters Acquity")]
            .iter().map(|(k, v)| (k.to_string(), v.to_string())).collect());
        assert_eq!(text, ". orbitrap .");
        let text = clean_spectrum_instrument_info(&[("COMMENT", "instrument=Orbitrap")]
            .iter().map(|(k, v)| (k.to_string(), v.to_string())).collect());
        assert_eq!(text, ". instrument orbitrap .");
    }

    #[test]
    fn tripletof_est_un_qtof() {
        assert!(clean_field("TripleTOF 6600").contains("qtof"));
        assert!(clean_field("Triple Quad 5500").contains("qqq"));
        assert_eq!(clean_field("GC-EI-TOF"), "gc ei tof");
        assert_eq!(clean_field("LC-ESI-Q-TOF"), "lc esi qtof");
        assert_eq!(clean_field("LC-ESI-IT-TOF"), "lc esi ittof");
    }
}

