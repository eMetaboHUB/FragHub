// src/normalizer/repair_mol_descriptors.rs
use std::collections::HashMap;

/// Répare la valeur 'INCHI' en s'assurant qu'elle commence par "InChI=".
///
/// Pour un développeur Python : Cette fonction remplace un bloc de texte avec une expression régulière.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec le champ INCHI réparé.
pub fn repair_inchi(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    // `if let Some(inchi)` gère proprement le cas où la clé "INCHI" n'existerait pas dans le dictionnaire.
    // `.get()` retourne un `Option<&String>`, on extrait donc une référence vers le texte.
    if let Some(inchi) = metadata_dict.get("INCHI") {
        if !inchi.is_empty() {
            // Remplacement via regex équivalent au `re.sub` de Python.
            // `.into_owned()` convertit le type `Cow<str>` (Copy-On-Write) retourné par `replace` 
            // en une vraie `String` indépendante. C'est très optimisé en Rust.
            let repaired = crate::globals_vars::REPAIR_INCHI_PATTERN
                .replace(inchi, "InChI=")
                .into_owned();
                
            // `.insert()` remplace l'ancienne valeur si la clé existe déjà.
            // En Rust, "INCHI" est une vue statique (`&str`), on doit appeler `.to_string()` pour l'allouer sur le tas (`String`).
            metadata_dict.insert("INCHI".to_string(), repaired);
        }
    }
    metadata_dict
}

/// Corrige les erreurs de placement entre SMILES, INCHI et INCHIKEY.
/// 
/// Pour un développeur Python : Remplace les logiques de déplacement de champs. Au lieu de copier 
/// `dict["A"] = dict["B"]`, on clone (`.clone()`) car Rust interdit à deux variables d'être propriétaires
/// de la même chaîne de caractères en mémoire.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec les descripteurs moléculaires aux bons endroits.
pub fn repair_mol_descriptors(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    // `.cloned()` crée une vraie copie de la String pour nous permettre de la manipuler.
    // `.unwrap_or_default()` renvoie une chaîne vide `""` si la clé n'existait pas (au lieu de crasher).
    let smiles = metadata_dict.get("SMILES").cloned().unwrap_or_default();
    let inchi = metadata_dict.get("INCHI").cloned().unwrap_or_default();
    let inchikey = metadata_dict.get("INCHIKEY").cloned().unwrap_or_default();

    // --- Vérification initiale ---
    if crate::globals_vars::SMILES_PATTERN.is_match(&smiles) &&
        crate::globals_vars::INCHI_PATTERN.is_match(&inchi) &&
        crate::globals_vars::INCHIKEY_PATTERN.is_match(&inchikey)
    {
        return repair_inchi(metadata_dict); // Sortie prématurée si tout est bon.
    }

    // --- Logique de réparation croisée (Cross-Field Repair) ---
    // En Rust, on modifie directement le dictionnaire via `.insert()` au lieu de passer par l'opérateur `[]` de Python.

    // 1. SMILES trouvé dans le champ INCHI
    // `&inchi` emprunte la variable pour que l'expression régulière la lise sans la détruire.
    if crate::globals_vars::SMILES_PATTERN.is_match(&inchi) &&
        !crate::globals_vars::INCHI_PATTERN.is_match(&inchi) &&
        !crate::globals_vars::INCHIKEY_FULL_PATTERN.is_match(&inchi)
    {
        metadata_dict.insert("SMILES".to_string(), inchi.clone());
        metadata_dict.insert("INCHI".to_string(), String::new());
    }

    // 2. SMILES trouvé dans le champ INCHIKEY
    if crate::globals_vars::SMILES_PATTERN.is_match(&inchikey) &&
        !crate::globals_vars::INCHI_PATTERN.is_match(&inchikey) &&
        !crate::globals_vars::INCHIKEY_PATTERN.is_match(&inchikey)
    {
        metadata_dict.insert("SMILES".to_string(), inchikey.clone());
        metadata_dict.insert("INCHIKEY".to_string(), String::new());
    }

    // 3. INCHI trouvé dans le champ SMILES
    if crate::globals_vars::INCHI_PATTERN.is_match(&smiles) {
        metadata_dict.insert("INCHI".to_string(), smiles.clone());
        metadata_dict.insert("SMILES".to_string(), String::new());
    }

    // 4. INCHI trouvé dans le champ INCHIKEY
    if crate::globals_vars::INCHI_PATTERN.is_match(&inchikey) {
        metadata_dict.insert("INCHI".to_string(), inchikey.clone());
        metadata_dict.insert("INCHIKEY".to_string(), String::new());
    }

    // 5. INCHIKEY trouvé dans le champ INCHI
    if crate::globals_vars::INCHIKEY_FULL_PATTERN.is_match(&inchi) {
        metadata_dict.insert("INCHIKEY".to_string(), inchi.clone());
        metadata_dict.insert("INCHI".to_string(), String::new());
    }

    // 6. INCHIKEY trouvé dans le champ SMILES
    if crate::globals_vars::INCHIKEY_FULL_PATTERN.is_match(&smiles) {
        metadata_dict.insert("INCHIKEY".to_string(), smiles.clone());
        metadata_dict.insert("SMILES".to_string(), String::new());
    }

    // Étape finale : on sécurise le préfixe de l'InChI
    repair_inchi(metadata_dict)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn meta(smiles: &str, inchi: &str, inchikey: &str) -> HashMap<String, String> {
        let mut m = HashMap::new();
        m.insert("SMILES".to_string(), smiles.to_string());
        m.insert("INCHI".to_string(), inchi.to_string());
        m.insert("INCHIKEY".to_string(), inchikey.to_string());
        m
    }

    #[test]
    fn lipid_smiles_is_not_taken_for_an_inchikey() {
        // Acide palmitique : 16 lettres consécutives, InChIKey absent.
        let out = repair_mol_descriptors(meta("CCCCCCCCCCCCCCCC(=O)O", "", ""));
        assert_eq!(out["SMILES"], "CCCCCCCCCCCCCCCC(=O)O");
        assert_eq!(out["INCHIKEY"], "");
    }

    #[test]
    fn pure_carbon_chain_smiles_is_kept() {
        // Tétradécane : exactement 14 "C", forme d'un InChIKey court.
        let out = repair_mol_descriptors(meta("CCCCCCCCCCCCCC", "", ""));
        assert_eq!(out["SMILES"], "CCCCCCCCCCCCCC");
        assert_eq!(out["INCHIKEY"], "");
    }

    #[test]
    fn stereo_smiles_is_not_taken_for_an_inchi() {
        // Acide oléique (double liaison cis notée avec "/C=C\\").
        let smiles = "CCCCCCCC/C=C\\CCCCCCCC(=O)O";
        let out = repair_mol_descriptors(meta(smiles, "", ""));
        assert_eq!(out["SMILES"], smiles);
        assert_eq!(out["INCHI"], "");
        // Autre forme stéréo avec "/C/".
        let smiles2 = "C/C=C/C(=O)O";
        let out2 = repair_mol_descriptors(meta(smiles2, "", ""));
        assert_eq!(out2["SMILES"], smiles2);
        assert_eq!(out2["INCHI"], "");
    }

    #[test]
    fn real_misplacements_are_still_repaired() {
        let inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3";
        let key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N";
        // InChI rangé dans SMILES
        let out = repair_mol_descriptors(meta(inchi, "", key));
        assert_eq!(out["INCHI"], inchi);
        assert_eq!(out["SMILES"], "");
        // InChIKey rangé dans SMILES
        let out = repair_mol_descriptors(meta(key, inchi, ""));
        assert_eq!(out["INCHIKEY"], key);
        assert_eq!(out["SMILES"], "");
        // SMILES rangé dans INCHI
        let out = repair_mol_descriptors(meta("", "CCO(C)CCN", key));
        assert_eq!(out["SMILES"], "CCO(C)CCN");
        assert_eq!(out["INCHI"], "");
    }

    #[test]
    fn inchi_prefix_is_repaired() {
        let out = repair_mol_descriptors(meta("CCO", "1S/C2H6O/c1-2-3/h3H,2H2,1H3", ""));
        assert_eq!(out["INCHI"], "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
        let out = repair_mol_descriptors(meta("CCO", "InChI=InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", ""));
        assert_eq!(out["INCHI"], "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
    }

    #[test]
    fn inchikey_field_validation() {
        use crate::globals_vars::INCHIKEY_PATTERN as P;
        assert!(P.is_match("LFQSCWFLJHTTHZ-UHFFFAOYSA-N"));
        assert!(P.is_match(" LFQSCWFLJHTTHZ-UHFFFAOYSA-N "));
        assert!(P.is_match("LFQSCWFLJHTTHZ"));
        assert!(!P.is_match("CCCCCCCCCCCCCCCC(=O)O"));
        assert!(!P.is_match("NOT FOUND LFQSCWFLJHTTHZ-UHFFFAOYSA-N"));
    }
}
