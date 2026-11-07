// src/normalizer/values_normalizer.rs
use std::collections::HashMap;

/// Point d'entrée principal (orchestrateur) pour toutes les étapes de normalisation.
/// 
/// Pour un développeur Python : C'est ici que l'on chaîne l'ensemble de notre pipeline.
/// Remarquez l'usage de `Option<HashMap>`. Si une étape (comme `delete_no_smiles_no_inchi` 
/// ou `check_for_bad_adduct`) décide que le spectre est invalide, elle renvoie `None`. 
/// L'orchestrateur attrape ce `None` via `.is_none()`, et stoppe la chaîne immédiatement (équivalent à `return None`).
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Dictionnaire contenant les métadonnées du spectre.
/// * `deletion_reason` (&mut Option<String>) : Pointeur mutable pour stocker la raison de suppression.
/// * `context` (&NormalizerContext) : Contexte global en lecture seule (dictionnaires chargés en RAM).
///
/// # Returns
/// * `Option<HashMap<String, String>>` : Le dictionnaire normalisé, ou `None` si le spectre a été rejeté.
pub fn normalize_values(
    mut metadata_dict: HashMap<String, String>,
    deletion_reason: &mut Option<String>,
    context: &super::NormalizerContext
) -> Option<HashMap<String, String>> {

    // 1. Nettoie les valeurs "null", "nan"
    metadata_dict = super::normalize_empties::normalize_empties(metadata_dict);
    
    // 2. Corrige les inversions SMILES/INCHI
    metadata_dict = super::repair_mol_descriptors::repair_mol_descriptors(metadata_dict);

    // 3. Supprime si on n'a ni SMILES ni INCHI ni INCHIKEY
    let metadata_opt = super::delete_no_smiles_no_inchi::delete_no_smiles_no_inchi_no_inchikey(metadata_dict, deletion_reason);
    if metadata_opt.is_none() { return None; } // Interruption prématurée si le spectre est mauvais
    metadata_dict = metadata_opt.unwrap();

    // 4. Standardise le statut in-silico (prédit). AVANT la normalisation des instruments, pour lire INSTRUMENTTYPE
    //    tel qu'il est dans le fichier (ex. LipidBlast : "in-silico QTOF"). Les étapes suivantes ne lisent ni ne
    //    modifient COMMENT, FILENAME et PREDICTED ; seule l'étape 8 modifie NAME (fin du nom : énergie, InChIKey, adduit).
    metadata_dict = super::normalize_predicted::normalize_predicted(metadata_dict);
    // Liste interne de toutes les lignes "Synon:", lue seulement par normalize_predicted
    metadata_dict.remove(crate::globals_vars::ALL_SYNONS_KEY);

    // 5. Normalise la méthode d'ionisation
    metadata_dict = super::normalize_ionization::normalize_ionization(metadata_dict);

    // 6. Détermine et normalise l'instrument et la résolution à partir des dictionnaires JSON en RAM
    metadata_dict = super::normalize_instruments_and_resolution::normalize_instruments_and_resolution(metadata_dict, context);

    // 7. Standardise le mode d'ionisation ("positive" / "negative"). AVANT l'adduit : un adduit inconnu
    //    écrit sans charge la prend du mode d'ionisation. Ni l'étape 5 ni l'étape 6 ne lisent IONMODE.
    metadata_dict = super::normalize_ionmode::normalize_ion_mode(metadata_dict);

    // 8. Nom : retire l'énergie et l'InChIKey écrites à la fin (MS-DIAL : "Piperine; CE0; MXXW…-N").
    //    Adduit vide : lu à la fin du nom s'il y est écrit (MGF GNPS : "Hectochlorin M+Na"), puis
    //    normalisation du type de précurseur (adduct)
    metadata_dict = super::normalize_name::clean_name_suffixes(metadata_dict);
    metadata_dict = super::normalize_adduct::adduct_from_name_if_missing(metadata_dict, context);
    metadata_dict = super::normalize_adduct::normalize_adduct(metadata_dict, context);

    // 9. Mode d'ionisation absent ou illisible : déduit du signe de l'adduit ("[M+H]+" -> positive)
    metadata_dict = super::normalize_ionmode::ion_mode_from_adduct_if_missing(metadata_dict);

    // 10. Complète l'adduit (in silico sans adduit, forme courte "M") AVANT le recalcul du m/z, qui en a besoin
    metadata_dict = super::check_for_bad_adduct::complete_adduct(metadata_dict);

    // 11. Recalcule le PRECURSORMZ manquant ou invalide avec RDKit (C++)
    metadata_dict = super::missing_precursormz_re_calculation::missing_precursormz_re_calculation(metadata_dict, context);

    // 12. Valide la cohérence Adduct / Mode Ionisation, et la présence du mode ; supprime sinon
    let metadata_opt_2 = super::check_for_bad_adduct::check_for_bad_adduct(metadata_dict, deletion_reason, context);
    if metadata_opt_2.is_none() {
        return None;
    }
    metadata_dict = metadata_opt_2.unwrap();

    // 13. Normalise le niveau MS (ex: "MS2")
    metadata_dict = super::normalize_ms_level::normalize_ms_level(metadata_dict);

    // 14. Normalise les unités de temps de rétention en minutes
    metadata_dict = super::normalize_retentiontime::normalize_retention_time(metadata_dict);

    // On enveloppe le dictionnaire valide et entièrement normalisé dans `Some` pour le renvoyer.
    Some(metadata_dict)
}