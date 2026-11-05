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
    //    tel qu'il est dans le fichier (ex. LipidBlast : "in-silico QTOF"). Les étapes 5 à 9 ne lisent ni ne modifient
    //    COMMENT, NAME, FILENAME et PREDICTED : déplacer cette étape ne change rien d'autre.
    metadata_dict = super::normalize_predicted::normalize_predicted(metadata_dict);

    // 5. Normalise la méthode d'ionisation
    metadata_dict = super::normalize_ionization::normalize_ionization(metadata_dict);

    // 6. Détermine et normalise l'instrument et la résolution à partir des dictionnaires JSON en RAM
    metadata_dict = super::normalize_instruments_and_resolution::normalize_instruments_and_resolution(metadata_dict, context);

    // 7. Normalise le type de précurseur (adduct)
    metadata_dict = super::normalize_adduct::normalize_adduct(metadata_dict, context);

    // 8. Recalculate missing or invalid PRECURSORMZ via RDKit (en passant par pyo3)
    metadata_dict = super::missing_precursormz_re_calculation::missing_precursormz_re_calculation(metadata_dict, context);

    // 9. Standardise le mode d'ionisation ("positive" / "negative")
    metadata_dict = super::normalize_ionmode::normalize_ion_mode(metadata_dict);

    // 10. Valide la cohérence Adduct / Mode Ionisation ; supprime si incohérent
    let metadata_opt_2 = super::check_for_bad_adduct::check_for_bad_adduct(metadata_dict, deletion_reason, context);
    if metadata_opt_2.is_none() {
        return None;
    }
    metadata_dict = metadata_opt_2.unwrap();

    // 11. Normalise le niveau MS (ex: "MS2")
    metadata_dict = super::normalize_ms_level::normalize_ms_level(metadata_dict);

    // 12. Normalise les unités de temps de rétention en minutes
    metadata_dict = super::normalize_retentiontime::normalize_retention_time(metadata_dict);

    // On enveloppe le dictionnaire valide et entièrement normalisé dans `Some` pour le renvoyer.
    Some(metadata_dict)
}