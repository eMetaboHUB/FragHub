use std::collections::HashMap;
use regex::Regex;
use once_cell::sync::Lazy;

static IE_EI_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)\b(IE|EI)\b").unwrap());

/// Masse exacte de la molécule (InChI de préférence, sinon SMILES), calculée par RDKit C++ (src/rdkit_native).
///
/// # Arguments
/// * `inchi` (&str) : Chaîne de caractères InChI.
/// * `smiles` (&str) : Chaîne de caractères SMILES.
///
/// # Returns
/// * `Option<f64>` : Masse exacte calculée par RDKit, ou `None` si échec.
fn get_exact_mass_from_rdkit(inchi: &str, smiles: &str) -> Option<f64> {
    // Sans verrou Python, donc en parallèle sur tous les cœurs.
    crate::rdkit_native::exact_mass(inchi, smiles)
}

/// Recalcule la valeur PRECURSORMZ si elle est manquante en utilisant la masse de RDKit et le delta de l'adduct.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
/// * `context` (&super::NormalizerContext) : Contexte global pour récupérer le delta de masse de l'adduct.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec le `PRECURSORMZ` recalculé si nécessaire.
pub fn missing_precursormz_re_calculation(
    mut metadata_dict: HashMap<String, String>,
    context: &super::NormalizerContext
) -> HashMap<String, String> {

    // 1. Vérifier si on a besoin de recalculer
    let pmz_str = metadata_dict.get("PRECURSORMZ").cloned().unwrap_or_default();
    let mut needs_recalc = false;

    if pmz_str.is_empty() {
        needs_recalc = true;
    } else if let Some(caps) = crate::globals_vars::FLOAT_CHECK_PATTERN.captures(&pmz_str) {
        // En Rust, pas besoin de `try/except ValueError` lourd, on utilise `parse::<f64>()` qui renvoie un `Result`.
        let matched_str = caps.get(1).unwrap().as_str().replace(',', ".");
        if let Ok(val) = matched_str.parse::<f64>() {
            if val <= 0.0 { needs_recalc = true; }
        } else { needs_recalc = true; }
    } else {
        needs_recalc = true;
    }

    if !needs_recalc {
        return metadata_dict; // Pas besoin de calcul, on renvoie directement
    }

    // 2. Extraire les infos pour la différence de masse
    let adduct = metadata_dict.get("PRECURSORTYPE").cloned().unwrap_or_default();
    let inst_type = metadata_dict.get("INSTRUMENTTYPE").cloned().unwrap_or_default();

    if IE_EI_PATTERN.is_match(&inst_type) && adduct.is_empty() {
        return metadata_dict;
    }

    // Chaînage `.or_else()` : cherche dans le dico Positif, et si absent, cherche dans le dico Négatif.
    let mass_diff = context.adduct_massdiff_pos.get(&adduct)
        .or_else(|| context.adduct_massdiff_neg.get(&adduct))
        .copied();

    // 3. Calculer avec RDKit (via PyO3) si une différence de masse est trouvée
    if let Some(diff) = mass_diff {
        let inchi = metadata_dict.get("INCHI").cloned().unwrap_or_default();
        let smiles = metadata_dict.get("SMILES").cloned().unwrap_or_default();

        if let Some(em) = get_exact_mass_from_rdkit(&inchi, &smiles) {
            // m/z = (n·M + Δ) / z : gère les multimères ([2M+H]+) et les ions multichargés ([M+2H]2+).
            let new_pmz = crate::adduct_mass_calculator::precursor_mz_from_exact_mass(em, diff, &adduct);
            metadata_dict.insert("PRECURSORMZ".to_string(), new_pmz.to_string());
        }
    }

    metadata_dict
}