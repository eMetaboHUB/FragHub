use std::collections::HashMap;
use pyo3::prelude::*;
use regex::Regex;
use once_cell::sync::Lazy;

static IE_EI_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)\b(IE|EI)\b").unwrap());

/// Appel direct au RDKit de Python depuis Rust !
///
/// Pour un développeur Python : C'est ici qu'intervient la magie de `pyo3`.
/// La fonction `Python::with_gil` acquiert le Global Interpreter Lock (GIL) de Python,
/// ce qui permet à Rust d'importer et d'utiliser la librairie `rdkit` comme si c'était du Python natif.
///
/// # Arguments
/// * `inchi` (&str) : Chaîne de caractères InChI.
/// * `smiles` (&str) : Chaîne de caractères SMILES.
///
/// # Returns
/// * `Option<f64>` : Masse exacte calculée par RDKit, ou `None` si échec.
fn get_exact_mass_from_rdkit(inchi: &str, smiles: &str) -> Option<f64> {
    // Par défaut : RDKit C++ (src/rdkit_native), sans verrou Python, donc en parallèle sur tous les cœurs.
    // FRAGHUB_RDKIT_BACKEND=python : ancienne version ci-dessous (RDKit Python), même résultat.
    if !crate::rdkit_native::use_python_backend() {
        return crate::rdkit_native::exact_mass(inchi, smiles);
    }
    Python::with_gil(|py| {
        // Désactiver les logs rdkit comme en Python
        if let Ok(rdlogger) = py.import_bound("rdkit.RDLogger") {
            let _ = rdlogger.call_method1("DisableLog", ("rdApp.*",)); // L'équivalent de rdLogger.DisableLog('rdApp.*')
        }

        // Importation dynamique des modules
        let chem = py.import_bound("rdkit.Chem").ok()?;
        let descriptors = py.import_bound("rdkit.Chem.Descriptors").ok()?;

        // Instanciation de la molécule (MolFromInchi ou MolFromSmiles)
        let mol = if !inchi.is_empty() && inchi.contains("InChI=") {
            chem.call_method1("MolFromInchi", (inchi,)).ok()
        } else if !smiles.is_empty() {
            chem.call_method1("MolFromSmiles", (smiles,)).ok()
        } else {
            None
        };

        if let Some(mol) = mol {
            if !mol.is_none() {
                // Calcul de la masse exacte (ExactMolWt)
                if let Ok(mass_obj) = descriptors.call_method1("ExactMolWt", (&mol,)) {
                    // `.extract::<f64>()` convertit l'objet Python renvoyé en un `float` Rust !
                    return mass_obj.extract::<f64>().ok();
                }
            }
        }
        None
    })
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