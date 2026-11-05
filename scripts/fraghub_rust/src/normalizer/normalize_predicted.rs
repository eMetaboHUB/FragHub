use std::collections::HashMap;

/// Détecte si le spectre est "in silico" (prédit par ordinateur).
///
/// Pour un développeur Python : La logique if/else if est stricte en Rust.
/// On doit comparer des `String` en utilisant `==` et `.to_lowercase()`.
///
/// # Arguments
/// * `metadata_dict` (HashMap<String, String>) : Le dictionnaire des métadonnées du spectre.
///
/// # Returns
/// * `HashMap<String, String>` : Le dictionnaire avec PREDICTED à "true" ou "false".
pub fn normalize_predicted(mut metadata_dict: HashMap<String, String>) -> HashMap<String, String> {
    let predicted = metadata_dict.get("PREDICTED").cloned().unwrap_or_default();
    if predicted.to_lowercase() == "false" {
        return metadata_dict;
    }

    let comment = metadata_dict.get("COMMENT").cloned().unwrap_or_default();
    let filename = metadata_dict.get("FILENAME").cloned().unwrap_or_default();
    let name = metadata_dict.get("NAME").cloned().unwrap_or_default();
    let instrument_type = metadata_dict.get("INSTRUMENTTYPE").cloned().unwrap_or_default();
    // Toutes les lignes "Synon:" du fichier (SYNON n'en garde que la dernière), voir raw_metadata.rs
    let all_synons = metadata_dict.get(crate::globals_vars::ALL_SYNONS_KEY).cloned().unwrap_or_default();

    // On retire d'abord les formulations qui ne concernent pas le spectre ("predicted rt", "theoretical m/z"…)
    let comment_clean = crate::globals_vars::IN_SILICO_EXCLUDE_PATTERN.replace_all(&comment, " ");

    let mut is_predicted = false;

    // Le `||` est le OR logique (équivalent à `or` en Python).
    if crate::globals_vars::IN_SILICO_PATTERN.is_match(&comment_clean) || predicted.to_lowercase() == "true" {
        is_predicted = true;
    } else if crate::globals_vars::IN_SILICO_STRONG_PATTERN.is_match(&instrument_type) {
        // Type d'instrument valeur brute (cette étape passe avant la normalisation des instruments)
        is_predicted = true;
    } else if crate::globals_vars::IN_SILICO_STRONG_PATTERN.is_match(&all_synons) {
        // MoNA : "Synon: $:07 In-Silico-Spectrum", souvent suivie d'une autre ligne Synon
        is_predicted = true;
    } else if !filename.contains("MSMS_Public") {
        // `format!` est l'équivalent parfait des f-strings Python : f"{filename} {name}"
        let combined = format!("{} {}", filename, name);
        let combined_clean = crate::globals_vars::IN_SILICO_EXCLUDE_PATTERN.replace_all(&combined, " ");
        if crate::globals_vars::IN_SILICO_PATTERN.is_match(&combined_clean) {
            is_predicted = true;
        }
    }

    // L'expression `if condition { A } else { B }` peut être utilisée pour retourner une valeur directement.
    // C'est l'équivalent de `A if condition else B` en Python.
    metadata_dict.insert("PREDICTED".to_string(), if is_predicted { "true".to_string() } else { "false".to_string() });

    metadata_dict
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(fields: &[(&str, &str)]) -> String {
        let m: HashMap<String, String> = fields.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect();
        normalize_predicted(m)["PREDICTED"].clone()
    }

    // --- Spectres expérimentaux qui étaient classés in silico à tort ---
    #[test]
    fn predicted_retention_time_is_not_in_silico() {
        assert_eq!(run(&[("COMMENT", r#""SMILES=CCO" "predicted rt rp=5.12" "msn observed rt=5.3""#)]), "false");
        assert_eq!(run(&[("COMMENT", r#""dataset id=MSV000" "predicted rt hilic=3.4""#)]), "false");
        assert_eq!(run(&[("COMMENT", "Predicted retention time: 12.3 min")]), "false");
    }

    #[test]
    fn theoretical_mz_of_annotation_is_not_in_silico() {
        assert_eq!(run(&[("COMMENT", r#""annotation [C6H5-H]+, rule of HR True=1" "theoretical m/z 77.0386, mass diff 0.0002 (2.6 ppm)""#)]), "false");
    }

    #[test]
    fn massbank_predicted_annotation_is_not_in_silico() {
        assert_eq!(run(&[("COMMENT", "DB#=MSBNK-mFam-MC000003; origin=MassBank; Annotation Predicted")]), "false");
        assert_eq!(run(&[("COMMENT", "DB#=MSBNK-mFam-MC000003; origin=MassBank; Annotation Predicted; ID 12")]), "false");
        assert_eq!(run(&[("COMMENT", "DB#=MSBNK-mFam-MC000010; origin=MassBank; Predicted")]), "false");
    }

    #[test]
    fn sirius_formula_and_annotation_level_are_not_in_silico() {
        assert_eq!(run(&[("COMMENT", "predicted molecular formula SIRIUS: C10H12O")]), "false");
        assert_eq!(run(&[("COMMENT", "DB#=MSBNK-RIKEN-PR000001; origin=MassBank; Annotation level-3")]), "false");
    }

    // --- Vrais spectres in silico : toujours reconnus ---
    #[test]
    fn real_predicted_spectra_are_still_in_silico() {
        assert_eq!(run(&[("COMMENT", r#""origin=Predicted_GC-MS_Spectra_MoNA.mgf""#)]), "true");
        assert_eq!(run(&[("COMMENT", r#""origin=in-silico.msp""#)]), "true");
        assert_eq!(run(&[("COMMENT", "theoretical spectrum")]), "true");
        assert_eq!(run(&[("PREDICTED", "true")]), "true");
        assert_eq!(run(&[("FILENAME", "UNPD_ISDB_R_p01_insilico.mgf"), ("NAME", "x")]), "true");
    }

    #[test]
    fn in_silico_instrument_type_is_detected() {
        // LipidBlast : "Instrument_type: in-silico QTOF"
        assert_eq!(run(&[("INSTRUMENTTYPE", "in-silico QTOF"), ("FILENAME", "LipidBlast.msp")]), "true");
        // "predicted" seul dans INSTRUMENTTYPE ne suffit pas
        assert_eq!(run(&[("INSTRUMENTTYPE", "LC-ESI-QTOF"), ("FILENAME", "a.msp")]), "false");
    }

    #[test]
    fn spectrum_prediction_tools_are_in_silico() {
        for c in ["predicted with CFM-ID 4.0", "cfmid", "CFM ID", "theoretical spectrum (QCEIMS)", "QCxMS",
                  "NEIMS model", "MassFormer", "ICEBERG", "SCARF", "3DMolMS", "GRAFF-MS", "FraGNNet", "LipidBlast"] {
            assert_eq!(run(&[("COMMENT", c)]), "true", "{c}");
        }
        assert_eq!(run(&[("FILENAME", "UNPD_ISDB_R_p01.mgf"), ("NAME", "x")]), "true");
        assert_eq!(run(&[("INSTRUMENTTYPE", "QCEIMS"), ("FILENAME", "a.msp")]), "true");
    }

    #[test]
    fn identification_tools_and_lookalikes_are_not_in_silico() {
        // Outils d'identification cités sur des spectres expérimentaux
        assert_eq!(run(&[("COMMENT", "Ranking 1;; PROPOSED_COMPOUND SIRIUS5_Confidence:biodatabase 0.9")]), "false");
        assert_eq!(run(&[("COMMENT", "data transformation=COMMENT MSFinder 3.5 used to identify precursor")]), "false");
        assert_eq!(run(&[("COMMENT", "MS-FINDER:MsfinderDB (Ranking 1), MetFrag:MsfinderDB (Ranking 2)")]), "false");
        assert_eq!(run(&[("NAME", "Candidate Carnitine (delta mass 14.01; sirius atomic difference: 1C, 2H)"), ("FILENAME", "ALL_GNPS.json")]), "false");
        // "ISDB" au milieu d'un InChIKey, "iceberg" en minuscules, "scarface"
        assert_eq!(run(&[("COMMENT", "INCHIKEY=AOCCRKXUBBAOQI-GOSISDBHSA-N")]), "false");
        assert_eq!(run(&[("NAME", "iceberg lettuce extract"), ("FILENAME", "a.msp")]), "false");
        assert_eq!(run(&[("NAME", "SCARFACE"), ("FILENAME", "a.msp")]), "false");
    }

    #[test]
    fn identification_clauses_citing_prediction_tools_are_not_in_silico() {
        // MassBank : CFM-ID utilisé pour IDENTIFIER un spectre mesuré
        assert_eq!(run(&[("COMMENT", "Spectral Match to DL-Arginine from NIST14; Ranking 1;; PROPOSED_COMPOUND CFM-ID:MsfinderDB 1; ODKSFYDXXFIFQN; [M+H]+; L-Arginine")]), "false");
        assert_eq!(run(&[("COMMENT", "PROPOSED_MF SIRIUS5_Zodiac 1; C6H14N4O2; [M+H]+")]), "false");
        assert_eq!(run(&[("COMMENT", "SMB Level 2a; Proposed by SIRIUS5_Confidence (FDR < 0.05); Similar")]), "false");
        assert_eq!(run(&[("COMMENT", "Annotation Tentative candidate (Level 3); SMB Level 3c;  ); Supported by CFM-ID:MsfinderDB (Ranking 1), MS-FINDER:MsfinderDB (Ranking 2); Biological context")]), "false");
        assert_eq!(run(&[("COMMENT", "Similarity by SIMILARITY:GNPS (Dot product > 0.7, 3 spectra); Supported by CFM-ID:MsfinderDB (Ranking 1)")]), "false");
        // MoNA
        assert_eq!(run(&[("COMMENT", r#""annotation=MS-FINDER ver 3.52" "biological matrix=leaf""#)]), "false");
        assert_eq!(run(&[("NAME", "Communesin A (identified trough in silico generated spectrum)"), ("FILENAME", "MoNA-export-Experimental_Spectra.msp")]), "false");
        // Mais un spectre GÉNÉRÉ par CFM-ID reste in silico
        assert_eq!(run(&[("COMMENT", "spectrum predicted with CFM-ID 4.0")]), "true");
        assert_eq!(run(&[("COMMENT", "origin=CFM-ID; energy=10")]), "true");
    }

    #[test]
    fn explicit_false_is_kept() {
        assert_eq!(run(&[("PREDICTED", "false"), ("COMMENT", "in silico")]), "false");
    }

    #[test]
    fn synon_in_silico_reconnu_meme_s_il_n_est_pas_le_dernier() {
        // Cas réel MoNA (DMED-FAHFA) : la ligne "$:07 In-Silico-Spectrum" est suivie de "$:00in-source"
        let sep = crate::globals_vars::INTERNAL_LIST_SEP;
        let synons = ["$:00 ms2", "$:05 30V CID", "$:07 In-Silico-Spectrum", "$:00in-source"].join(sep);
        assert_eq!(run(&[(crate::globals_vars::ALL_SYNONS_KEY, synons.as_str())]), "true");
        // Synonymes ordinaires : pas in silico
        let synons = ["Cinchonain Ia, 6TMS", "$:00in-source"].join(sep);
        assert_eq!(run(&[(crate::globals_vars::ALL_SYNONS_KEY, synons.as_str())]), "false");
        // "predicted" seul dans un synonyme ne suffit pas (motif fort uniquement)
        assert_eq!(run(&[(crate::globals_vars::ALL_SYNONS_KEY, "predicted metabolite")]), "false");
    }
}
