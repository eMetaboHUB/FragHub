// src/convertors/raw_metadata.rs
use std::collections::HashMap;

/// Champs dont le NOM dit que la valeur est en secondes (format MGF : RTINSECONDS).
/// Sans ce marquage, le renommage RTINSECONDS -> RT (key_to_convert.csv) perdait l'unité, et le nombre seul
/// passait pour des minutes (règle par défaut, qui reste vraie pour tous les autres champs).
const SECONDS_KEYS: &[&str] = &["rtinseconds"];

/// Garde l'unité portée par le nom du champ : "70.99" d'un champ RTINSECONDS devient "70.99 s",
/// que normalize_retention_time convertit ensuite en minutes. Toute autre valeur est rendue telle quelle.
pub fn keep_unit_from_key(raw_key: &str, value: String) -> String {
    if SECONDS_KEYS.contains(&raw_key) {
        let t = value.trim();
        if t.replace(',', ".").parse::<f64>().map_or(false, |x| x.is_finite()) {
            return format!("{t} s");
        }
    }
    value
}

/// Métadonnées brutes d'un spectre (avant renommage vers les clés FragHub), qui retiennent
/// l'ordre d'apparition des champs dans le fichier.
///
/// Pourquoi : plusieurs champs bruts peuvent donner la même clé FragHub
/// (ex. INSTRUMENT_TYPE et SOURCE_INSTRUMENT → INSTRUMENTTYPE, voir key_to_convert.csv).
/// Avant, les champs étaient rangés dans une HashMap, dont l'ordre de parcours est tiré au hasard
/// à chaque lancement : la valeur gardée changeait d'un run à l'autre.
///
/// Règles (identiques pour MSP, MGF et JSON) :
/// - deux champs DIFFÉRENTS qui donnent la même clé FragHub : le premier du fichier gagne, SAUF si sa
///   valeur est vide ou un « vide » écrit (N/A, nan, unknown... : EMPTY_PATTERN, la règle de
///   normalize_empties) : un champ suivant avec une vraie valeur la remplace. Ex. GNPS JSON :
///   INCHI_AUX = "N/A" ne bloque plus InChIKey_smiles ;
/// - le MÊME champ présent plusieurs fois (ex. plusieurs lignes "Synon:") : la dernière valeur
///   gagne, comme avant. Pour "Synon", toutes les valeurs sont en plus gardées dans la clé interne
///   ALL_SYNONS (jamais écrite), lue par la détection in silico ;
/// - nom de fichier : le lecteur met en tête du spectre le nom et l'empreinte (FILEHASH) du fichier
///   d'entrée. FILENAME prend le nom du fichier de la base d'origine si le spectre en donne un (champ
///   FILENAME du spectre, ex. MGF MSnLib "FILENAME=…mzML", MGF GNPS "FILENAME=…mgf", ou un synonyme :
///   GNPS JSON "source_file"...), sinon le nom du fichier d'entrée. Le nom du fichier d'entrée est
///   toujours gardé dans la clé interne INPUT_FILENAME (jamais écrite), lue par les règles qui
///   reconnaissent un fichier à son nom. FILEHASH reste celui du fichier d'entrée ;
/// - énergies de collision numérotées (Spectraverse : COLLISION_ENERGY_1, _2, _3) : réunies en une
///   seule valeur "e1, e2, e3" (valeurs vides ou "nan" ignorées), lue comme COLLISION_ENERGY_1.
#[derive(Default)]
pub struct RawMetadata {
    values: HashMap<String, String>,
    order: Vec<String>,
    synons: Vec<String>,
    /// Premier champ FILENAME rempli venant du spectre lui-même (après celui du lecteur)
    source_filename: Option<String>,
}

impl RawMetadata {
    pub fn with_capacity(n: usize) -> Self {
        RawMetadata { values: HashMap::with_capacity(n), order: Vec::with_capacity(n), synons: Vec::new(), source_filename: None }
    }

    /// Ajoute (ou remplace) un champ. Un champ déjà présent garde sa place d'origine.
    /// FILENAME / FILEHASH : la première valeur (celle du lecteur) est gardée ; un FILENAME venant
    /// ensuite du spectre est mis de côté comme nom du fichier de la base d'origine.
    pub fn insert(&mut self, key: String, value: String) {
        if key == "synon" {
            self.synons.push(value.clone());
        }
        if self.values.contains_key(&key) {
            if key == "filename" && self.source_filename.is_none() && !is_empty_value(&value) {
                self.source_filename = Some(value.trim().to_string());
            }
            if READER_KEYS.contains(&key.as_str()) {
                return;
            }
        } else {
            self.order.push(key.clone());
        }
        self.values.insert(key, value);
    }

    /// Ajoute un champ seulement s'il manque, ou si sa valeur actuelle est vide (EMPTY_PATTERN).
    /// Sert aux valeurs de repli, ex. les valeurs « computed » de MoNA (JSON), qui ne doivent pas
    /// remplacer la valeur d'origine.
    pub fn insert_if_absent(&mut self, key: String, value: String) {
        let absent = self.values.get(&key).map_or(true, |v| is_empty_value(v));
        if absent {
            self.insert(key, value);
        }
    }

    /// Réunit COLLISION_ENERGY_1, _2, _3... (Spectraverse) en une seule valeur dans COLLISION_ENERGY_1.
    fn merge_numbered_collision_energies(&mut self) {
        let mut parts = Vec::new();
        for i in 1..=9 {
            match self.values.get(&format!("collision_energy_{i}")) {
                Some(v) if !is_empty_value(v) => parts.push(v.trim().to_string()),
                Some(_) => {}
                None => break,
            }
        }
        if !parts.is_empty() {
            self.values.insert("collision_energy_1".to_string(), parts.join(", "));
        }
    }

    pub fn get(&self, key: &str) -> Option<&String> {
        self.values.get(key)
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Renomme les champs vers les clés FragHub et les range dans `target`.
    /// `keep` dit si une clé FragHub fait partie des clés à garder.
    /// Si plusieurs champs donnent la même clé, le premier du fichier gagne.
    pub fn map_into<F: Fn(&str) -> bool>(
        mut self,
        keys_dict: &HashMap<String, String>,
        keep: F,
        target: &mut HashMap<String, String>,
    ) {
        if !self.synons.is_empty() && keys_dict.contains_key("synon") {
            let all = self.synons.join(crate::globals_vars::INTERNAL_LIST_SEP);
            target.insert(crate::globals_vars::ALL_SYNONS_KEY.to_string(), all);
        }
        if self.values.contains_key("collision_energy_1") {
            self.merge_numbered_collision_energies();
        }
        // Clés FragHub remplies par ce spectre (pour ne remplacer qu'une valeur venue de lui)
        let mut filled: std::collections::HashSet<String> = std::collections::HashSet::new();
        // Nom du fichier de la base d'origine : champ FILENAME du spectre, sinon premier synonyme rempli
        let mut original_filename = self.source_filename.take();
        // Identifiants de bibliothèque lus à part (voir library_spectrum_id)
        let usi = self.values.get("usi").cloned();
        let title = self.values.get("title").cloned();
        for key in self.order {
            if let (Some(mapped), Some(value)) = (keys_dict.get(&key), self.values.remove(&key)) {
                if mapped == "FILENAME" && key != "filename" {
                    if original_filename.is_none() && !is_empty_value(&value) {
                        original_filename = Some(value.trim().to_string());
                    }
                    continue;
                }
                if mapped == "FILENAME" {
                    // Champ du lecteur : nom du fichier d'entrée
                    target.insert(crate::globals_vars::INPUT_FILENAME_KEY.to_string(), value.clone());
                }
                if keep(mapped) {
                    let value = keep_unit_from_key(&key, value);
                    match target.get(mapped) {
                        None => {
                            target.insert(mapped.clone(), value);
                            filled.insert(mapped.clone());
                        }
                        // Premier champ vide ou « N/A » : remplacé par une vraie valeur d'un champ suivant
                        Some(old) if filled.contains(mapped) && is_empty_value(old) && !is_empty_value(&value) => {
                            target.insert(mapped.clone(), value);
                        }
                        Some(_) => {}
                    }
                }
            }
        }
        if let Some(original) = original_filename {
            target.insert("FILENAME".to_string(), original);
        }
        if keep("SPECTRUMID") {
            let current = target.get("SPECTRUMID").map(|s| s.trim().to_string()).unwrap_or_default();
            let numbered = !current.is_empty() && current.chars().all(|c| c.is_ascii_digit());
            if is_empty_value(&current) || numbered {
                if let Some(id) = library_spectrum_id(usi.as_deref().filter(|_| numbered), title.as_deref()) {
                    target.insert("SPECTRUMID".to_string(), id);
                }
            }
        }
    }
}

/// Identifiant du spectre dans sa bibliothèque, quand SPECTRUMID manque ou n'est qu'un numéro d'ordre :
/// - USI d'un seul spectre, quand SPECTRUMID n'est qu'un numéro (MSnLib 2025 : SPECTRUMID=1 et
///   USI "mzspec:MSV000094528:20250828_targetmol_hts_np_pos_ms2:1", le spectre dans la bibliothèque).
///   Sans SPECTRUMID (MSnLib 2024), l'USI désigne les scans bruts d'où vient le spectre (une liste pour un
///   spectre fusionné), pas le spectre de la bibliothèque : ignoré ;
/// - TITLE de Spectraverse ("SPECTRAVERSE0000000001"). Ailleurs, TITLE est un texte libre : ignoré.
fn library_spectrum_id(usi: Option<&str>, title: Option<&str>) -> Option<String> {
    if let Some(u) = usi.map(str::trim) {
        if u.starts_with("mzspec:") && !u.contains(',') && !u.contains('[') && !u.contains(char::is_whitespace) {
            return Some(u.to_string());
        }
    }
    if let Some(t) = title.map(str::trim) {
        if t.len() > "SPECTRAVERSE".len() && t.starts_with("SPECTRAVERSE") && t["SPECTRAVERSE".len()..].chars().all(|c| c.is_ascii_digit()) {
            return Some(t.to_string());
        }
    }
    None
}

/// Champs mis par le lecteur en tête de chaque spectre : un champ du même nom dans le fichier est ignoré.
const READER_KEYS: &[&str] = &["filename", "filehash"];

/// Valeur vide ou « vide » écrit (N/A, nan, unknown, 0...) : même règle que normalize_empties.
fn is_empty_value(v: &str) -> bool {
    crate::globals_vars::EMPTY_PATTERN.is_match(v.trim())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn keys() -> HashMap<String, String> {
        [("instrument_type", "INSTRUMENTTYPE"), ("source_instrument", "INSTRUMENTTYPE"),
         ("synon", "SYNON"), ("name", "NAME"), ("inconnu_hors_liste", "HORS_LISTE")]
            .iter().map(|(a, b)| (a.to_string(), b.to_string())).collect()
    }

    fn map(fields: &[(&str, &str)]) -> HashMap<String, String> {
        let mut raw = RawMetadata::default();
        for (k, v) in fields {
            raw.insert(k.to_string(), v.to_string());
        }
        let mut out = HashMap::new();
        raw.map_into(&keys(), |k| k != "HORS_LISTE", &mut out);
        out
    }

    #[test]
    fn rtinseconds_garde_son_unite() {
        let mut k = keys();
        k.insert("rtinseconds".into(), "RT".into());
        k.insert("retentiontime".into(), "RT".into());
        let mut raw = RawMetadata::default();
        raw.insert("rtinseconds".into(), "70.99".into());
        let mut out = HashMap::new();
        raw.map_into(&k, |_| true, &mut out);
        assert_eq!(out["RT"], "70.99 s");
        // Autres champs de temps : nombre seul inchangé (minutes par défaut)
        let mut raw = RawMetadata::default();
        raw.insert("retentiontime".into(), "5.2".into());
        let mut out = HashMap::new();
        raw.map_into(&k, |_| true, &mut out);
        assert_eq!(out["RT"], "5.2");
        // Valeur vide ou non numérique : inchangée
        assert_eq!(keep_unit_from_key("rtinseconds", "".into()), "");
        assert_eq!(keep_unit_from_key("rtinseconds", "NaN".into()), "NaN");
        assert_eq!(keep_unit_from_key("rtinseconds", " 12,5 ".into()), "12,5 s");
    }

    #[test]
    fn le_premier_champ_du_fichier_gagne() {
        let a = map(&[("instrument_type", "Orbitrap"), ("source_instrument", "Orbitrap ID-X")]);
        assert_eq!(a["INSTRUMENTTYPE"], "Orbitrap");
        let b = map(&[("source_instrument", "Orbitrap ID-X"), ("instrument_type", "Orbitrap")]);
        assert_eq!(b["INSTRUMENTTYPE"], "Orbitrap ID-X");
    }

    #[test]
    fn meme_champ_repete_la_derniere_valeur_gagne_comme_avant() {
        let a = map(&[("synon", "premier"), ("name", "x"), ("synon", "dernier")]);
        assert_eq!(a["SYNON"], "dernier");
        assert_eq!(a["NAME"], "x");
    }

    #[test]
    fn toutes_les_lignes_synon_gardees_dans_la_cle_interne() {
        let a = map(&[("synon", "$:07 In-Silico-Spectrum"), ("name", "x"), ("synon", "$:00in-source")]);
        assert_eq!(a["SYNON"], "$:00in-source");
        let sep = crate::globals_vars::INTERNAL_LIST_SEP;
        assert_eq!(a[crate::globals_vars::ALL_SYNONS_KEY], format!("$:07 In-Silico-Spectrum{sep}$:00in-source"));
        // Sans ligne Synon : pas de clé interne
        assert!(!map(&[("name", "x")]).contains_key(crate::globals_vars::ALL_SYNONS_KEY));
    }

    #[test]
    fn resultat_identique_a_chaque_appel() {
        // Avant, ce cas dépendait du tirage aléatoire de la HashMap : on le répète pour s'en assurer.
        for _ in 0..200 {
            let a = map(&[("instrument_type", "Orbitrap"), ("name", "x"), ("source_instrument", "Orbitrap ID-X")]);
            assert_eq!(a["INSTRUMENTTYPE"], "Orbitrap");
        }
    }

    #[test]
    fn premier_champ_vide_remplace_par_un_synonyme_rempli() {
        // GNPS JSON : INCHI_AUX ("N/A") puis InChIKey_smiles, tous deux -> INCHIKEY
        let mut k = keys();
        k.insert("inchi_aux".into(), "INCHIKEY".into());
        k.insert("inchikey_smiles".into(), "INCHIKEY".into());
        let run = |fields: &[(&str, &str)]| {
            let mut raw = RawMetadata::default();
            for (a, b) in fields { raw.insert(a.to_string(), b.to_string()); }
            let mut out = HashMap::new();
            raw.map_into(&k, |_| true, &mut out);
            out
        };
        assert_eq!(run(&[("inchi_aux", "N/A"), ("inchikey_smiles", "IYDKWWDUBYWQGF-NNAZGLEUSA-N")])["INCHIKEY"],
                   "IYDKWWDUBYWQGF-NNAZGLEUSA-N");
        assert_eq!(run(&[("inchi_aux", ""), ("inchikey_smiles", "AAA")])["INCHIKEY"], "AAA");
        // Premier champ rempli : il gagne, comme avant
        assert_eq!(run(&[("inchi_aux", "CGZK"), ("inchikey_smiles", "AUTRE")])["INCHIKEY"], "CGZK");
        // Les deux vides : le premier reste (aucune vraie valeur)
        assert_eq!(run(&[("inchi_aux", "N/A"), ("inchikey_smiles", "nan")])["INCHIKEY"], "N/A");
    }

    #[test]
    fn filename_de_la_base_d_origine_garde_fichier_d_entree_en_interne() {
        let mut k = keys();
        k.insert("filename".into(), "FILENAME".into());
        k.insert("source_file".into(), "FILENAME".into());
        k.insert("filehash".into(), "FILEHASH".into());
        let input = crate::globals_vars::INPUT_FILENAME_KEY;
        let run = |fields: &[(&str, &str)]| {
            let mut raw = RawMetadata::default();
            for (a, b) in fields { raw.insert(a.to_string(), b.to_string()); }
            let mut out = HashMap::new();
            raw.map_into(&k, |_| true, &mut out);
            out
        };
        // Champ FILENAME dans le spectre (MGF MSnLib / GNPS) : nom d'origine ; fichier d'entrée en interne
        let out = run(&[("filename", "entree.mgf"), ("filehash", "h"), ("name", "x"), ("filename", "acquisition.mzML"), ("filehash", "autre")]);
        assert_eq!(out["FILENAME"], "acquisition.mzML");
        assert_eq!(out[input], "entree.mgf");
        assert_eq!(out["FILEHASH"], "h");
        // Synonyme (GNPS JSON source_file)
        let out = run(&[("filename", "ALL_GNPS.json"), ("source_file", "130618_Ger_Jenia.mgf")]);
        assert_eq!(out["FILENAME"], "130618_Ger_Jenia.mgf");
        assert_eq!(out[input], "ALL_GNPS.json");
        // Champ du spectre vide ou « N/A » : nom du fichier d'entrée
        let out = run(&[("filename", "entree.mgf"), ("filename", "N/A"), ("source_file", "")]);
        assert_eq!(out["FILENAME"], "entree.mgf");
        // Aucun nom d'origine : nom du fichier d'entrée
        let out = run(&[("filename", "entree.msp"), ("name", "x")]);
        assert_eq!(out["FILENAME"], "entree.msp");
        assert_eq!(out[input], "entree.msp");
    }

    #[test]
    fn valeur_de_repli_ne_remplace_pas_l_originale() {
        let mut raw = RawMetadata::default();
        raw.insert("smiles".into(), "CO/N=C".into());
        raw.insert_if_absent("smiles".into(), "O=C(O)".into());
        raw.insert_if_absent("inchi".into(), "InChI=1S/X".into());
        raw.insert("formula".into(), "N/A".into());
        raw.insert_if_absent("formula".into(), "C16H16".into());
        assert_eq!(raw.get("smiles").unwrap(), "CO/N=C");
        assert_eq!(raw.get("inchi").unwrap(), "InChI=1S/X");
        assert_eq!(raw.get("formula").unwrap(), "C16H16");
    }

    #[test]
    fn energies_de_collision_numerotees_reunies() {
        let mut k = keys();
        k.insert("collision_energy_1".into(), "COLLISIONENERGY".into());
        let run = |fields: &[(&str, &str)]| {
            let mut raw = RawMetadata::default();
            for (a, b) in fields { raw.insert(a.to_string(), b.to_string()); }
            let mut out = HashMap::new();
            raw.map_into(&k, |_| true, &mut out);
            out.get("COLLISIONENERGY").cloned().unwrap_or_default()
        };
        assert_eq!(run(&[("collision_energy_1", "11.326204"), ("collision_energy_2", "nan"), ("collision_energy_3", "nan")]), "11.326204");
        assert_eq!(run(&[("collision_energy_1", "5.36"), ("collision_energy_2", "8.04"), ("collision_energy_3", "10.72")]), "5.36, 8.04, 10.72");
        assert_eq!(run(&[("collision_energy_1", "nan"), ("collision_energy_2", "nan")]), "nan");
    }

    #[test]
    fn identifiant_de_bibliotheque_usi_ou_title_spectraverse() {
        let mut k = keys();
        k.insert("spectrumid".into(), "SPECTRUMID".into());
        k.insert("title".into(), "TITLE".into());
        let run = |fields: &[(&str, &str)]| {
            let mut raw = RawMetadata::default();
            for (a, b) in fields { raw.insert(a.to_string(), b.to_string()); }
            let mut out = HashMap::new();
            raw.map_into(&k, |m| m != "TITLE", &mut out);
            out.get("SPECTRUMID").cloned().unwrap_or_default()
        };
        // MSnLib 2025 : SPECTRUMID n'est qu'un numéro, l'USI identifie le spectre
        assert_eq!(run(&[("spectrumid", "1"), ("usi", "mzspec:MSV000094528:20250828_targetmol_hts_np_pos_ms2:1")]),
                   "mzspec:MSV000094528:20250828_targetmol_hts_np_pos_ms2:1");
        // MSnLib 2024 : USI des scans bruts (liste ou un seul), sans SPECTRUMID : pas d'identifiant
        assert_eq!(run(&[("usi", "[mzspec:MSV1:a:176, mzspec:MSV1:a:177]")]), "");
        assert_eq!(run(&[("usi", "mzspec:MSV1:a_id_MSn_negative:146")]), "");
        // Spectraverse
        assert_eq!(run(&[("title", "SPECTRAVERSE0000000001")]), "SPECTRAVERSE0000000001");
        // TITLE libre : ignoré
        assert_eq!(run(&[("title", "File: run1.raw Scan 123")]), "");
        // Vrai identifiant déjà présent (GNPS) : gardé
        assert_eq!(run(&[("spectrumid", "CCMSLIB00000001547"), ("usi", "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000001547")]),
                   "CCMSLIB00000001547");
    }

    #[test]
    fn cles_hors_liste_ignorees() {
        let a = map(&[("inconnu_hors_liste", "v"), ("champ_sans_correspondance", "w")]);
        assert!(a.is_empty());
    }
}
