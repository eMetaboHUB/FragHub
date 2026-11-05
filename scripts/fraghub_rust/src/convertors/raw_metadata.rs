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
/// - deux champs DIFFÉRENTS qui donnent la même clé FragHub : le premier du fichier gagne ;
/// - le MÊME champ présent plusieurs fois (ex. plusieurs lignes "Synon:") : la dernière valeur
///   gagne, comme avant. Pour "Synon", toutes les valeurs sont en plus gardées dans la clé interne
///   ALL_SYNONS (jamais écrite), lue par la détection in silico.
#[derive(Default)]
pub struct RawMetadata {
    values: HashMap<String, String>,
    order: Vec<String>,
    synons: Vec<String>,
}

impl RawMetadata {
    pub fn with_capacity(n: usize) -> Self {
        RawMetadata { values: HashMap::with_capacity(n), order: Vec::with_capacity(n), synons: Vec::new() }
    }

    /// Ajoute (ou remplace) un champ. Un champ déjà présent garde sa place d'origine.
    pub fn insert(&mut self, key: String, value: String) {
        if key == "synon" {
            self.synons.push(value.clone());
        }
        if !self.values.contains_key(&key) {
            self.order.push(key.clone());
        }
        self.values.insert(key, value);
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
        for key in self.order {
            if let (Some(mapped), Some(value)) = (keys_dict.get(&key), self.values.remove(&key)) {
                if keep(mapped) {
                    target.entry(mapped.clone()).or_insert_with(|| keep_unit_from_key(&key, value));
                }
            }
        }
    }
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
    fn cles_hors_liste_ignorees() {
        let a = map(&[("inconnu_hors_liste", "v"), ("champ_sans_correspondance", "w")]);
        assert!(a.is_empty());
    }
}
