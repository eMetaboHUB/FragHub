// src/convertors/raw_metadata.rs
use std::collections::HashMap;

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
///   gagne, comme avant.
#[derive(Default)]
pub struct RawMetadata {
    values: HashMap<String, String>,
    order: Vec<String>,
}

impl RawMetadata {
    pub fn with_capacity(n: usize) -> Self {
        RawMetadata { values: HashMap::with_capacity(n), order: Vec::with_capacity(n) }
    }

    /// Ajoute (ou remplace) un champ. Un champ déjà présent garde sa place d'origine.
    pub fn insert(&mut self, key: String, value: String) {
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
        for key in self.order {
            if let (Some(mapped), Some(value)) = (keys_dict.get(&key), self.values.remove(&key)) {
                if keep(mapped) {
                    target.entry(mapped.clone()).or_insert(value);
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
