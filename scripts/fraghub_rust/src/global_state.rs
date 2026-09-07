use once_cell::sync::Lazy;
use std::collections::HashMap;
use std::sync::RwLock;

/// État global de l'application stocké en RAM.
///
/// Pour un développeur Python : En Python, on définirait simplement des variables globales.
/// En Rust, les variables globales mutables sont interdites par défaut (car très dangereuses en multithreading).
/// On utilise donc `Lazy` (initialisation au premier appel) et `RwLock` 
/// (Read-Write Lock) pour s'assurer que plusieurs cœurs peuvent lire la RAM
/// simultanément en toute sécurité.
pub struct GlobalState {
    pub pubchem_datas: HashMap<String, HashMap<String, String>>,
    pub ontologies_datas: HashMap<String, HashMap<String, String>>,
    pub adduct_dict_pos: HashMap<String, String>,
    pub adduct_massdiff_dict_pos: HashMap<String, f64>,
    pub adduct_dict_neg: HashMap<String, String>,
    pub adduct_massdiff_dict_neg: HashMap<String, f64>,
    pub keys_dict: HashMap<String, String>,
    pub keys_list: Vec<String>,
    pub instrument_tree: serde_json::Value,
}

impl GlobalState {
    pub fn new() -> Self {
        GlobalState {
            pubchem_datas: HashMap::new(),
            ontologies_datas: HashMap::new(),
            adduct_dict_pos: HashMap::new(),
            adduct_massdiff_dict_pos: HashMap::new(),
            adduct_dict_neg: HashMap::new(),
            adduct_massdiff_dict_neg: HashMap::new(),
            keys_dict: HashMap::new(),
            keys_list: vec![
                "FILENAME".to_string(), "FILEHASH".to_string(), "PREDICTED".to_string(), "SPLASH".to_string(),
                "SPECTRUMID".to_string(), "RESOLUTION".to_string(), "SYNON".to_string(), "IONIZATION".to_string(),
                "MSLEVEL".to_string(), "FRAGMENTATIONMETHODE".to_string(), "NAME".to_string(), "PRECURSORMZ".to_string(),
                "EXACTMASS".to_string(), "AVERAGEMASS".to_string(), "PRECURSORTYPE".to_string(), "INSTRUMENTTYPE".to_string(),
                "INSTRUMENT".to_string(), "SMILES".to_string(), "INCHI".to_string(), "INCHIKEY".to_string(),
                "COLLISIONENERGY".to_string(), "FORMULA".to_string(), "RT".to_string(), "IONMODE".to_string(),
                "COMMENT".to_string(), "ENTROPY".to_string(), "CLASSYFIRE_SUPERCLASS".to_string(), "CLASSYFIRE_CLASS".to_string(),
                "CLASSYFIRE_SUBCLASS".to_string(), "NPCLASS_PATHWAY".to_string(), "NPCLASS_SUPERCLASS".to_string(),
                "NPCLASS_CLASS".to_string(), "NUM PEAKS".to_string(), "PEAKS_LIST".to_string()
            ],
            instrument_tree: serde_json::Value::Null,
        }
    }
}

pub static STATE: Lazy<RwLock<GlobalState>> = Lazy::new(|| RwLock::new(GlobalState::new()));
