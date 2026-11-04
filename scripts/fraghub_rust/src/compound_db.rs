// src/compound_db.rs
//! Base interne des composés (PubChem + ClassyFire + NPClassifier), gardée en RAM sous forme compacte.
//!
//! Pour un développeur Python : avant, chaque composé était un petit dictionnaire Python-like
//! (`HashMap<String, HashMap<String, String>>`) : les 13 noms de colonnes étaient recopiés
//! 7 millions de fois, les nombres et les classes convertis en texte, avec des dizaines de petites
//! allocations par ligne (≈ 14 Go au total).
//!
//! Ici on garde directement les colonnes Polars lues dans les fichiers parquet (format Arrow :
//! le texte est stocké d'un bloc, les classes chimiques sont des catégories stockées une seule fois,
//! les masses restent des f64), et on ajoute un index INCHIKEY → numéro de ligne.
//!
//! Le résultat d'une recherche est EXACTEMENT celui de l'ancienne version :
//! - la valeur est convertie en texte au moment de la lecture, avec la même règle qu'avant
//!   (`valeur_en_texte`, copie de l'ancien code) ;
//! - si un INCHIKEY apparaît plusieurs fois, c'est la dernière occurrence lue qui gagne,
//!   comme avec l'ancien `HashMap::insert` / `extend` ;
//! - une colonne absente d'un fichier renvoie `None`, une valeur nulle renvoie `Some("")`.
//!
//! Ce module ne dépend pas de PyO3, pour pouvoir être testé seul.

use hashbrown::HashMap;
use polars::prelude::*;
use std::fs;
use std::path::Path;

/// Longueur d'un InChIKey standard (ex. "PGWQMEOSXZWQQQ-UHFFFAOYSA-N").
const STANDARD_KEY_LEN: usize = 27;

/// Conversion d'une cellule Polars en texte : même règle que l'ancien `read_parquet_to_dict_of_dicts`.
pub fn valeur_en_texte(val: AnyValue) -> String {
    match val {
        AnyValue::String(s) => s.to_string(),
        AnyValue::StringOwned(s) => s.to_string(),
        AnyValue::Null => "".to_string(),
        _ => val.to_string(),
    }
}

/// Les colonnes d'un fichier parquet, alignées sur la liste globale des colonnes (`CompoundDb::columns`).
/// `None` = colonne absente de ce fichier.
#[derive(Default)]
struct Frame {
    columns: Vec<Option<Series>>,
}

/// Emplacement d'un composé dans la base (renvoyé par `locate`, utilisé par `value`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RowLocation {
    pub frame: usize,
    pub row: usize,
}

#[derive(Default)]
pub struct CompoundDb {
    /// Noms des colonnes (hors colonne clé), dans l'ordre où on les a rencontrées.
    columns: Vec<String>,
    frames: Vec<Frame>,
    /// Numéro global de la première ligne de chaque fichier (pour retrouver le fichier d'une ligne).
    frame_starts: Vec<u32>,
    /// Index des InChIKey standards (27 caractères) : clé stockée en place, sans allocation.
    index_standard: HashMap<[u8; STANDARD_KEY_LEN], u32>,
    /// Index des autres clés (longueur inhabituelle), au cas où.
    index_other: HashMap<Box<str>, u32>,
    /// Nombre total de lignes chargées (tous fichiers confondus).
    total_rows: usize,
}

impl CompoundDb {
    /// Charge un dossier de fichiers parquet (ou un fichier unique), comme l'ancien `load_internal_databases`.
    /// Les fichiers illisibles sont ignorés silencieusement, comme avant.
    pub fn load(parquet_path: &Path, key_col: &str) -> CompoundDb {
        let mut db = CompoundDb::default();
        if parquet_path.is_dir() {
            if let Ok(entries) = fs::read_dir(parquet_path) {
                // Même ordre de lecture qu'avant (ordre de read_dir) : important quand une clé est en double.
                for entry in entries.flatten() {
                    let path = entry.path();
                    if path.extension().and_then(|s| s.to_str()) == Some("parquet") {
                        db.add_file(&path, key_col);
                    }
                }
            }
        } else {
            db.add_file(parquet_path, key_col);
        }
        db.finish();
        db
    }

    fn add_file(&mut self, path: &Path, key_col: &str) {
        let df = match fs::File::open(path) {
            Ok(mut file) => match ParquetReader::new(&mut file).finish() {
                Ok(df) => df,
                Err(_) => return,
            },
            Err(_) => return,
        };
        self.add_dataframe(df, key_col);
    }

    /// Ajoute les lignes d'un DataFrame (un fichier parquet).
    pub fn add_dataframe(&mut self, df: DataFrame, key_col: &str) {
        // Sans colonne clé, l'ancienne version n'insérait aucune ligne : on fait pareil.
        let key_series = match df.column(key_col) {
            Ok(s) => s.clone(),
            Err(_) => return,
        };
        let n_rows = df.height();
        let start = match self.total_rows.checked_add(n_rows).and_then(|end| u32::try_from(end).ok()) {
            Some(_) => self.total_rows as u32,
            None => {
                eprintln!("⚠ Base interne trop grande pour l'index (plus de 4 milliards de lignes) : fichier ignoré.");
                return;
            }
        };

        // Colonnes du fichier, rangées selon la liste globale.
        let mut frame = Frame { columns: vec![None; self.columns.len()] };
        for series in df.get_columns() {
            if series.name() == key_col {
                continue;
            }
            let id = match self.columns.iter().position(|c| c == series.name()) {
                Some(id) => id,
                None => {
                    self.columns.push(series.name().to_string());
                    frame.columns.push(None);
                    self.columns.len() - 1
                }
            };
            frame.columns[id] = Some(series.clone());
        }

        // Index : la dernière occurrence d'une clé écrase les précédentes (comme HashMap::insert avant).
        let mut insert = |key: &str, local_row: usize| {
            if key.is_empty() {
                return; // l'ancienne version ignorait les lignes sans clé
            }
            let global = start + local_row as u32;
            if key.len() == STANDARD_KEY_LEN {
                let mut k = [0u8; STANDARD_KEY_LEN];
                k.copy_from_slice(key.as_bytes());
                self.index_standard.insert(k, global);
            } else {
                self.index_other.insert(key.into(), global);
            }
        };
        match key_series.str() {
            // Cas normal : colonne texte, lecture rapide. None (nul) donne "", comme valeur_en_texte.
            Ok(ca) => {
                for (i, key) in ca.into_iter().enumerate() {
                    insert(key.unwrap_or(""), i);
                }
            }
            // Autre type : même conversion que pour les valeurs.
            Err(_) => {
                for i in 0..n_rows {
                    if let Ok(v) = key_series.get(i) {
                        insert(&valeur_en_texte(v), i);
                    }
                }
            }
        }

        self.frame_starts.push(start);
        self.frames.push(frame);
        self.total_rows += n_rows;
    }

    /// Aligne tous les fichiers sur la liste finale des colonnes et libère la mémoire en trop.
    fn finish(&mut self) {
        let n = self.columns.len();
        for frame in &mut self.frames {
            frame.columns.resize(n, None);
        }
        self.index_standard.shrink_to_fit();
        self.index_other.shrink_to_fit();
    }

    /// Nombre de composés distincts dans l'index.
    pub fn len(&self) -> usize {
        self.index_standard.len() + self.index_other.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Numéro de colonne à passer à `value` (à calculer une fois, avant une boucle).
    pub fn column_id(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|c| c == name)
    }

    /// Cherche un composé par son INCHIKEY.
    pub fn locate(&self, key: &str) -> Option<RowLocation> {
        let global = if key.len() == STANDARD_KEY_LEN {
            let mut k = [0u8; STANDARD_KEY_LEN];
            k.copy_from_slice(key.as_bytes());
            *self.index_standard.get(&k)?
        } else {
            *self.index_other.get(key)?
        };
        // Fichier contenant cette ligne : le dernier dont la première ligne est <= global.
        let frame = self.frame_starts.partition_point(|&s| s <= global) - 1;
        Some(RowLocation { frame, row: (global - self.frame_starts[frame]) as usize })
    }

    /// Valeur d'une colonne pour un composé, en texte, comme l'ancien `row_dict.get(col)` :
    /// `None` si la colonne n'existe pas dans le fichier du composé, `Some("")` si la valeur est nulle.
    pub fn value(&self, loc: RowLocation, column_id: usize) -> Option<String> {
        let series = self.frames.get(loc.frame)?.columns.get(column_id)?.as_ref()?;
        series.get(loc.row).ok().map(valeur_en_texte)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap as StdHashMap;

    /// Copie EXACTE de l'ancienne fonction de loading_db.rs (avant le passage au format compact),
    /// gardée ici comme référence : le nouveau format doit donner les mêmes réponses.
    fn ancienne_lecture(filepath: &str, key_col: &str) -> StdHashMap<String, StdHashMap<String, String>> {
        let mut map = StdHashMap::new();
        if let Ok(mut file) = fs::File::open(filepath) {
            if let Ok(df) = ParquetReader::new(&mut file).finish() {
                let headers: Vec<String> = df.get_column_names().iter().map(|s| s.to_string()).collect();
                let n_rows = df.height();
                let columns = df.get_columns();
                for i in 0..n_rows {
                    let mut row_dict = StdHashMap::new();
                    let mut key_val = String::new();
                    for (col_idx, header) in headers.iter().enumerate() {
                        if let Ok(val) = columns[col_idx].get(i) {
                            let val_clean = match val {
                                AnyValue::String(s) => s.to_string(),
                                AnyValue::StringOwned(s) => s.to_string(),
                                AnyValue::Null => "".to_string(),
                                _ => val.to_string(),
                            };
                            if header == key_col {
                                key_val = val_clean.clone();
                            }
                            row_dict.insert(header.clone(), val_clean);
                        }
                    }
                    if !key_val.is_empty() {
                        map.insert(key_val, row_dict);
                    }
                }
            }
        }
        map
    }

    /// Ancien chargement complet d'un dossier (même ordre de lecture que CompoundDb::load).
    fn ancien_chargement(dir: &Path) -> StdHashMap<String, StdHashMap<String, String>> {
        let mut all = StdHashMap::new();
        for entry in fs::read_dir(dir).unwrap().flatten() {
            let path = entry.path();
            if path.extension().and_then(|s| s.to_str()) == Some("parquet") {
                all.extend(ancienne_lecture(&path.to_string_lossy(), "INCHIKEY"));
            }
        }
        all
    }

    /// Vérifie que la nouvelle base répond exactement comme l'ancienne, pour toutes les clés et colonnes.
    pub fn compare(dir: &Path) -> usize {
        let old = ancien_chargement(dir);
        let db = CompoundDb::load(dir, "INCHIKEY");
        assert_eq!(db.len(), old.len(), "nombre de composés différent");
        let mut all_cols: Vec<String> = old.values().flat_map(|r| r.keys().cloned()).collect();
        all_cols.sort();
        all_cols.dedup();
        all_cols.push("COLONNE_INEXISTANTE".to_string());
        let mut checked = 0;
        for (key, row) in &old {
            let loc = db.locate(key).unwrap_or_else(|| panic!("clé absente : {key}"));
            for col in &all_cols {
                if col == "INCHIKEY" {
                    continue; // jamais demandée par le pipeline
                }
                let new_val = db.column_id(col).and_then(|id| db.value(loc, id));
                assert_eq!(new_val.as_deref(), row.get(col).map(|s| s.as_str()), "clé {key}, colonne {col}");
                checked += 1;
            }
        }
        assert!(db.locate("CLE-QUI-NEXISTE-PAS-XXXXXXX").is_none());
        checked
    }

    fn write_parquet(dir: &Path, name: &str, mut df: DataFrame) {
        let file = fs::File::create(dir.join(name)).unwrap();
        ParquetWriter::new(file).finish(&mut df).unwrap();
    }

    #[test]
    fn meme_resultat_que_l_ancienne_version() {
        let dir = std::env::temp_dir().join(format!("fraghub_compound_db_test_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        // Fichier 1 : clé en double dans le fichier, clé nulle, clé vide, clé non standard,
        // masse nulle et NaN, classe en catégorie.
        let classes = Series::new("NPCLASS_PATHWAY", &[Some("Alkaloids"), None, Some("Terpenoids"), Some("Alkaloids"), Some("X"), Some("Y")])
            .cast(&DataType::Categorical(None, Default::default())).unwrap();
        let df1 = DataFrame::new(vec![
            Series::new("INCHIKEY", &[Some("AAAAAAAAAAAAAA-BBBBBBBBBB-N"), Some("CCCCCCCCCCCCCC-DDDDDDDDDD-N"), Some("AAAAAAAAAAAAAA-BBBBBBBBBB-N"), None, Some(""), Some("CLE-COURTE")]),
            Series::new("NAME", &[Some("premier"), Some("autre"), Some("doublon gagnant"), Some("sans cle"), Some("cle vide"), None]),
            Series::new("EXACTMASS", &[Some(190.120509056), None, Some(f64::NAN), Some(1.0), Some(2.0), Some(1e-7)]),
            classes,
        ]).unwrap();
        // Fichier 2 : écrase une clé du fichier 1, colonne NAME absente, nouvelle colonne.
        let df2 = DataFrame::new(vec![
            Series::new("NPCLASS_PATHWAY", &[Some("Polyketides"), Some("Fatty acids")]),
            Series::new("INCHIKEY", &[Some("CCCCCCCCCCCCCC-DDDDDDDDDD-N"), Some("EEEEEEEEEEEEEE-FFFFFFFFFF-N")]),
            Series::new("SMILES", &[Some("CCO"), None]),
        ]).unwrap();
        write_parquet(&dir, "part_0000.parquet", df1);
        write_parquet(&dir, "part_0001.parquet", df2);
        fs::write(dir.join("notes.txt"), "pas un parquet").unwrap();

        let checked = compare(&dir);
        assert!(checked > 0);

        // Quelques réponses attendues, écrites à la main.
        let db = CompoundDb::load(&dir, "INCHIKEY");
        assert_eq!(db.len(), 4);
        let a = db.locate("AAAAAAAAAAAAAA-BBBBBBBBBB-N").unwrap();
        assert_eq!(db.value(a, db.column_id("NAME").unwrap()).as_deref(), Some("doublon gagnant"));
        let c = db.locate("CCCCCCCCCCCCCC-DDDDDDDDDD-N").unwrap();
        assert_eq!(db.value(c, db.column_id("NAME").unwrap()), None); // colonne absente du fichier 2
        assert_eq!(db.value(c, db.column_id("SMILES").unwrap()).as_deref(), Some("CCO"));
        assert!(db.locate("CLE-COURTE").is_some());

        let _ = fs::remove_dir_all(&dir);
    }

    /// Mémoire occupée par le processus (Linux uniquement), pour mesurer le gain.
    fn memoire_mo() -> String {
        fs::read_to_string("/proc/self/status").ok()
            .and_then(|s| s.lines().find(|l| l.starts_with("VmRSS")).map(|l| l.to_string()))
            .unwrap_or_else(|| "mesure indisponible".to_string())
    }

    /// Même contrôle que `compare`, mais fichier par fichier, pour ne jamais avoir l'ancienne base
    /// complète en mémoire (≈ 14 Go). Pour chaque fichier, dans l'ordre de lecture :
    /// - les clés dont ce fichier est la dernière occurrence doivent donner exactement les mêmes valeurs ;
    /// - les autres doivent être écrasées par un fichier lu plus tard (comme avec l'ancien `extend`).
    fn compare_par_fichier(dir: &Path) -> (usize, usize) {
        let db = CompoundDb::load(dir, "INCHIKEY");
        println!("Nouvelle base chargée : {} composés, {}", db.len(), memoire_mo());
        let mut cols = db.columns.clone();
        cols.push("COLONNE_INEXISTANTE".to_string());
        let (mut gagnants, mut valeurs) = (0usize, 0usize);
        let mut frame = 0usize;
        for entry in fs::read_dir(dir).unwrap().flatten() {
            let path = entry.path();
            if path.extension().and_then(|s| s.to_str()) != Some("parquet") {
                continue;
            }
            // Même règle que CompoundDb::add_file pour savoir si ce fichier a donné un "frame".
            let lisible = fs::File::open(&path).ok()
                .and_then(|mut f| ParquetReader::new(&mut f).finish().ok())
                .map(|df| df.column("INCHIKEY").is_ok())
                .unwrap_or(false);
            let old = ancienne_lecture(&path.to_string_lossy(), "INCHIKEY");
            if !lisible {
                assert!(old.is_empty());
                continue;
            }
            for (key, row) in &old {
                let loc = db.locate(key).unwrap_or_else(|| panic!("clé absente : {key}"));
                if loc.frame == frame {
                    gagnants += 1;
                    for col in &cols {
                        let new_val = db.column_id(col).and_then(|id| db.value(loc, id));
                        assert_eq!(new_val.as_deref(), row.get(col).map(|s| s.as_str()), "clé {key}, colonne {col}");
                        valeurs += 1;
                    }
                } else {
                    assert!(loc.frame > frame, "clé {key} : un fichier plus ancien a gagné");
                }
            }
            println!("  {:?} : {} lignes contrôlées", path.file_name().unwrap(), old.len());
            frame += 1;
        }
        // Chaque composé de la nouvelle base a été retrouvé exactement une fois comme "gagnant".
        assert_eq!(gagnants, db.len(), "nombre de composés différent");
        (gagnants, valeurs)
    }

    /// Comparaison sur les vraies données :
    /// FRAGHUB_DB_DIR=.../datas/internal_databases.parquet cargo test --release -- --ignored --nocapture
    #[test]
    #[ignore]
    fn meme_resultat_sur_les_vraies_donnees() {
        let dir = std::env::var("FRAGHUB_DB_DIR").expect("définir FRAGHUB_DB_DIR");
        let (composes, valeurs) = compare_par_fichier(Path::new(&dir));
        println!("{composes} composés et {valeurs} valeurs comparés : tous identiques");
    }
}
