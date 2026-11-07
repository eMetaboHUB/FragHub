// src/database_renormalization.rs
//! Re-normalisation de la base interne (datas/internal_databases.parquet) avec EXACTEMENT le traitement
//! RDKit appliqué aux spectres (rdkit_native::process_single) : lecture, neutralisation, InChI et InChIKey
//! standard de la molécule neutralisée, SMILES normalisé (tautomère canonique s'il garde le même InChI).
//!
//! But : que l'InChIKey d'un spectre et celle de la base soient calculées de la même façon, pour que la
//! même molécule se retrouve (étape PubChem, ontologies), et que les SMILES soient écrits de la même façon.
//!
//! Règles :
//! - colonnes réécrites : INCHI, INCHIKEY, SMILES, FORMULA, EXACTMASS, AVERAGEMASS (celles qui existent) ;
//!   NAME et ontologies sont recopiés tels quels ;
//! - la molécule de départ est choisie comme dans le pipeline : l'INCHI s'il existe, sinon le SMILES ;
//! - pas de dessalage (pas de m/z précurseur dans la base) : comme un spectre sans m/z ;
//! - une molécule illisible par RDKit garde ses valeurs d'origine ;
//! - deux lignes peuvent tomber sur la même InChIKey : elles sont comptées ; à la lecture, FragHub garde la
//!   dernière (même règle que pour les doublons déjà présents) ;
//! - un fichier sans colonne de structure est recopié octet pour octet ;
//! - les fichiers sont écrits dans un NOUVEAU dossier : la base d'origine n'est jamais modifiée ;
//! - chaque fichier est d'abord écrit sous un nom temporaire puis renommé : un fichier présent dans le
//!   dossier de sortie est complet. Un traitement interrompu (Ctrl+C) reprend donc au fichier suivant.
use polars::prelude::*;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::fs;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use crate::globals_vars::INDIGO_SMILES_CORRECTION_PATTERN;

/// Compteurs d'une re-normalisation.
#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct RenormStats {
    pub rows: usize,
    /// Lignes avec une structure (INCHI ou SMILES)
    pub with_structure: usize,
    /// InChIKey différente après normalisation
    pub key_changed: usize,
    /// SMILES différent après normalisation
    pub smiles_changed: usize,
    /// Molécule illisible par RDKit : valeurs d'origine gardées
    pub failed: usize,
    /// InChIKey présentes sur plusieurs lignes après normalisation (calculé sur toute la base)
    pub duplicate_keys: usize,
}

impl RenormStats {
    fn add(&mut self, o: &RenormStats) {
        self.rows += o.rows; self.with_structure += o.with_structure; self.key_changed += o.key_changed;
        self.smiles_changed += o.smiles_changed; self.failed += o.failed;
    }
}

/// Colonne texte d'un tableau, ou None si elle manque ou n'est pas du texte (ex. colonne entièrement vide).
fn text_column<'a>(df: &'a DataFrame, name: &str) -> Option<&'a StringChunked> {
    df.column(name).ok().and_then(|s| s.str().ok())
}

fn usable(v: Option<&str>) -> Option<&str> {
    v.filter(|s| !s.is_empty() && *s != "nan")
}

const TEXT_COLS: [&str; 4] = ["INCHI", "INCHIKEY", "SMILES", "FORMULA"];
const MASS_COLS: [&str; 2] = ["EXACTMASS", "AVERAGEMASS"];

/// Re-normalise un tableau. Renvoie None si le tableau n'a pas de colonne de structure.
/// `tick(molécules_faites, molécules_à_faire)` est appelé toutes les 0,5 s ; s'il renvoie false, on s'arrête (erreur).
pub fn renormalize_frame(df: &DataFrame, tick: &dyn Fn(usize, usize) -> bool) -> PolarsResult<Option<(DataFrame, RenormStats)>> {
    let smiles = text_column(df, "SMILES");
    let inchis = text_column(df, "INCHI");
    if smiles.is_none() && inchis.is_none() { return Ok(None); }
    let n = df.height();
    let keys = text_column(df, "INCHIKEY");

    // Molécule de départ de chaque ligne : même choix que rdkit_bridge.rs (InChI d'abord, sinon SMILES
    // sans la partie « | ... ») ; None si la ligne n'a pas de structure.
    let targets: Vec<Option<String>> = (0..n).map(|i| {
        let inchi = inchis.and_then(|c| c.get(i));
        let smi = smiles.and_then(|c| c.get(i));
        match (usable(inchi), usable(smi)) {
            (Some(x), _) => Some(x.to_string()),
            (None, Some(s)) => Some(if s.contains("InChI=") { s.to_string() } else { INDIGO_SMILES_CORRECTION_PATTERN.replace_all(s, "").to_string() }),
            (None, None) => None,
        }
    }).collect();

    // Calcul RDKit avec une file d'attente partagée, comme dans le pipeline (rdkit_native::compute_cache) :
    // chaque thread prend la molécule suivante dès qu'il a fini la sienne, les plus longues d'abord.
    // (Avant : paquets de 5 000 molécules ; chaque paquet attendait sa molécule la plus lente pendant que
    // les autres threads ne faisaient rien.)
    let mut order: Vec<usize> = (0..n).filter(|&i| targets[i].is_some()).collect();
    order.sort_by(|&a, &b| targets[b].as_ref().unwrap().len().cmp(&targets[a].as_ref().unwrap().len()));
    let todo = order.len();
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let cancel = AtomicBool::new(false);
    let mut stopped = false;
    let mut computed: Vec<(usize, Vec<(&'static str, String)>)> = Vec::new();
    std::thread::scope(|scope| {
        let worker = scope.spawn(|| {
            let workers = rayon::current_num_threads().max(1);
            (0..workers).into_par_iter().flat_map_iter(|_| {
                let mut mine = Vec::new();
                while !cancel.load(Ordering::Relaxed) {
                    let k = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&i) = order.get(k) else { break };
                    mine.push((i, crate::rdkit_native::process_single(targets[i].as_ref().unwrap(), "", "")));
                    done.fetch_add(1, Ordering::Relaxed);
                }
                mine
            }).collect::<Vec<_>>()
        });
        // Ce thread-ci donne la progression (et peut demander l'arrêt) pendant que les autres calculent
        while !worker.is_finished() {
            std::thread::sleep(Duration::from_millis(500));
            if !stopped && !tick(done.load(Ordering::Relaxed), todo) {
                stopped = true;
                cancel.store(true, Ordering::Relaxed);
            }
        }
        match worker.join() {
            Ok(r) => computed = r,
            Err(panic) => std::panic::resume_unwind(panic),
        }
    });
    if stopped || !tick(todo, todo) {
        return Err(PolarsError::ComputeError("interrompu".into()));
    }
    // Sans InChIKey, la molécule n'a pas été lue correctement : on ne touche pas à la ligne
    let mut results: Vec<Option<Vec<(&'static str, String)>>> = vec![None; n];
    for (i, out) in computed {
        results[i] = Some(if out.iter().any(|(k, _)| *k == "INCHIKEY") { out } else { Vec::new() });
    }

    let mut stats = RenormStats { rows: n, ..Default::default() };
    for (i, r) in results.iter().enumerate() {
        match r {
            None => {}
            Some(v) if v.is_empty() => { stats.with_structure += 1; stats.failed += 1; }
            Some(v) => {
                stats.with_structure += 1;
                let g = |k: &str| v.iter().find(|x| x.0 == k).map(|x| x.1.as_str());
                if g("INCHIKEY") != keys.and_then(|c| c.get(i)) { stats.key_changed += 1; }
                if g("SMILES") != smiles.and_then(|c| c.get(i)) { stats.smiles_changed += 1; }
            }
        }
    }

    let new_value = |i: usize, col: &str| -> Option<&str> {
        results[i].as_ref().and_then(|v| v.iter().find(|x| x.0 == col).map(|x| x.1.as_str()))
    };
    let mut out = df.clone();
    for col in TEXT_COLS {
        let Some(old) = text_column(df, col) else { continue };
        let values: Vec<Option<String>> = (0..n).map(|i| match &results[i] {
            Some(v) if !v.is_empty() => new_value(i, col).map(|s| s.to_string()).or_else(|| old.get(i).map(|s| s.to_string())),
            _ => old.get(i).map(|s| s.to_string()),
        }).collect();
        out.replace(col, Series::new(col, values))?;
    }
    for col in MASS_COLS {
        let Ok(old) = df.column(col).and_then(|s| s.f64().map(|c| c.clone())) else { continue };
        let values: Vec<Option<f64>> = (0..n).map(|i| match &results[i] {
            Some(v) if !v.is_empty() => new_value(i, col).and_then(|s| s.parse::<f64>().ok()).or(old.get(i)),
            _ => old.get(i),
        }).collect();
        out.replace(col, Series::new(col, values))?;
    }
    Ok(Some((out, stats)))
}

/// Re-normalise tous les fichiers *.parquet de `input_dir` dans `output_dir` (mêmes noms de fichiers).
/// `progress(fichier, i, nb_fichiers, lignes_faites, lignes_totales)` renvoie false pour arrêter.
/// `only` : noms des fichiers à traiter (None = tous) ; sert à répartir la base entre plusieurs processus.
pub fn renormalize_database(
    input_dir: &Path,
    output_dir: &Path,
    only: Option<&[String]>,
    log: &dyn Fn(&str),
    progress: &dyn Fn(&str, usize, usize, usize, usize) -> bool,
) -> Result<RenormStats, String> {
    if input_dir == output_dir { return Err("le dossier de sortie doit être différent du dossier d'entrée".into()); }
    fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;
    let mut parts: Vec<_> = fs::read_dir(input_dir).map_err(|e| e.to_string())?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().map_or(false, |x| x == "parquet"))
        .filter(|p| only.map_or(true, |o| o.iter().any(|n| Some(n.as_str()) == p.file_name().and_then(|f| f.to_str()))))
        .collect();
    parts.sort();
    if parts.is_empty() { return Err(format!("aucun fichier .parquet dans {}", input_dir.display())); }

    let mut total = RenormStats::default();
    let t0 = Instant::now();
    for (i, path) in parts.iter().enumerate() {
        let name = path.file_name().unwrap();
        let dest = output_dir.join(name);
        if dest.exists() {
            // Fichier déjà traité lors d'un lancement précédent (écrit en entier, puis renommé)
            log(&format!("[{}/{}] {} : déjà traité, ignoré", i + 1, parts.len(), name.to_string_lossy()));
            continue;
        }
        let tmp = output_dir.join(format!("{}.en_cours", name.to_string_lossy()));
        let t = Instant::now();
        let mut file = fs::File::open(path).map_err(|e| e.to_string())?;
        let df = ParquetReader::new(&mut file).finish().map_err(|e| format!("{}: {e}", path.display()))?;
        let fname = name.to_string_lossy().to_string();
        let tick = |done: usize, total_rows: usize| progress(&fname, i + 1, parts.len(), done, total_rows);
        match renormalize_frame(&df, &tick).map_err(|e| e.to_string())? {
            None => {
                fs::copy(path, &tmp).map_err(|e| e.to_string())?;
                fs::rename(&tmp, &dest).map_err(|e| e.to_string())?;
                log(&format!("[{}/{}] {} : pas de structure, recopié tel quel ({} lignes)", i + 1, parts.len(), name.to_string_lossy(), df.height()));
                total.rows += df.height();
            }
            Some((mut out, st)) => {
                let f = fs::File::create(&tmp).map_err(|e| e.to_string())?;
                ParquetWriter::new(f).with_compression(ParquetCompression::Zstd(None)).finish(&mut out).map_err(|e| e.to_string())?;
                fs::rename(&tmp, &dest).map_err(|e| e.to_string())?;
                log(&format!("[{}/{}] {} : {} lignes, {} avec structure, {} InChIKey modifiées, {} SMILES modifiés, {} illisibles (gardées) - {:.0} s",
                    i + 1, parts.len(), name.to_string_lossy(), st.rows, st.with_structure, st.key_changed, st.smiles_changed, st.failed, t.elapsed().as_secs_f64()));
                total.add(&st);
            }
        }
    }
    // InChIKey en double après normalisation (sur toute la base)
    let mut seen = std::collections::HashMap::<String, usize>::new();
    for path in &parts {
        let dest = output_dir.join(path.file_name().unwrap());
        if let Ok(mut f) = fs::File::open(&dest) {
            if let Ok(df) = ParquetReader::new(&mut f).finish() {
                if let Some(k) = text_column(&df, "INCHIKEY") {
                    for key in k.into_iter().flatten() { *seen.entry(key.to_string()).or_insert(0) += 1; }
                }
            }
        }
    }
    total.duplicate_keys = seen.values().filter(|&&c| c > 1).count();
    log(&format!("TOTAL : {} lignes, {} avec structure, {} InChIKey modifiées, {} SMILES modifiés, {} illisibles (gardées), {} InChIKey présentes plusieurs fois - {:.0} min",
        total.rows, total.with_structure, total.key_changed, total.smiles_changed, total.failed, total.duplicate_keys, t0.elapsed().as_secs_f64() / 60.0));
    Ok(total)
}

/// Appel depuis Python : `fraghub_rust.renormalize_internal_database(dossier_entree, dossier_sortie)`.
/// `only` : liste de noms de fichiers à traiter (par défaut, tous).
/// `prefix` : texte ajouté devant chaque ligne ; s'il n'est pas vide (plusieurs processus dans la même
/// console), la progression est écrite sur une nouvelle ligne toutes les minutes au lieu d'une ligne
/// qui se met à jour.
#[pyfunction]
#[pyo3(signature = (input_dir, output_dir, only=None, prefix=String::new()))]
pub fn renormalize_internal_database(py: Python, input_dir: &str, output_dir: &str, only: Option<Vec<String>>, prefix: String) -> PyResult<(usize, usize, usize, usize, usize, usize)> {
    use std::io::Write;
    let t0 = Instant::now();
    let interrupted: std::sync::Mutex<Option<PyErr>> = std::sync::Mutex::new(None);
    // Début du fichier en cours (pour la vitesse et le temps restant) et dernière ligne écrite
    let file_start: std::sync::Mutex<(String, Instant)> = std::sync::Mutex::new((String::new(), Instant::now()));
    let last_line = std::sync::Mutex::new(Instant::now());
    let one_line = prefix.is_empty();
    let r = py.allow_threads(|| renormalize_database(
        Path::new(input_dir),
        Path::new(output_dir),
        only.as_deref(),
        &|m| println!("\r{prefix}{m}{}", " ".repeat(30)),
        // Progression, et Ctrl+C vérifié toutes les 0,5 s (Python reprend la main un instant)
        &|name, i, nb, done, total| {
            let mut fs_ = file_start.lock().unwrap();
            if fs_.0 != name { *fs_ = (name.to_string(), Instant::now()); }
            let minutes = fs_.1.elapsed().as_secs_f64() / 60.0;
            let rate = if minutes > 0.0 { done as f64 / minutes } else { 0.0 };
            let left = if rate > 0.0 { (total - done) as f64 / rate } else { f64::NAN };
            let text = format!("{prefix}[{i}/{nb}] {name} : {done}/{total} molécules, {rate:.0} molécules/min, reste ~{left:.0} min pour ce fichier ({:.0} min écoulées au total)",
                t0.elapsed().as_secs_f64() / 60.0);
            if one_line {
                print!("\r{text}   ");
            } else {
                let mut last = last_line.lock().unwrap();
                if last.elapsed() >= Duration::from_secs(60) { println!("{text}"); *last = Instant::now(); }
            }
            let _ = std::io::stdout().flush();
            match Python::with_gil(|py| py.check_signals()) {
                Ok(()) => true,
                Err(e) => { *interrupted.lock().unwrap() = Some(e); false }
            }
        },
    ));
    if let Some(e) = interrupted.into_inner().unwrap() {
        println!("\nArrêté. Les fichiers déjà terminés sont gardés : relancer le script reprend au fichier suivant.");
        return Err(e);
    }
    match r {
        Ok(s) => Ok((s.rows, s.with_structure, s.key_changed, s.smiles_changed, s.failed, s.duplicate_keys)),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compound_db::CompoundDb;

    fn frame() -> DataFrame {
        let classes = Series::new("NPCLASS_PATHWAY", &[Some("Alkaloids"), None, Some("Terpenoids"), Some("X")])
            .cast(&DataType::Categorical(None, Default::default())).unwrap();
        DataFrame::new(vec![
            // 1 : forme iminol -> SMILES amide, même clé ; 2 : éthanol inchangé ;
            // 3 : carboxylate chargé -> neutralisé, clé et formule changent ; 4 : illisible ; 5 : sans structure
            Series::new("INCHIKEY", &[Some("OKJIRPAQVSHGFK-UHFFFAOYSA-N"), Some("LFQSCWFLJHTTHZ-UHFFFAOYSA-N"), Some("QTBSBXVTEAMEQO-UHFFFAOYSA-M"), Some("XXXXXXXXXXXXXX-UHFFFAOYSA-N"), None]),
            Series::new("INCHI", &[None::<&str>, Some("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"), None, None, None]),
            Series::new("SMILES", &[Some("CC(O)=NCC(=O)O"), Some("CCO"), Some("CC(=O)[O-]"), Some("pas un smiles"), None]),
            Series::new("FORMULA", &[Some("C4H7NO3"), Some("C2H6O"), Some("C2H3O2-"), Some("?"), None]),
            Series::new("NAME", &[Some("acetylglycine"), Some("ethanol"), Some("acetate"), Some("illisible"), Some("vide")]),
            Series::new("EXACTMASS", &[Some(117.0), Some(46.0), Some(59.0), Some(f64::NAN), Some(1.0)]),
            classes.extend_constant(AnyValue::Null, 1).unwrap(),
        ]).unwrap()
    }

    #[test]
    fn colonnes_de_structure_recalculees_le_reste_recopie() {
        let df = frame();
        let (out, st) = renormalize_frame(&df, &|_, _| true).unwrap().unwrap();
        assert_eq!((st.rows, st.with_structure, st.failed), (5, 4, 1));
        let col = |c: &str| out.column(c).unwrap().str().unwrap().clone();
        let (k, s, f) = (col("INCHIKEY"), col("SMILES"), col("FORMULA"));
        // Même traitement que les spectres
        let same = |m: &str, c: &str| crate::rdkit_native::process_single(m, "", "").into_iter().find(|x| x.0 == c).unwrap().1;
        assert_eq!(s.get(0), Some("CC(=O)NCC(=O)O"));
        assert_eq!(k.get(0), Some(same("CC(O)=NCC(=O)O", "INCHIKEY").as_str()));
        assert_eq!(s.get(1), Some("CCO"));
        assert_eq!(k.get(2), Some(same("CC(=O)[O-]", "INCHIKEY").as_str())); // acide acétique neutralisé
        assert_eq!(f.get(2), Some("C2H4O2"));
        // illisible et sans structure : inchangés
        assert_eq!((k.get(3), s.get(3), f.get(3)), (Some("XXXXXXXXXXXXXX-UHFFFAOYSA-N"), Some("pas un smiles"), Some("?")));
        assert_eq!((k.get(4), s.get(4)), (None, None));
        let m = out.column("EXACTMASS").unwrap().f64().unwrap().clone();
        assert!((m.get(2).unwrap() - 60.021129).abs() < 1e-5);
        assert_eq!(m.get(4), Some(1.0));
        for c in ["NAME", "NPCLASS_PATHWAY"] {
            assert!(out.column(c).unwrap().equals_missing(df.column(c).unwrap()), "colonne {c} modifiée");
        }
    }

    #[test]
    fn base_relue_par_fraghub() {
        let dir = std::env::temp_dir().join(format!("fraghub_renorm_{}", std::process::id()));
        let (src, dst) = (dir.join("in"), dir.join("out"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&src).unwrap();
        let mut df = frame();
        ParquetWriter::new(fs::File::create(src.join("part_0000.parquet")).unwrap()).finish(&mut df).unwrap();
        let mut sans = DataFrame::new(vec![Series::new("INCHIKEY", &[Some("ZZZZZZZZZZZZZZ-UHFFFAOYSA-N")]), Series::new("EXACTMASS", &[Some(3.0)])]).unwrap();
        ParquetWriter::new(fs::File::create(src.join("part_0001.parquet")).unwrap()).finish(&mut sans).unwrap();
        let st = renormalize_database(&src, &dst, None, &|_| {}, &|_, _, _, _, _| true).unwrap();
        assert_eq!(st.duplicate_keys, 0);
        assert_eq!(fs::read(src.join("part_0001.parquet")).unwrap(), fs::read(dst.join("part_0001.parquet")).unwrap());
        let db = CompoundDb::load(&dst, "INCHIKEY");
        let key = crate::rdkit_native::process_single("CC(=O)[O-]", "", "").into_iter().find(|x| x.0 == "INCHIKEY").unwrap().1;
        let loc = db.locate(&key).expect("clé neutralisée trouvée");
        assert_eq!(db.value(loc, db.column_id("NAME").unwrap()).as_deref(), Some("acetate"));
        assert!(db.locate("ZZZZZZZZZZZZZZ-UHFFFAOYSA-N").is_some());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn interruption_puis_reprise() {
        let dir = std::env::temp_dir().join(format!("fraghub_renorm_stop_{}", std::process::id()));
        let (src, dst) = (dir.join("in"), dir.join("out"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&src).unwrap();
        for name in ["part_0000.parquet", "part_0001.parquet"] {
            let mut df = frame();
            ParquetWriter::new(fs::File::create(src.join(name)).unwrap()).finish(&mut df).unwrap();
        }
        // Arrêt pendant le 2e fichier : le 1er est complet, le 2e n'existe pas (seulement le temporaire)
        let r = renormalize_database(&src, &dst, None, &|_| {}, &|_, i, _, _, _| i < 2);
        assert!(r.is_err());
        assert!(dst.join("part_0000.parquet").exists());
        assert!(!dst.join("part_0001.parquet").exists());
        // Reprise : le 1er est ignoré, le 2e est fait ; résultat identique à un traitement d'un coup
        let st = renormalize_database(&src, &dst, None, &|_| {}, &|_, _, _, _, _| true).unwrap();
        assert_eq!(st.rows, 5);
        assert_eq!(fs::read(dst.join("part_0000.parquet")).unwrap(), fs::read(dst.join("part_0001.parquet")).unwrap());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn sous_ensemble_de_fichiers() {
        let dir = std::env::temp_dir().join(format!("fraghub_renorm_only_{}", std::process::id()));
        let (src, dst) = (dir.join("in"), dir.join("out"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&src).unwrap();
        for name in ["part_0000.parquet", "part_0001.parquet"] {
            let mut df = frame();
            ParquetWriter::new(fs::File::create(src.join(name)).unwrap()).finish(&mut df).unwrap();
        }
        let only = vec!["part_0001.parquet".to_string()];
        renormalize_database(&src, &dst, Some(&only), &|_| {}, &|_, _, _, _, _| true).unwrap();
        assert!(!dst.join("part_0000.parquet").exists());
        assert!(dst.join("part_0001.parquet").exists());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn refuse_d_ecrire_dans_le_dossier_d_entree() {
        let d = std::env::temp_dir();
        assert!(renormalize_database(&d, &d, None, &|_| {}, &|_, _, _, _, _| true).is_err());
    }
}

/// Essai sur une vraie base : FRAGHUB_RENORM_IN=... FRAGHUB_RENORM_OUT=... cargo test --release -- --ignored vraie_base --nocapture
#[cfg(test)]
#[test]
#[ignore]
fn vraie_base() {
    let (i, o) = (std::env::var("FRAGHUB_RENORM_IN").unwrap(), std::env::var("FRAGHUB_RENORM_OUT").unwrap());
    renormalize_database(Path::new(&i), Path::new(&o), None, &|m| println!("{m}"), &|_, _, _, _, _| true).unwrap();
}

