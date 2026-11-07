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
//!   dossier de sortie est complet ;
//! - chaque molécule calculée est ajoutée à un fichier de reprise (<fichier>.reprise) : après un arrêt
//!   (Ctrl+C, panne), le lancement suivant reprend à la molécule près ;
//! - les molécules qui prennent plus d'une minute sont notées dans _molecules_lentes.tsv (durée en s,
//!   fichier, ligne, molécule).
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

/// Résultat RDKit d'une ligne : None = pas de structure ; Some(vide) = molécule illisible (ligne gardée telle
/// quelle) ; Some(valeurs) = valeurs recalculées.
type RowResult = Option<Vec<(&'static str, String)>>;

/// Molécule de départ de chaque ligne : même choix que rdkit_bridge.rs (InChI d'abord, sinon SMILES sans la
/// partie « | ... »). None (pour le tableau entier) s'il n'y a ni colonne INCHI ni colonne SMILES.
fn targets_of(df: &DataFrame) -> Option<Vec<Option<String>>> {
    let smiles = text_column(df, "SMILES");
    let inchis = text_column(df, "INCHI");
    if smiles.is_none() && inchis.is_none() { return None; }
    Some((0..df.height()).map(|i| {
        let inchi = inchis.and_then(|c| c.get(i));
        let smi = smiles.and_then(|c| c.get(i));
        match (usable(inchi), usable(smi)) {
            (Some(x), _) => Some(x.to_string()),
            (None, Some(s)) => Some(if s.contains("InChI=") { s.to_string() } else { INDIGO_SMILES_CORRECTION_PATTERN.replace_all(s, "").to_string() }),
            (None, None) => None,
        }
    }).collect())
}

/// Calcul RDKit d'une molécule (sans InChIKey, la molécule n'a pas été lue correctement : résultat vide).
fn compute_one(target: &str) -> Vec<(&'static str, String)> {
    let out = crate::rdkit_native::process_single(target, "", "");
    if out.iter().any(|(k, _)| *k == "INCHIKEY") { out } else { Vec::new() }
}

/// Applique les résultats au tableau : colonnes réécrites et compteurs.
fn apply_results(df: &DataFrame, results: &[RowResult]) -> PolarsResult<(DataFrame, RenormStats)> {
    let n = df.height();
    let keys = text_column(df, "INCHIKEY");
    let smiles = text_column(df, "SMILES");
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
    Ok((out, stats))
}

/// Re-normalise un tableau en une fois (utilisé par les tests). None si le tableau n'a pas de structure.
pub fn renormalize_frame(df: &DataFrame) -> PolarsResult<Option<(DataFrame, RenormStats)>> {
    let Some(targets) = targets_of(df) else { return Ok(None) };
    let results: Vec<RowResult> = targets.par_iter().map(|t| t.as_deref().map(compute_one)).collect();
    apply_results(df, &results).map(Some)
}

// ------------------------------------------------------------------------------------------------
// Fichier de reprise : chaque molécule calculée y est ajoutée (une ligne JSON) ; un arrêt (Ctrl+C,
// panne, redémarrage) ne perd donc rien : au lancement suivant, les molécules déjà présentes ne sont pas
// recalculées. Le fichier est supprimé quand le fichier .parquet correspondant est écrit.
// ------------------------------------------------------------------------------------------------

const OUTPUT_KEYS: [&str; 6] = ["INCHI", "INCHIKEY", "SMILES", "FORMULA", "EXACTMASS", "AVERAGEMASS"];

fn static_key(k: &str) -> Option<&'static str> {
    OUTPUT_KEYS.iter().copied().find(|x| *x == k)
}

fn checkpoint_line(row: usize, out: &[(&'static str, String)]) -> String {
    let pairs: Vec<serde_json::Value> = out.iter().map(|(k, v)| serde_json::json!([k, v])).collect();
    format!("{}\n", serde_json::json!({"i": row, "r": pairs}))
}

/// Lit un fichier de reprise ; une dernière ligne coupée (arrêt brutal) est ignorée.
fn read_checkpoint(path: &Path, results: &mut [RowResult]) -> usize {
    let Ok(text) = fs::read_to_string(path) else { return 0 };
    let mut n = 0;
    for line in text.lines() {
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else { continue };
        let Some(i) = v["i"].as_u64().map(|x| x as usize) else { continue };
        let Some(arr) = v["r"].as_array() else { continue };
        if i >= results.len() { continue; }
        let mut out = Vec::new();
        for p in arr {
            if let (Some(k), Some(val)) = (p[0].as_str().and_then(static_key), p[1].as_str()) {
                out.push((k, val.to_string()));
            }
        }
        if results[i].is_none() { n += 1; }
        results[i] = Some(out);
    }
    n
}

/// Avancement envoyé à l'affichage.
pub struct Progress {
    pub done: usize,
    pub total: usize,
    pub files_done: usize,
    pub files_total: usize,
    /// Molécules en cours depuis plus d'une minute : (secondes, fichier, molécule)
    pub slow: Vec<(f64, String, String)>,
}

/// Seuil au-delà duquel une molécule est signalée comme lente (affichage et _molecules_lentes.tsv).
const SLOW_SECS: f64 = 60.0;

/// Au-delà de ce nombre d'atomes lourds (tous les atomes sauf H), une molécule est gardée pour la fin :
/// le tautomère des très grosses molécules (porphyrines en chaîne, nanographènes, foldamères de plus de
/// 200 atomes...) peut prendre des heures, et douze d'entre elles occupaient tous les cœurs pendant que le
/// reste de la base attendait. Réglable par la variable d'environnement FRAGHUB_RENORM_BIG_ATOMS.
/// Ne change que l'ORDRE du calcul : le résultat de chaque molécule est le même.
const BIG_ATOMS_DEFAULT: usize = 100;

fn big_atoms_threshold() -> usize {
    std::env::var("FRAGHUB_RENORM_BIG_ATOMS").ok().and_then(|v| v.trim().parse().ok()).unwrap_or(BIG_ATOMS_DEFAULT)
}

/// Nombre d'atomes lourds d'une molécule, lu dans la formule de l'InChI ("InChI=1S/C216H184N16O8/..." -> 240 ;
/// "C10H12.2ClH" -> 10 + 2). Pour un SMILES : nombre de symboles d'atomes autres que H (estimation suffisante
/// pour trier).
fn heavy_atoms(target: &str) -> usize {
    if let Some(rest) = target.strip_prefix("InChI=") {
        let formula = rest.split('/').nth(1).unwrap_or("");
        return formula.split('.').map(|comp| {
            let digits: String = comp.chars().take_while(|c| c.is_ascii_digit()).collect();
            let mult: usize = digits.parse().unwrap_or(1);
            let body = &comp[digits.len()..];
            let mut n = 0usize;
            let chars: Vec<char> = body.chars().collect();
            let mut i = 0;
            while i < chars.len() {
                if chars[i].is_ascii_uppercase() {
                    let mut sym = chars[i].to_string();
                    i += 1;
                    while i < chars.len() && chars[i].is_ascii_lowercase() { sym.push(chars[i]); i += 1; }
                    let mut num = String::new();
                    while i < chars.len() && chars[i].is_ascii_digit() { num.push(chars[i]); i += 1; }
                    if sym != "H" { n += num.parse::<usize>().unwrap_or(1); }
                } else {
                    i += 1;
                }
            }
            n * mult
        }).sum();
    }
    // SMILES : majuscules (C, N, O, Cl, Br...) et atomes aromatiques (c, n, o, s, p, b), sauf H
    let b = target.as_bytes();
    let mut n = 0;
    for (i, &c) in b.iter().enumerate() {
        let prev_upper = i > 0 && (b[i - 1] == b'C' || b[i - 1] == b'B'); // "Cl", "Br" : une seule unité
        match c {
            b'A'..=b'Z' if c != b'H' => n += 1,
            b'c' | b'n' | b'o' | b's' | b'p' | b'b' if !prev_upper => n += 1,
            _ => {}
        }
    }
    n
}

struct Part {
    name: String,
    path: std::path::PathBuf,
    dest: std::path::PathBuf,
    checkpoint: std::path::PathBuf,
    targets: Vec<Option<String>>,
    results: std::sync::Mutex<Vec<RowResult>>,
    remaining: AtomicUsize,
    writer: std::sync::Mutex<Option<std::io::BufWriter<fs::File>>>,
    written: AtomicBool,
}

/// Données partagées entre le thread principal et le thread de calcul.
struct Shared {
    parts: Vec<Part>,
    /// (fichier, ligne) dans l'ordre de calcul
    queue: Vec<(usize, usize)>,
    next: AtomicUsize,
    done: AtomicUsize,
    cancel: AtomicBool,
    /// Molécule en cours sur chaque thread : (début, fichier, ligne)
    running: std::sync::Mutex<Vec<Option<(Instant, usize, usize)>>>,
    slow_log: std::sync::Mutex<fs::File>,
}

/// Re-normalise tous les fichiers *.parquet de `input_dir` dans `output_dir` (mêmes noms de fichiers).
///
/// Une seule file d'attente pour toute la base : les fichiers sont pris dans l'ordre, et dans chaque fichier
/// les molécules les plus longues d'abord ; un thread libre prend toujours la molécule suivante, même si elle
/// est dans le fichier d'après (plus de cœurs inoccupés pendant les dernières molécules d'un fichier).
/// Un fichier est écrit dès que toutes ses molécules sont calculées.
///
/// `progress` est appelé toutes les 0,5 s ; s'il renvoie false, le calcul s'arrête (ce qui est calculé est
/// gardé dans les fichiers de reprise). `only` : noms des fichiers à traiter (None = tous).
pub fn renormalize_database(
    input_dir: &Path,
    output_dir: &Path,
    only: Option<&[String]>,
    log: &dyn Fn(&str),
    progress: &dyn Fn(&Progress) -> bool,
) -> Result<RenormStats, String> {
    use std::io::Write;
    if input_dir == output_dir { return Err("le dossier de sortie doit être différent du dossier d'entrée".into()); }
    fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;
    let mut paths: Vec<_> = fs::read_dir(input_dir).map_err(|e| e.to_string())?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().map_or(false, |x| x == "parquet"))
        .filter(|p| only.map_or(true, |o| o.iter().any(|n| Some(n.as_str()) == p.file_name().and_then(|f| f.to_str()))))
        .collect();
    paths.sort();
    if paths.is_empty() { return Err(format!("aucun fichier .parquet dans {}", input_dir.display())); }
    let files_total = paths.len();
    let t0 = Instant::now();
    let mut total = RenormStats::default();

    // --- Préparation : fichiers déjà faits ignorés, fichiers sans structure recopiés, reprises relues ---
    let mut parts: Vec<Part> = Vec::new();
    let mut files_done = 0;
    for path in &paths {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        let dest = output_dir.join(&name);
        if dest.exists() {
            log(&format!("{name} : déjà traité, ignoré"));
            files_done += 1;
            continue;
        }
        let mut file = fs::File::open(path).map_err(|e| e.to_string())?;
        let df = ParquetReader::new(&mut file).finish().map_err(|e| format!("{}: {e}", path.display()))?;
        let Some(targets) = targets_of(&df) else {
            let tmp = output_dir.join(format!("{name}.en_cours"));
            fs::copy(path, &tmp).map_err(|e| e.to_string())?;
            fs::rename(&tmp, &dest).map_err(|e| e.to_string())?;
            log(&format!("{name} : pas de structure, recopié tel quel ({} lignes)", df.height()));
            total.rows += df.height();
            files_done += 1;
            continue;
        };
        drop(df);
        let checkpoint = output_dir.join(format!("{name}.reprise"));
        let mut results: Vec<RowResult> = vec![None; targets.len()];
        let resumed = read_checkpoint(&checkpoint, &mut results);
        let remaining = (0..targets.len()).filter(|&i| targets[i].is_some() && results[i].is_none()).count();
        if resumed > 0 {
            log(&format!("{name} : reprise, {resumed} molécules déjà calculées, {remaining} à faire"));
        }
        let writer = fs::OpenOptions::new().create(true).append(true).open(&checkpoint)
            .map(|f| std::io::BufWriter::with_capacity(1 << 16, f)).map_err(|e| e.to_string())?;
        parts.push(Part {
            name, path: path.clone(), dest, checkpoint, targets,
            results: std::sync::Mutex::new(results),
            remaining: AtomicUsize::new(remaining),
            writer: std::sync::Mutex::new(Some(writer)),
            written: AtomicBool::new(false),
        });
    }

    // --- File d'attente : d'abord toutes les molécules ordinaires (fichiers dans l'ordre, les plus longues
    //     d'abord dans chaque fichier), puis les très grosses molécules de toute la base, de la plus petite à la
    //     plus grosse (les fichiers qui en ont peu sont écrits plus tôt) ---
    let big = big_atoms_threshold();
    let mut queue: Vec<(usize, usize)> = Vec::new();
    let mut big_queue: Vec<(usize, usize, usize)> = Vec::new();
    for (p, part) in parts.iter().enumerate() {
        let res = part.results.lock().unwrap();
        let mut rows: Vec<usize> = Vec::new();
        for i in (0..part.targets.len()).filter(|&i| part.targets[i].is_some() && res[i].is_none()) {
            let atoms = heavy_atoms(part.targets[i].as_ref().unwrap());
            if atoms > big { big_queue.push((atoms, p, i)); } else { rows.push(i); }
        }
        rows.sort_by(|&a, &b| part.targets[b].as_ref().unwrap().len().cmp(&part.targets[a].as_ref().unwrap().len()));
        queue.extend(rows.into_iter().map(|r| (p, r)));
    }
    if !big_queue.is_empty() {
        log(&format!("{} très grosses molécules (plus de {big} atomes lourds) gardées pour la fin", big_queue.len()));
    }
    big_queue.sort();
    queue.extend(big_queue.into_iter().map(|(_, p, r)| (p, r)));
    let todo = queue.len();
    let workers = rayon::current_num_threads().max(1);
    // État partagé avec le thread de calcul. Ce thread n'est pas attendu à l'arrêt (Ctrl+C) : une molécule
    // géante en cours peut durer des heures ; on enregistre ce qui est calculé et on rend la main tout de suite.
    let shared = std::sync::Arc::new(Shared {
        parts, queue,
        next: AtomicUsize::new(0),
        done: AtomicUsize::new(0),
        cancel: AtomicBool::new(false),
        running: std::sync::Mutex::new(vec![None; workers + 1]),
        slow_log: std::sync::Mutex::new(
            fs::OpenOptions::new().create(true).append(true).open(output_dir.join("_molecules_lentes.tsv")).map_err(|e| e.to_string())?),
    });
    let parts = &shared.parts;

    // Écrit un fichier terminé (sur le thread principal)
    let write_part = |part: &Part, total: &mut RenormStats| -> Result<(), String> {
        if let Some(mut w) = part.writer.lock().unwrap().take() { let _ = w.flush(); }
        let t = Instant::now();
        let mut file = fs::File::open(&part.path).map_err(|e| e.to_string())?;
        let df = ParquetReader::new(&mut file).finish().map_err(|e| e.to_string())?;
        let results = std::mem::take(&mut *part.results.lock().unwrap());
        let (mut out, st) = apply_results(&df, &results).map_err(|e| e.to_string())?;
        let tmp = output_dir.join(format!("{}.en_cours", part.name));
        let f = fs::File::create(&tmp).map_err(|e| e.to_string())?;
        ParquetWriter::new(f).with_compression(ParquetCompression::Zstd(None)).finish(&mut out).map_err(|e| e.to_string())?;
        fs::rename(&tmp, &part.dest).map_err(|e| e.to_string())?;
        let _ = fs::remove_file(&part.checkpoint);
        part.written.store(true, Ordering::Relaxed);
        log(&format!("{} : écrit - {} lignes, {} avec structure, {} InChIKey modifiées, {} SMILES modifiés, {} illisibles (gardées) ({:.0} s d'écriture)",
            part.name, st.rows, st.with_structure, st.key_changed, st.smiles_changed, st.failed, t.elapsed().as_secs_f64()));
        total.add(&st);
        Ok(())
    };

    let mut stopped = false;
    let mut error: Option<String> = None;
    let worker = {
        let sh = shared.clone();
        std::thread::spawn(move || {
            (0..workers).into_par_iter().for_each(|_| {
                let slot = rayon::current_thread_index().unwrap_or(workers);
                while !sh.cancel.load(Ordering::Relaxed) {
                    let k = sh.next.fetch_add(1, Ordering::Relaxed);
                    let Some(&(p, row)) = sh.queue.get(k) else { break };
                    let part = &sh.parts[p];
                    let target = part.targets[row].as_ref().unwrap();
                    let start = Instant::now();
                    sh.running.lock().unwrap()[slot] = Some((start, p, row));
                    let out = compute_one(target);
                    sh.running.lock().unwrap()[slot] = None;
                    let secs = start.elapsed().as_secs_f64();
                    if secs >= SLOW_SECS {
                        let _ = writeln!(sh.slow_log.lock().unwrap(), "{secs:.0}\t{}\t{row}\t{target}", part.name);
                    }
                    if let Some(w) = part.writer.lock().unwrap().as_mut() {
                        let _ = w.write_all(checkpoint_line(row, &out).as_bytes());
                    }
                    part.results.lock().unwrap()[row] = Some(out);
                    part.remaining.fetch_sub(1, Ordering::Relaxed);
                    sh.done.fetch_add(1, Ordering::Relaxed);
                }
            });
        })
    };
    // Thread principal : écrit les fichiers terminés, affiche l'avancement, gère l'arrêt
    loop {
        let finished = worker.is_finished();
        for part in parts {
            if !part.written.load(Ordering::Relaxed) && part.remaining.load(Ordering::Relaxed) == 0 && error.is_none() {
                if let Err(e) = write_part(part, &mut total) { error = Some(e); shared.cancel.store(true, Ordering::Relaxed); }
            }
            // Le fichier de reprise est vidé sur le disque régulièrement
            if let Some(w) = part.writer.lock().unwrap().as_mut() { let _ = w.flush(); }
        }
        let slow: Vec<(f64, String, String)> = shared.running.lock().unwrap().iter().flatten()
            .filter(|(s, _, _)| s.elapsed().as_secs_f64() >= SLOW_SECS)
            .map(|(s, p, r)| (s.elapsed().as_secs_f64(), parts[*p].name.clone(), parts[*p].targets[*r].clone().unwrap_or_default()))
            .collect();
        let written = files_done + parts.iter().filter(|p| p.written.load(Ordering::Relaxed)).count();
        let info = Progress { done: shared.done.load(Ordering::Relaxed), total: todo, files_done: written, files_total, slow };
        if !progress(&info) {
            // Arrêt demandé : on n'attend pas les molécules en cours (le thread de calcul s'arrêtera de lui-même
            // après elles, ou avec la fin du programme) ; tout ce qui est calculé est dans les fichiers de reprise.
            stopped = true;
            shared.cancel.store(true, Ordering::Relaxed);
            break;
        }
        if finished || error.is_some() { break; }
        std::thread::sleep(Duration::from_millis(500));
    }
    if !stopped && error.is_none() {
        if let Err(panic) = worker.join() { std::panic::resume_unwind(panic); }
    }
    for part in parts {
        // Arrêt : on vide le fichier de reprise sans le fermer (le thread de calcul peut encore y écrire)
        if let Some(w) = part.writer.lock().unwrap().as_mut() { let _ = w.flush(); }
    }
    if let Some(e) = error { return Err(e); }
    if stopped { return Err("interrompu".into()); }

    // InChIKey en double après normalisation (sur toute la base)
    let mut seen = std::collections::HashMap::<String, usize>::new();
    for path in &paths {
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
    let last_line = std::sync::Mutex::new(Instant::now());
    let last_slow = std::sync::Mutex::new(Instant::now());
    let one_line = prefix.is_empty();
    let r = py.allow_threads(|| renormalize_database(
        Path::new(input_dir),
        Path::new(output_dir),
        only.as_deref(),
        &|m| println!("\r{prefix}{m}{}", " ".repeat(30)),
        // Progression, molécules lentes, et Ctrl+C vérifié toutes les 0,5 s (Python reprend la main un instant)
        &|p: &Progress| {
            let minutes = t0.elapsed().as_secs_f64() / 60.0;
            let rate = if minutes > 0.0 { p.done as f64 / minutes } else { 0.0 };
            let left = if rate > 0.0 { (p.total - p.done) as f64 / rate } else { f64::NAN };
            let mut text = format!("{prefix}{}/{} fichiers écrits, {}/{} molécules, {rate:.0} molécules/min, reste ~{left:.0} min ({minutes:.0} min écoulées)",
                p.files_done, p.files_total, p.done, p.total);
            if !p.slow.is_empty() {
                let longest = p.slow.iter().map(|s| s.0).fold(0.0, f64::max);
                text.push_str(&format!(", {} molécule(s) longue(s) en cours (la plus longue depuis {:.0} min)", p.slow.len(), longest / 60.0));
            }
            if one_line {
                print!("\r{text}   ");
            } else {
                let mut last = last_line.lock().unwrap();
                if last.elapsed() >= Duration::from_secs(60) { println!("{text}"); *last = Instant::now(); }
            }
            // Toutes les 10 minutes : détail des molécules longues en cours
            let mut ls = last_slow.lock().unwrap();
            if !p.slow.is_empty() && ls.elapsed() >= Duration::from_secs(600) {
                println!();
                for (secs, file, mol) in &p.slow {
                    let short: String = mol.chars().take(100).collect();
                    println!("{prefix}  en cours depuis {:.0} min ({file}) : {short}{}", secs / 60.0, if mol.chars().count() > 100 { "..." } else { "" });
                }
                *ls = Instant::now();
            }
            let _ = std::io::stdout().flush();
            match Python::with_gil(|py| py.check_signals()) {
                Ok(()) => true,
                Err(e) => { *interrupted.lock().unwrap() = Some(e); false }
            }
        },
    ));
    if let Some(e) = interrupted.into_inner().unwrap() {
        println!("\nArrêté. Ce qui est calculé est gardé (fichiers .reprise) : relancer le script reprend là où il s'était arrêté.");
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
        let (out, st) = renormalize_frame(&df).unwrap().unwrap();
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
        let st = renormalize_database(&src, &dst, None, &|_| {}, &|_: &Progress| true).unwrap();
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
    fn reprise_a_la_molecule_pres() {
        let dir = std::env::temp_dir().join(format!("fraghub_renorm_stop_{}", std::process::id()));
        let (src, dst, ref_) = (dir.join("in"), dir.join("out"), dir.join("ref"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&src).unwrap();
        for name in ["part_0000.parquet", "part_0001.parquet"] {
            let mut df = frame();
            ParquetWriter::new(fs::File::create(src.join(name)).unwrap()).finish(&mut df).unwrap();
        }
        // Référence : traitement d'un coup
        renormalize_database(&src, &ref_, None, &|_| {}, &|_: &Progress| true).unwrap();

        // Arrêt demandé dès le début : erreur « interrompu », aucun fichier .en_cours laissé à moitié
        let r = renormalize_database(&src, &dst, None, &|_| {}, &|_: &Progress| false);
        assert!(r.is_err());
        assert!(fs::read_dir(&dst).unwrap().all(|e| !e.unwrap().file_name().to_string_lossy().ends_with(".en_cours")));
        // Fichier de reprise fabriqué : la ligne 0 du 2e fichier est « déjà calculée » avec un SMILES repère ;
        // il doit être repris tel quel (preuve que la reprise évite le recalcul)
        let _ = fs::remove_file(dst.join("part_0001.parquet"));
        fs::write(dst.join("part_0001.parquet.reprise"),
            format!("{}{{coupé", checkpoint_line(0, &[("INCHIKEY", "OKJIRPAQVSHGFK-UHFFFAOYSA-N".into()), ("SMILES", "REPERE".into())]))).unwrap();
        let st = renormalize_database(&src, &dst, None, &|_| {}, &|_: &Progress| true).unwrap();
        assert!(st.rows <= 10);
        assert!(!dst.join("part_0001.parquet.reprise").exists(), "fichier de reprise supprimé après écriture");
        let read = |p: &Path| { let mut f = fs::File::open(p).unwrap(); ParquetReader::new(&mut f).finish().unwrap() };
        let smi = read(&dst.join("part_0001.parquet")).column("SMILES").unwrap().str().unwrap().get(0).map(|s| s.to_string());
        assert_eq!(smi.as_deref(), Some("REPERE"));
        // Le 1er fichier, lui, est identique à la référence
        assert!(read(&dst.join("part_0000.parquet")).equals_missing(&read(&ref_.join("part_0000.parquet"))));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn arret_sans_attendre_une_molecule_geante() {
        // Macrocycle C99H90O54 de la base : environ 270 s de calcul du tautomère. L'arrêt doit rendre la main
        // tout de suite, sans attendre la fin de ce calcul.
        let dir = std::env::temp_dir().join(format!("fraghub_renorm_geant_{}", std::process::id()));
        let (src, dst) = (dir.join("in"), dir.join("out"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&src).unwrap();
        let mut df = DataFrame::new(vec![
            Series::new("INCHIKEY", &[Some("AAAAAAAAAAAAAA-UHFFFAOYSA-N")]),
            Series::new("INCHI", &[Some("InChI=1S/C99H90O54/c100-82(101)28-136-64-12-48-2-50-14-70(142-34-88(112)113)52(16-68(50)140-32-86(108)109)4-54-18-74(146-38-92(120)121)56(20-72(54)144-36-90(116)117)6-58-22-78(150-42-96(128)129)60(24-76(58)148-40-94(124)125)8-62-26-81(153-45-99(134)135)63(27-80(62)152-44-98(132)133)9-61-25-77(149-41-95(126)127)59(23-79(61)151-43-97(130)131)7-57-21-73(145-37-91(118)119)55(19-75(57)147-39-93(122)123)5-53-17-69(141-33-87(110)111)51(15-71(53)143-35-89(114)115)3-49-13-65(137-29-83(102)103)47(11-67(49)139-31-85(106)107)1-46(64)10-66(48)138-30-84(104)105/h10-27H,1-9,28-45H2,(H,100,101)(H,102,103)(H,104,105)(H,106,107)(H,108,109)(H,110,111)(H,112,113)(H,114,115)(H,116,117)(H,118,119)(H,120,121)(H,122,123)(H,124,125)(H,126,127)(H,128,129)(H,130,131)(H,132,133)(H,134,135)")]),
            Series::new("SMILES", &[None::<&str>]),
        ]).unwrap();
        ParquetWriter::new(fs::File::create(src.join("part_0000.parquet")).unwrap()).finish(&mut df).unwrap();
        let t = Instant::now();
        let calls = AtomicUsize::new(0);
        let r = renormalize_database(&src, &dst, None, &|_| {}, &|_: &Progress| calls.fetch_add(1, Ordering::Relaxed) < 4);
        assert!(r.is_err());
        assert!(t.elapsed() < Duration::from_secs(20), "arrêt en {:?}", t.elapsed());
        assert!(!dst.join("part_0000.parquet").exists());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn ordre_et_resultats_identiques_au_calcul_direct() {
        // La file d'attente commune à plusieurs fichiers donne exactement le même résultat que le calcul direct
        let dir = std::env::temp_dir().join(format!("fraghub_renorm_same_{}", std::process::id()));
        let (src, dst) = (dir.join("in"), dir.join("out"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&src).unwrap();
        let mut df = frame();
        ParquetWriter::new(fs::File::create(src.join("part_0000.parquet")).unwrap()).finish(&mut df).unwrap();
        renormalize_database(&src, &dst, None, &|_| {}, &|_: &Progress| true).unwrap();
        let (direct, _) = renormalize_frame(&frame()).unwrap().unwrap();
        let mut f = fs::File::open(dst.join("part_0000.parquet")).unwrap();
        assert!(ParquetReader::new(&mut f).finish().unwrap().equals_missing(&direct));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn atomes_lourds_et_grosses_molecules_a_la_fin() {
        assert_eq!(heavy_atoms("InChI=1S/C216H184N16O8/c1-137"), 240);
        assert_eq!(heavy_atoms("InChI=1S/C10H12.2ClH/c1"), 12);
        assert_eq!(heavy_atoms("InChI=1S/H2O/h1H2"), 1);
        assert_eq!(heavy_atoms("InChI=1S/C5H12O5/c6-1"), 10);
        assert_eq!(heavy_atoms("OCC(O)C(O)C(O)CO"), 10);
        assert_eq!(heavy_atoms("c1ccccc1Cl"), 7);
        assert_eq!(heavy_atoms("[Na+].[Cl-]"), 2);
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
        renormalize_database(&src, &dst, Some(&only), &|_| {}, &|_: &Progress| true).unwrap();
        assert!(!dst.join("part_0000.parquet").exists());
        assert!(dst.join("part_0001.parquet").exists());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn refuse_d_ecrire_dans_le_dossier_d_entree() {
        let d = std::env::temp_dir();
        assert!(renormalize_database(&d, &d, None, &|_| {}, &|_: &Progress| true).is_err());
    }
}

/// Essai sur une vraie base : FRAGHUB_RENORM_IN=... FRAGHUB_RENORM_OUT=... cargo test --release -- --ignored vraie_base --nocapture
#[cfg(test)]
#[test]
#[ignore]
fn vraie_base() {
    let (i, o) = (std::env::var("FRAGHUB_RENORM_IN").unwrap(), std::env::var("FRAGHUB_RENORM_OUT").unwrap());
    renormalize_database(Path::new(&i), Path::new(&o), None, &|m| println!("{m}"), &|_: &Progress| true).unwrap();
}

