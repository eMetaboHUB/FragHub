// src/rdkit_native/mod.rs
//! Calcul RDKit directement en C++ (bibliothèques précompilées de conda-forge), sans Python.
//!
//! Pour un développeur Python : avant, Rust envoyait toutes les molécules à un script Python
//! (`rdkit_worker.py`, supprimé depuis), qui les répartissait sur ~20 processus Python (multiprocessing).
//! Ici, Rust appelle RDKit en C++ via une petite passerelle (shim.cc) qui reproduit EXACTEMENT les calculs
//! de cet ancien script (sorties vérifiées identiques), et répartit le travail sur des threads (rayon) :
//! un seul processus, pas de copie des molécules entre processus.
//!
//! `run_parallel` renvoie le dictionnaire :
//! clé "molécule|adduit|m/z" -> {INCHI, INCHIKEY, SMILES, FORMULA, EXACTMASS, AVERAGEMASS}
//! (une clé est absente si le calcul correspondant a échoué).
//!
//! Deux optimisations qui ne changent aucun résultat :
//! - une molécule sans « . » (ou sans m/z lisible) donne le même résultat quels que soient l'adduit et le m/z
//!   (ils ne servent qu'à choisir un fragment) : elle n'est calculée qu'une fois ;
//! - répartition par FILE D'ATTENTE PARTAGÉE : chaque thread prend la molécule suivante dès qu'il est libre.
//!   (Avant : `par_iter` de rayon, qui découpe la liste en gros blocs ; mesuré sur 189 752 molécules et
//!   12 threads, des threads finissaient seuls leur bloc : 311 s au lieu des 201 s possibles.) Les molécules
//!   les plus longues partent toujours en premier.

use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Once;
use std::time::{Duration, Instant};

#[cxx::bridge(namespace = "fraghub")]
mod ffi {
    #[derive(Default)]
    struct MolResult {
        has_inchi: bool,
        inchi: String,
        has_inchikey: bool,
        inchikey: String,
        has_smiles: bool,
        smiles: String,
        has_formula: bool,
        formula: String,
        has_exactmass: bool,
        exactmass: f64,
        has_averagemass: bool,
        averagemass: f64,
    }

    unsafe extern "C++" {
        include!("fraghub_rust/src/rdkit_native/shim.h");
        fn init_rdkit();
        fn process_single(
            clean_mol: &str,
            ptype: &str,
            precursormz: f64,
            has_mz: bool,
            adduct_known: bool,
            adduct_mult: i32,
            adduct_delta: f64,
            adduct_charge: i32,
        ) -> MolResult;
        fn exact_mass(mol: &str, is_inchi: bool) -> f64;
    }
}

static INIT: Once = Once::new();

/// Équivalent de `str(x)` en Python pour un nombre décimal (mêmes chiffres, même notation).
pub fn py_float_str(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    // Rust ({:?}) et Python (repr) donnent tous deux le plus court nombre de chiffres qui relit
    // exactement la même valeur. Python écrit en notation fixe si 1e-4 <= |x| < 1e16.
    let a = x.abs();
    if x == 0.0 || (1e-4..1e16).contains(&a) {
        let s = format!("{:?}", x);
        return if s.contains('e') { python_fixed(&s) } else { s };
    }
    let s = format!("{:e}", x);
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    format!("{}e{}{:02}", mant, if exp < 0 { '-' } else { '+' }, exp.abs())
}

fn python_fixed(sci: &str) -> String {
    let (mant, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let neg = mant.starts_with('-');
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let point = 1 + exp;
    let mut out = String::new();
    if point <= 0 {
        out.push_str("0.");
        out.push_str(&"0".repeat((-point) as usize));
        out.push_str(&digits);
    } else if point as usize >= digits.len() {
        out.push_str(&digits);
        out.push_str(&"0".repeat(point as usize - digits.len()));
        out.push_str(".0");
    } else {
        out.push_str(&digits[..point as usize]);
        out.push('.');
        out.push_str(&digits[point as usize..]);
    }
    if neg { format!("-{out}") } else { out }
}

/// Même règle que `str.strip()` de Python pour les espaces.
fn py_strip(s: &str) -> &str {
    s.trim_matches(|c: char| c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c))
}

/// `float(precursormz_str)` en Python (FragHub fournit "" ou un nombre formaté par Rust).
fn py_float(s: &str) -> Option<f64> {
    s.trim().parse::<f64>().ok()
}

/// Calcul d'une molécule (équivalent exact de l'ancien `rdkit_worker.process_single`) : (clé, valeur).
pub fn process_single(clean_mol: &str, precursortype: &str, precursormz_str: &str) -> Vec<(&'static str, String)> {
    INIT.call_once(ffi::init_rdkit);
    let mz = py_float(precursormz_str);
    let ptype = py_strip(precursortype);
    // Adduit pour le dessalage : paramètres lus dans adduct_to_convert.csv (formes FragHub) ; à défaut,
    // le code C++ cherche dans sa propre table.
    let (known, mult, delta, charge) = adduct_params(ptype).map_or((false, 0, 0.0, 0), |(n, d, z)| (true, n, d, z));
    let r = ffi::process_single(clean_mol, ptype, mz.unwrap_or(f64::NAN), mz.is_some(), known, mult, delta, charge);
    let mut out = Vec::with_capacity(6);
    if r.has_inchi { out.push(("INCHI", r.inchi)); }
    if r.has_inchikey { out.push(("INCHIKEY", r.inchikey)); }
    if r.has_smiles { out.push(("SMILES", r.smiles)); }
    if r.has_formula { out.push(("FORMULA", r.formula)); }
    if r.has_exactmass { out.push(("EXACTMASS", py_float_str(r.exactmass))); }
    if r.has_averagemass { out.push(("AVERAGEMASS", py_float_str(r.averagemass))); }
    out
}

/// (n molécules, différence de masse totale, charge) d'un adduit au format FragHub, d'après les tables
/// chargées depuis adduct_to_convert.csv. None si l'adduit n'y est pas.
fn adduct_params(ptype: &str) -> Option<(i32, f64, i32)> {
    if ptype.is_empty() { return None; }
    let state = crate::global_state::STATE.read().ok()?;
    let delta = state.adduct_massdiff_dict_pos.get(ptype).or_else(|| state.adduct_massdiff_dict_neg.get(ptype)).copied()?;
    let (n, z) = crate::adduct_mass_calculator::adduct_nmol_and_charge(ptype)?;
    Some((n as i32, delta, z as i32))
}

/// Masse exacte d'une molécule, même choix que l'ancienne version Python de
/// `missing_precursormz_re_calculation` : l'InChI s'il contient "InChI=", sinon le SMILES s'il n'est pas vide.
/// `None` si la molécule ne se lit pas (Python : MolFromInchi / MolFromSmiles renvoie None).
pub fn exact_mass(inchi: &str, smiles: &str) -> Option<f64> {
    INIT.call_once(ffi::init_rdkit);
    let m = if !inchi.is_empty() && inchi.contains("InChI=") {
        ffi::exact_mass(inchi, true)
    } else if !smiles.is_empty() {
        ffi::exact_mass(smiles, false)
    } else {
        return None;
    };
    if m.is_nan() { None } else { Some(m) }
}

/// Fonction de vérification appelable depuis Python : `fraghub_rust.native_exact_mass(inchi, smiles)`.
/// Sert à comparer le C++ avec `Descriptors.ExactMolWt` de RDKit Python (voir scratch/).
#[pyfunction]
pub fn native_exact_mass(inchi: &str, smiles: &str) -> Option<f64> {
    exact_mass(inchi, smiles)
}

/// Le résultat dépend-il de l'adduit et du m/z ? Seulement pour une molécule à plusieurs fragments
/// (« . ») avec un m/z lisible : c'est la seule branche de process_single qui les utilise.
fn depends_on_precursor(clean_mol: &str, precursormz_str: &str) -> bool {
    clean_mol.contains('.') && py_float(precursormz_str).is_some()
}

/// Calcule toutes les molécules et renvoie le dictionnaire (voir en tête de fichier).
/// La barre de progression avance en nombre de molécules d'entrée (comme avant).
/// Si `progress_callback` lève une exception (bouton stop), le calcul s'arrête et l'erreur est renvoyée.
pub fn run_parallel(
    py: Python,
    mols_list: Vec<(String, String, String)>,
    progress_callback: Option<PyObject>,
) -> PyResult<HashMap<String, HashMap<String, String>>> {
    let total = mols_list.len();
    let cache = compute_cache(
        mols_list,
        |done| match &progress_callback {
            Some(cb) => cb.call1(py, (done,)).map(|_| ()),
            None => Ok(()),
        },
        || py.allow_threads(|| std::thread::sleep(Duration::from_millis(100))),
    )?;
    if let Some(cb) = &progress_callback {
        cb.call1(py, (total,))?;
    }
    Ok(cache)
}

/// Cœur du calcul, sans Python : `tick(n)` est appelé régulièrement avec le nombre de molécules
/// d'entrée traitées (une erreur arrête le calcul) ; `wait()` fait patienter entre deux appels.
fn compute_cache<E>(
    mols_list: Vec<(String, String, String)>,
    mut tick: impl FnMut(usize) -> Result<(), E>,
    wait: impl Fn(),
) -> Result<HashMap<String, HashMap<String, String>>, E> {
    INIT.call_once(ffi::init_rdkit);

    // Unités de calcul : une par molécule, sauf quand l'adduit et le m/z comptent (mélanges de fragments)
    let mut unit_index: HashMap<(&str, &str, &str), usize> = HashMap::new();
    let mut units: Vec<(&str, &str, &str)> = Vec::new();
    let mut weights: Vec<usize> = Vec::new();
    let mut unit_of: Vec<usize> = Vec::with_capacity(mols_list.len());
    for (m, pt, mz) in &mols_list {
        let key = if depends_on_precursor(m, mz) { (m.as_str(), pt.as_str(), mz.as_str()) } else { (m.as_str(), "", "") };
        let id = *unit_index.entry(key).or_insert_with(|| {
            units.push((m.as_str(), pt.as_str(), mz.as_str()));
            weights.push(0);
            units.len() - 1
        });
        weights[id] += 1;
        unit_of.push(id);
    }
    drop(unit_index);

    // Les plus longues d'abord (approximation de la difficulté), ordre stable sinon
    let mut order: Vec<usize> = (0..units.len()).collect();
    order.sort_by(|&a, &b| units[b].0.len().cmp(&units[a].0.len()));

    let done = AtomicUsize::new(0);
    let cancel = AtomicBool::new(false);
    // Mesure facultative (FRAGHUB_RDKIT_TIMING=<dossier>) : début et durée de chaque molécule, sans effet sur les résultats
    let timing_dir = std::env::var("FRAGHUB_RDKIT_TIMING").ok().filter(|d| !d.trim().is_empty());
    let t0 = Instant::now();
    let mut results: Vec<Option<Vec<(&'static str, String)>>> = Vec::new();
    let mut tick_error: Option<E> = None;

    std::thread::scope(|scope| {
        let worker = scope.spawn(|| {
            // File d'attente partagée : `next` est le rang de la prochaine molécule à calculer dans `order`.
            // Un travailleur par thread de rayon ; chacun prend la suivante dès qu'il a fini la sienne.
            let next = AtomicUsize::new(0);
            let workers = rayon::current_num_threads().max(1);
            let computed: Vec<(usize, Vec<(&'static str, String)>, Timing)> = (0..workers)
                .into_par_iter()
                .flat_map_iter(|_| {
                    let mut mine = Vec::new();
                    while !cancel.load(Ordering::Relaxed) {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(&id) = order.get(i) else { break };
                        let (m, pt, mz) = units[id];
                        let start = t0.elapsed();
                        let r = process_single(m, pt, mz);
                        let duration = t0.elapsed() - start;
                        done.fetch_add(weights[id], Ordering::Relaxed);
                        mine.push((id, r, Timing { thread: rayon::current_thread_index().unwrap_or(usize::MAX), start, duration }));
                    }
                    mine
                })
                .collect();
            if let Some(dir) = &timing_dir {
                write_timing(dir, &order, &units, &weights, &computed);
            }
            let mut by_unit: Vec<Option<Vec<(&'static str, String)>>> = vec![None; units.len()];
            for (id, r, _) in computed {
                by_unit[id] = Some(r);
            }
            by_unit
        });
        // Ce thread-ci (qui a le droit d'appeler Python) fait avancer la barre
        let mut last = usize::MAX;
        while !worker.is_finished() {
            wait();
            let d = done.load(Ordering::Relaxed);
            if d != last && tick_error.is_none() {
                last = d;
                if let Err(e) = tick(d) {
                    cancel.store(true, Ordering::Relaxed); // bouton stop : on arrête le calcul
                    tick_error = Some(e);
                }
            }
        }
        match worker.join() {
            Ok(r) => results = r,
            Err(panic) => std::panic::resume_unwind(panic),
        }
    });
    if let Some(e) = tick_error {
        return Err(e);
    }

    // Même dictionnaire que la version Python : clé "molécule|adduit|m/z"
    let mut cache = HashMap::with_capacity(mols_list.len());
    for ((m, pt, mz), id) in mols_list.iter().zip(unit_of) {
        let transforms: HashMap<String, String> = results[id]
            .as_ref()
            .map(|r| r.iter().map(|(k, v)| (k.to_string(), v.clone())).collect())
            .unwrap_or_default();
        cache.insert(format!("{}|{}|{}", m, pt, mz), transforms);
    }
    Ok(cache)
}

/// Mesure d'une molécule : thread rayon, début et durée (depuis le début du calcul).
struct Timing {
    thread: usize,
    start: Duration,
    duration: Duration,
}

static TIMING_PASS: AtomicUsize = AtomicUsize::new(0);

/// Écrit `<dossier>/rdkit_timing_<date>_passN.tsv` : une ligne par molécule calculée, par ordre de début.
/// Colonnes : rang dans l'ordre de traitement, thread, début (s), durée (s), nombre de spectres qui
/// utilisent ce résultat, longueur du texte, molécule. Une erreur d'écriture est signalée sans arrêter FragHub.
fn write_timing(dir: &str, order: &[usize], units: &[(&str, &str, &str)], weights: &[usize],
                computed: &[(usize, Vec<(&'static str, String)>, Timing)]) {
    use std::io::Write;
    let pass = TIMING_PASS.fetch_add(1, Ordering::Relaxed) + 1;
    let date = chrono::Local::now().format("%Y%m%d-%H%M%S");
    let path = std::path::Path::new(dir).join(format!("rdkit_timing_{date}_pass{pass}.tsv"));
    let mut rank = vec![0usize; units.len()];
    for (r, &id) in order.iter().enumerate() {
        rank[id] = r;
    }
    let mut rows: Vec<&(usize, Vec<(&'static str, String)>, Timing)> = computed.iter().collect();
    rows.sort_by_key(|(_, _, t)| t.start);
    let write = || -> std::io::Result<()> {
        std::fs::create_dir_all(dir)?;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
        writeln!(f, "# threads={} unites={}", rayon::current_num_threads(), units.len())?;
        writeln!(f, "rang\tthread\tdebut_s\tduree_s\tpoids\tlongueur\tmolecule")?;
        for (id, _, t) in rows {
            let m = units[*id].0.replace(&['\t', '\n', '\r'][..], " ");
            writeln!(f, "{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}", rank[*id], t.thread, t.start.as_secs_f64(),
                     t.duration.as_secs_f64(), weights[*id], units[*id].0.len(), m)?;
        }
        f.flush()
    };
    match write() {
        Ok(()) => println!("[RDKIT_TIMING] {}", path.display()),
        Err(e) => eprintln!("[RDKIT_TIMING] écriture impossible ({}) : {e}", path.display()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ecriture_des_nombres_comme_python() {
        for (x, py) in [(12.0, "12.0"), (190.120509056, "190.120509056"), (0.0, "0.0"), (-0.0, "-0.0"),
                        (1e-5, "1e-05"), (1e16, "1e+16"), (1.5e-7, "1.5e-07"), (123456789012345.6, "123456789012345.6"),
                        (0.0001, "0.0001"), (591.1719292839999, "591.1719292839999"), (f64::NAN, "nan")] {
            assert_eq!(py_float_str(x), py, "{x}");
        }
    }

    #[test]
    fn adduit_et_mz_ne_comptent_que_pour_les_melanges_lisibles() {
        assert!(!depends_on_precursor("CCO", "100.0000"));
        assert!(!depends_on_precursor("CCO.Cl", ""));
        assert!(depends_on_precursor("CCO.Cl", "100.0000"));
        assert!(depends_on_precursor("CCO.Cl", "NaN")); // float("NaN") réussit en Python
    }

    #[test]
    fn calcul_unique_meme_resultat_que_le_calcul_par_entree() {
        let mols = vec![
            ("CCO".to_string(), "[M+H]+".to_string(), "47.0491".to_string()),
            ("CCO".to_string(), "[M+Na]+".to_string(), "69.0311".to_string()),
            ("CCO".to_string(), "".to_string(), "".to_string()),
            ("CC(=O)[O-].[Na+]".to_string(), "[M+H]+".to_string(), "61.0284".to_string()),
            ("CC(=O)[O-].[Na+]".to_string(), "[M+Na]+".to_string(), "83.0104".to_string()),
            ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3".to_string(), "[M+H]+".to_string(), "".to_string()),
            ("pas une molécule".to_string(), "".to_string(), "".to_string()),
        ];
        let cache = compute_cache(mols.clone(), |_| Ok::<(), ()>(()), || std::thread::sleep(Duration::from_millis(5))).unwrap();
        assert_eq!(cache.len(), mols.len());
        for (m, pt, mz) in &mols {
            let attendu: HashMap<String, String> =
                process_single(m, pt, mz).into_iter().map(|(k, v)| (k.to_string(), v)).collect();
            assert_eq!(cache[&format!("{m}|{pt}|{mz}")], attendu, "{m} {pt} {mz}");
        }
        assert!(cache["pas une molécule||"].is_empty());
    }

    #[test]
    fn masse_exacte() {
        // Éthanol C2H6O : 46.041864814 (RDKit)
        let m = exact_mass("", "CCO").unwrap();
        assert!((m - 46.041864814).abs() < 1e-6, "{m}");
        assert_eq!(exact_mass("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "").unwrap(), m);
        assert_eq!(exact_mass("1S/C2H6O/c1-2-3", "CCO").unwrap(), m); // InChI sans préfixe : le SMILES est utilisé
        assert_eq!(exact_mass("", "pas une molécule"), None);
        assert_eq!(exact_mass("InChI=1S/n'importe quoi", "CCO"), None); // comme Python : pas de repli sur le SMILES
        assert_eq!(exact_mass("", ""), None);
    }

    #[test]
    fn mesure_ecrite_si_demandee_sans_changer_les_resultats() {
        let dir = std::env::temp_dir().join(format!("fraghub_timing_test_{}", std::process::id()));
        let mols: Vec<_> = ["CCO", "c1ccccc1O", "CC(=O)Nc1ccc(O)cc1"].iter()
            .map(|m| (m.to_string(), String::new(), String::new())).collect();
        let sans = compute_cache(mols.clone(), |_| Ok::<(), ()>(()), || std::thread::sleep(Duration::from_millis(2))).unwrap();
        std::env::set_var("FRAGHUB_RDKIT_TIMING", &dir);
        let avec = compute_cache(mols, |_| Ok::<(), ()>(()), || std::thread::sleep(Duration::from_millis(2))).unwrap();
        std::env::remove_var("FRAGHUB_RDKIT_TIMING");
        assert_eq!(sans, avec);
        // Un fichier avec l'en-tête et une ligne par molécule (d'autres tests peuvent en écrire en parallèle)
        let ok = std::fs::read_dir(&dir).unwrap().any(|e| {
            let t = std::fs::read_to_string(e.unwrap().path()).unwrap();
            t.lines().nth(1) == Some("rang\tthread\tdebut_s\tduree_s\tpoids\tlongueur\tmolecule")
                && t.lines().count() == 5 && t.contains("\tCC(=O)Nc1ccc(O)cc1")
        });
        let _ = std::fs::remove_dir_all(&dir);
        assert!(ok);
    }

    #[test]
    fn arret_demande_pendant_le_calcul() {
        let mols: Vec<_> = (0..2000).map(|i| (format!("C{}O", "C".repeat(i % 40)), String::new(), String::new())).collect();
        let r = compute_cache(mols, |_| Err::<(), &str>("stop"), || std::thread::sleep(Duration::from_millis(1)));
        assert_eq!(r.unwrap_err(), "stop");
    }

    fn get<'a>(out: &'a [(&'static str, String)], k: &str) -> &'a str {
        out.iter().find(|x| x.0 == k).map(|x| x.1.as_str()).unwrap_or("")
    }

    #[test]
    fn inchikey_calculee_avant_le_tautomere_stereo_gardee() {
        // L-ornithine : avant, le tautomère canonique effaçait la stéréo (clé ...-UHFFFAOYSA-N)
        let out = process_single("NCCC[C@H](N)C(=O)O", "[M+H]+", "133.0972");
        assert_eq!(get(&out, "INCHIKEY"), "AHLPHDHHMVZTML-BYPYZUCNSA-N");
        assert_eq!(get(&out, "SMILES"), "NCCC[C@H](N)C(=O)O");
    }

    #[test]
    fn tautomere_qui_change_la_molecule_ignore_pour_le_smiles() {
        // Énone de stéroïde : RDKit déplaçait la double liaison (autre composé). Clé et SMILES de la molécule lue.
        let smi = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]43C)[C@@H]1CC[C@@H]2O";
        let out = process_single(smi, "[M+H]+", "289.2162");
        assert_eq!(get(&out, "INCHIKEY"), "MUMGGOZAMZWBJJ-VIUKOLAESA-N");
        assert_eq!(get(&out, "SMILES"), smi);
    }

    #[test]
    fn tautomere_entre_heteroatomes_normalise_le_smiles() {
        // Iminol -> amide : même InChI standard, le SMILES est écrit sous la forme normalisée
        let a = process_single("CC(O)=NCC(=O)O", "", "");
        let b = process_single("CC(=O)NCC(=O)O", "", "");
        assert_eq!(get(&a, "SMILES"), "CC(=O)NCC(=O)O");
        assert_eq!(get(&a, "SMILES"), get(&b, "SMILES"));
        assert_eq!(get(&a, "INCHIKEY"), get(&b, "INCHIKEY"));
    }

    #[test]
    fn sel_ecrit_en_inchi_dessale() {
        // Sulfate d'agmatine (InChI à deux composants), mesuré en [M+H]+ à m/z 131,13 : on garde l'agmatine
        let inchi = "InChI=1S/C5H14N4.H2O4S/c6-3-1-2-4-9-5(7)8;1-5(2,3)4/h1-4,6H2,(H4,7,8,9);(H2,1,2,3,4)";
        let out = process_single(inchi, "[M+H]+", "131.13");
        assert_eq!(get(&out, "INCHIKEY"), "QYPPJABKJHAVHS-UHFFFAOYSA-N");
        assert_eq!(get(&out, "FORMULA"), "C5H14N4");
        // Sans m/z précurseur : pas de dessalage (comportement inchangé)
        let out = process_single(inchi, "[M+H]+", "");
        assert!(get(&out, "FORMULA").contains('S'));
    }
}

