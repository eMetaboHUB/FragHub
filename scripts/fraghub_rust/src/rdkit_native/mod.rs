// src/rdkit_native/mod.rs
//! Calcul RDKit directement en C++ (bibliothèques précompilées de conda-forge), sans Python.
//!
//! Pour un développeur Python : avant, Rust envoyait toutes les molécules à `rdkit_worker.py`, qui les
//! répartissait sur ~20 processus Python (multiprocessing). Ici, Rust appelle RDKit en C++ via une petite
//! passerelle (shim.cc) qui reproduit EXACTEMENT `rdkit_worker.process_single`, et répartit le travail sur
//! des threads (rayon) : un seul processus, pas de copie des molécules entre processus.
//!
//! `run_parallel` renvoie le même dictionnaire que `rdkit_worker.run_parallel` :
//! clé "molécule|adduit|m/z" -> {INCHI, INCHIKEY, SMILES, FORMULA, EXACTMASS, AVERAGEMASS}
//! (une clé est absente si le calcul correspondant a échoué).
//!
//! Deux optimisations qui ne changent aucun résultat :
//! - une molécule sans « . » (ou sans m/z lisible) donne le même résultat quels que soient l'adduit et le m/z
//!   (ils ne servent qu'à choisir un fragment) : elle n'est calculée qu'une fois ;
//! - les molécules les plus longues (les plus lentes) sont lancées en premier, pour éviter qu'en fin de calcul
//!   un seul thread travaille pendant que les autres attendent.

use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Once;
use std::time::Duration;

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
        fn process_single(clean_mol: &str, ptype: &str, precursormz: f64, has_mz: bool) -> MolResult;
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

/// Équivalent exact de `rdkit_worker.process_single` : (clé, valeur) dans l'ordre du dictionnaire Python.
pub fn process_single(clean_mol: &str, precursortype: &str, precursormz_str: &str) -> Vec<(&'static str, String)> {
    INIT.call_once(ffi::init_rdkit);
    let mz = py_float(precursormz_str);
    let r = ffi::process_single(clean_mol, py_strip(precursortype), mz.unwrap_or(f64::NAN), mz.is_some());
    let mut out = Vec::with_capacity(6);
    if r.has_inchi { out.push(("INCHI", r.inchi)); }
    if r.has_inchikey { out.push(("INCHIKEY", r.inchikey)); }
    if r.has_smiles { out.push(("SMILES", r.smiles)); }
    if r.has_formula { out.push(("FORMULA", r.formula)); }
    if r.has_exactmass { out.push(("EXACTMASS", py_float_str(r.exactmass))); }
    if r.has_averagemass { out.push(("AVERAGEMASS", py_float_str(r.averagemass))); }
    out
}

/// Le résultat dépend-il de l'adduit et du m/z ? Seulement pour une molécule à plusieurs fragments
/// (« . ») avec un m/z lisible : c'est la seule branche de process_single qui les utilise.
fn depends_on_precursor(clean_mol: &str, precursormz_str: &str) -> bool {
    clean_mol.contains('.') && py_float(precursormz_str).is_some()
}

/// Remplace `rdkit_worker.run_parallel` : calcule toutes les molécules et renvoie le même dictionnaire.
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
    let mut results: Vec<Option<Vec<(&'static str, String)>>> = Vec::new();
    let mut tick_error: Option<E> = None;

    std::thread::scope(|scope| {
        let worker = scope.spawn(|| {
            let computed: Vec<(usize, Vec<(&'static str, String)>)> = order
                .par_iter()
                .filter_map(|&id| {
                    if cancel.load(Ordering::Relaxed) {
                        return None;
                    }
                    let (m, pt, mz) = units[id];
                    let r = process_single(m, pt, mz);
                    done.fetch_add(weights[id], Ordering::Relaxed);
                    Some((id, r))
                })
                .collect();
            let mut by_unit: Vec<Option<Vec<(&'static str, String)>>> = vec![None; units.len()];
            for (id, r) in computed {
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

/// Faut-il utiliser l'ancienne version Python (rdkit_worker.py) ? Oui si FRAGHUB_RDKIT_BACKEND=python.
pub fn use_python_backend() -> bool {
    std::env::var("FRAGHUB_RDKIT_BACKEND").map(|v| v.eq_ignore_ascii_case("python")).unwrap_or(false)
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
    fn arret_demande_pendant_le_calcul() {
        let mols: Vec<_> = (0..2000).map(|i| (format!("C{}O", "C".repeat(i % 40)), String::new(), String::new())).collect();
        let r = compute_cache(mols, |_| Err::<(), &str>("stop"), || std::thread::sleep(Duration::from_millis(1)));
        assert_eq!(r.unwrap_err(), "stop");
    }
}
