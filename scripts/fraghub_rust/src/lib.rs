// src/lib.rs
use pyo3::prelude::*;

/// Module principal ("Crate Root") de l'extension Rust pour Python (PyO3).
///
/// Pour un développeur Python : C'est l'équivalent du fichier `__init__.py`. 
/// Il déclare tous les sous-fichiers (`pub mod ...`) pour que le compilateur Rust les inclue.
/// En bas, la macro `#[pymodule]` crée le point d'entrée que Python appellera
/// quand vous ferez `import fraghub_rust`.

// Allocateur mémoire : mimalloc au lieu de celui du système.
// Il ne change aucun résultat, seulement la vitesse des allocations (très nombreuses
// et en parallèle pendant le parsing). Seule la mémoire réservée par Rust est concernée,
// Python garde son propre allocateur.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

// Déclaration de vos modules
pub mod spectrum;
pub mod globals_vars;
pub mod global_state;
pub mod compound_db;
pub mod loading_db;
pub mod loaders;
pub mod convertors;
pub mod splash_generator;
pub mod duplicatas_remover;
pub mod update_checker;
pub mod normalizer;
pub mod peaks_filters;
pub mod spectrum_cleaning;
pub mod complete_from_pubchem_datas;
pub mod ontologies_completion;
pub mod de_novo_calculation;
pub mod normalize_to_not_found;
pub mod splitter;
pub mod spectra_to_msp;
pub mod writers;
pub mod report;
pub mod set_projects;
pub mod deletion_report;
pub mod rdkit_bridge;
pub mod rdkit_native;
pub mod main_orchestrator;
pub mod mz_correction;
pub mod adduct_mass_calculator;

#[pymodule]
fn fraghub_rust(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // --- Orchestrator ---
    m.add_function(wrap_pyfunction!(main_orchestrator::main_orchestrator, m)?)?;
    // Vérification : masse exacte calculée par RDKit C++ (à comparer avec RDKit Python)
    m.add_function(wrap_pyfunction!(rdkit_native::native_exact_mass, m)?)?;

    // --- Partie 1 : Chargement des bases de données ---
    m.add_function(wrap_pyfunction!(loading_db::load_internal_databases, m)?)?;

    // --- Projets et rapports ---
    m.add_function(wrap_pyfunction!(set_projects::reset_updates, m)?)?;
    m.add_function(wrap_pyfunction!(set_projects::init_project, m)?)?;


    Ok(())
}
