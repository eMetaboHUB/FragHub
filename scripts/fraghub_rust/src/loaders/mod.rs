// src/loaders/mod.rs
//! Lecture des fichiers d'entrée, un fichier par format : chaque lecteur découpe le fichier en
//! spectres bruts (textes), ensuite analysés par les parseurs (module parsers : msp_parser.rs, mgf_parser.rs, json_parser.rs).
//! - `json_loader.rs`  : JSON en tableau (format GNPS) ;
//! - `jsonl_loader.rs` : JSONL (un spectre JSON par ligne) ;
//! - `msp_loader.rs`   : MSP ;
//! - `mgf_loader.rs`   : MGF.
//! (Le CSV est lu et analysé en une seule étape, dans parsers/csv_parser.rs.)

use pyo3::prelude::*;
use sha2::{Sha256, Digest};

mod json_loader;
mod jsonl_loader;
mod msp_loader;
mod mgf_loader;

pub use json_loader::load_spectrum_list_json;
pub use jsonl_loader::load_spectrum_list_json_2;
pub use msp_loader::load_spectrum_list_from_msp;
pub use mgf_loader::load_spectrum_list_from_mgf;

/// Clés ajoutées par les lecteurs JSON / JSONL en tête de chaque spectre : nom et empreinte du fichier.
/// Préfixe à part pour ne jamais se confondre avec une clé du fichier (avant : "filename" et "filehash",
/// remplacés par une clé "filename" présente dans le spectre).
pub const INJECTED_PREFIX: &str = "__fraghub_";
pub const INJECTED_FILENAME: &str = "__fraghub_filename";
pub const INJECTED_FILEHASH: &str = "__fraghub_filehash";

/// Gets the size of a file in bytes, converts it to a string, and returns the SHA-256 hash.
///
/// Pour un développeur Python : La macro `#[pyfunction]` indique à PyO3 que cette fonction
/// Rust doit être compilée et exposée comme une fonction appelable directement depuis Python.
///
/// # Arguments
/// * `file_path` (&str) : Le chemin du fichier.
///
/// # Returns
/// * `String` : Le hash SHA-256 de la taille du fichier.
#[pyfunction]
pub fn generate_file_hash(file_path: &str) -> String {
    if let Ok(metadata) = std::fs::metadata(file_path) {
        let size = metadata.len().to_string();
        let mut hasher = Sha256::new();
        hasher.update(size.as_bytes());
        format!("{:x}", hasher.finalize())
    } else {
        format!("Error: File not found at {}", file_path)
    }
}
