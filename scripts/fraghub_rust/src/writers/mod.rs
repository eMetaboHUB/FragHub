// src/writers/mod.rs
//! Écriture des fichiers de sortie de FragHub, un fichier par format :
//! - `csv_writer.rs`  : fichiers CSV (CSV/POS, CSV/NEG) ;
//! - `json_writer.rs` : fichiers JSON (JSON/POS, JSON/NEG) ;
//! - `msp_writer.rs`  : fichiers MSP (MSP/POS, MSP/NEG), à partir des textes préparés par spectra_to_msp.rs.
//!
//! (Le format mzSpecLib est écrit ensuite en Python, par scripts/mzspeclib_converter.py, à partir des MSP.)

mod csv_writer;
mod json_writer;
mod msp_writer;

pub use csv_writer::writting_csv_processing;
pub use json_writer::writting_json_processing;
pub use msp_writer::writting_msp_processing;

// ---------------------------------------------------------------------------------------------
// Écriture d'un fichier CSV ou JSON : mise en forme EN PARALLÈLE, écriture DANS L'ORDRE.
//
// Avant, chaque fichier était écrit par un seul thread (le parallélisme se faisait fichier par
// fichier) : POS_LC, qui contient l'essentiel des spectres, occupait un seul cœur pendant que les
// autres attendaient. Maintenant, les spectres d'un fichier sont découpés en paquets ; chaque paquet
// est mis en forme par un thread dans un tampon mémoire, puis les tampons sont écrits dans l'ordre.
// La mise en forme d'une ligne est EXACTEMENT celle d'avant (mêmes fonctions) : fichiers identiques.
// ---------------------------------------------------------------------------------------------

/// Nombre de spectres mis en forme par paquet (un paquet = une tâche pour un thread).
pub(super) const WRITE_CHUNK: usize = 2000;
/// Nombre de paquets mis en forme avant d'écrire sur le disque (limite la mémoire utilisée).
pub(super) const WRITE_CHUNKS_PER_ROUND: usize = 64;

/// Jeu de spectres de test commun aux tests des différents formats.
#[cfg(test)]
pub(crate) mod test_data {
    use crate::spectrum::Spectrum;
    use std::collections::HashMap;
    use std::fs;

    pub(crate) fn jeu(n: usize, graine: u64) -> (Vec<String>, Vec<Spectrum>) {
        let cols: Vec<String> = ["NAME", "PRECURSORMZ", "MSLEVEL", "RT", "COMMENT", "SMILES", "NUM PEAKS", "PEAKS_LIST", "TREE_X"].iter().map(|s| s.to_string()).collect();
        let mut x = graine;
        let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let mut specs = Vec::new();
        for i in 0..n {
            let mut m = HashMap::new();
            m.insert("NAME".into(), match i % 5 { 0 => format!("nom;avec;points-virgules {i}"), 1 => format!("nom \"guillemets\" {i}"), 2 => format!("ligne\nsuivante {i}"), 3 => "nan".into(), _ => format!("Composé é {i}") });
            m.insert("PRECURSORMZ".into(), format!("{}", (rnd() % 100000) as f64 / 97.0));
            m.insert("MSLEVEL".into(), if i % 7 == 0 { "MS2".into() } else { "2".into() });
            if i % 3 != 0 { m.insert("RT".into(), format!("{:.3}", (rnd() % 1000) as f64 / 7.0)); }
            m.insert("COMMENT".into(), if i % 4 == 0 { String::new() } else { format!("c{i}") });
            m.insert("SMILES".into(), "CC(=O)O".into());
            let np = (rnd() % 12) as usize;
            let peaks: Vec<(f64, f64)> = (0..np).map(|k| (50.0 + k as f64 * 13.371 + (rnd() % 1000) as f64 / 1e4, (rnd() % 1000) as f64 / 999.0)).collect();
            m.insert("NUM PEAKS".into(), np.to_string());
            if i % 6 == 0 && np > 0 {
                // pics annotés de novo (format avec formules)
                m.insert("PEAKS_LIST".into(), peaks.iter().map(|p| format!("{} {} C2H3O+", p.0, p.1)).collect::<Vec<_>>().join("\n"));
            } else if i % 6 == 1 {
                m.insert("PEAKS_LIST".into(), peaks.iter().map(|p| format!("{:.8} {:.8}", p.0, p.1)).collect::<Vec<_>>().join("\n"));
            }
            specs.push(Spectrum { metadata: m, peaks });
        }
        (cols, specs)
    }

    pub(crate) fn dossier(nom: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("fraghub_writers_{}_{}", nom, std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }
}
