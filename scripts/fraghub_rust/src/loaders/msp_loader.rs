// src/loaders/msp_loader.rs
//! Lecture des fichiers MSP.
use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::{Duration, Instant};
use regex::Regex;
use once_cell::sync::Lazy;

use super::generate_file_hash;

/// Ligne de pic : deux nombres séparés par des espaces, une tabulation ou « : » (suivis éventuellement d'une
/// annotation ou d'autres pics au format NIST).
static PEAK_LINE: Lazy<Regex> = Lazy::new(|| Regex::new(r"^-?\d+(?:[.,]\d+)?(?:[Ee][+-]?\d+)?(?:[ \t]+|:)-?\d+").unwrap());

/// Regroupe les blocs de lignes (séparés par des lignes vides) en spectres.
///
/// Un bloc qui ne contient QUE des lignes de pics est la suite du spectre précédent : certains fichiers
/// (ex. MoNA, MoNA_0015754) mettent une ligne vide entre les pics. Avant, le spectre s'arrêtait à la première
/// ligne vide (il perdait ses autres pics) et chaque pic suivant devenait un faux spectre, rejeté. Un fichier
/// sans ce défaut donne exactement les mêmes spectres qu'avant.
#[derive(Default)]
struct BlockJoiner {
    pending: String,
}

impl BlockJoiner {
    /// Reçoit le bloc qui vient de finir (et le vide) ; renvoie le spectre précédent s'il est terminé.
    fn push(&mut self, body: &mut String) -> Option<String> {
        if body.is_empty() {
            return None;
        }
        if !self.pending.is_empty() && body.split('\n').all(|l| PEAK_LINE.is_match(l)) {
            self.pending.push('\n');
            self.pending.push_str(body);
            body.clear();
            return None;
        }
        let done = (!self.pending.is_empty()).then(|| std::mem::take(&mut self.pending));
        std::mem::swap(&mut self.pending, body);
        body.clear();
        done
    }

    /// Dernier spectre du fichier.
    fn finish(self) -> Option<String> {
        (!self.pending.is_empty()).then_some(self.pending)
    }
}

// ----------------------------------------------------------------------
// MSP FORMAT
// ----------------------------------------------------------------------
/// Assemble le texte final d'un spectre MSP : en-tête FILENAME/FILEHASH + lignes du spectre.
///
/// Ancienne méthode : on écrivait "FILENAME: <nom>\n<lignes>\n", puis une regex remplaçait cette
/// première ligne par `replacement`, puis `trim_end()`. Comme la première ligne est toujours celle
/// qu'on vient d'écrire, on peut écrire `replacement` directement : même texte, sans parcourir
/// tout le spectre ni le recopier trois fois.
/// Exception : si le nom contient '$' (interprété par Regex::replace comme un groupe de capture)
/// ou un retour à la ligne, on garde l'ancienne méthode pour obtenir exactement le même texte.
fn assemble_msp_spectrum(filename: &str, replacement: &str, re: &Regex, body: &str) -> String {
    if filename.contains('$') || filename.contains('\n') {
        let spectrum = format!("FILENAME: {}\n{}\n", filename, body);
        return re.replace(&spectrum, replacement).trim_end().to_string();
    }
    let mut s = String::with_capacity(replacement.len() + body.len() + 1);
    s.push_str(replacement);
    s.push_str(body);
    s.push('\n');
    let len = s.trim_end().len();
    s.truncate(len);
    s
}

/// RUST NATIVE MSP LOADER
///
/// # Arguments
/// * `py` (Python) : Token PyO3.
/// * `msp_file_path` (&str) : Le chemin vers le fichier MSP.
/// * `progress_callback` (Option<PyObject>) : Callback pour la progression en MB.
/// * `total_items_callback` (Option<PyObject>) : Callback pour la taille totale en MB.
/// * `prefix_callback` (Option<PyObject>) : Callback du préfixe.
/// * `item_type_callback` (Option<PyObject>) : Callback du type d'élément ("MB").
///
/// # Returns
/// * `PyResult<Vec<String>>` : Une liste de spectres bruts sous forme de chaînes de caractères.
pub fn load_spectrum_list_from_msp(
    py: Python,
    msp_file_path: &str,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<String>> {
    let file_hash = generate_file_hash(msp_file_path);
    let path = Path::new(msp_file_path);
    let filename = path.file_name().unwrap_or_default().to_string_lossy().to_string();

    let total_bytes = std::fs::metadata(msp_file_path).map(|m| m.len()).unwrap_or(0);
    let total_mb = std::cmp::max(1, total_bytes / (1024 * 1024));

    if let Some(cb) = &total_items_callback { let _ = cb.call1(py, (total_mb, 0)); }
    if let Some(cb) = &prefix_callback { let _ = cb.call1(py, (format!("loading [{}]:", filename),)); }
    if let Some(cb) = &item_type_callback { let _ = cb.call1(py, ("MB",)); }

    let file = File::open(msp_file_path)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 4, file);
    let mut spectrum_list = Vec::new();
    // Lignes du spectre en cours, déjà séparées par "\n" (avant : une String par ligne puis join)
    let mut body = String::with_capacity(16 * 1024);
    let mut total_consumed: u64 = 0;
    let mut last_update = Instant::now();

    let re = Regex::new(r"(?i)FILENAME: .*\n").unwrap();
    let replacement = format!("FILENAME: {}\nFILEHASH: {}\n", filename, file_hash);

    let mut line = String::with_capacity(1024);
    let mut joiner = BlockJoiner::default();
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 { break; }
        total_consumed += bytes_read as u64;

        let trimmed = line.trim();
        if trimmed.is_empty() {
            if !body.is_empty() {
                if let Some(done) = joiner.push(&mut body) {
                    spectrum_list.push(assemble_msp_spectrum(&filename, &replacement, &re, &done));
                }

                if last_update.elapsed() >= Duration::from_millis(50) {
                    let current_mb = total_consumed / (1024 * 1024);
                    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (current_mb,)); }
                    py.allow_threads(|| std::thread::sleep(Duration::from_millis(1)));
                    last_update = Instant::now();
                }
            }
        } else {
            if !body.is_empty() {
                body.push('\n');
            }
            body.push_str(trimmed);
        }
    }

    if let Some(done) = joiner.push(&mut body) {
        spectrum_list.push(assemble_msp_spectrum(&filename, &replacement, &re, &done));
    }
    if let Some(done) = joiner.finish() {
        spectrum_list.push(assemble_msp_spectrum(&filename, &replacement, &re, &done));
    }

    if let Some(cb) = &progress_callback { let _ = cb.call1(py, (total_mb,)); }

    Ok(spectrum_list)
}

#[cfg(test)]
mod tests {
    use super::assemble_msp_spectrum;
    use super::BlockJoiner;
    use regex::Regex;

    /// Découpe un texte MSP comme le lecteur (lignes vides = fin de bloc), avec regroupement des blocs.
    fn spectra(text: &str) -> Vec<String> {
        let mut out = Vec::new();
        let mut j = BlockJoiner::default();
        let mut body = String::new();
        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() {
                if let Some(d) = j.push(&mut body) { out.push(d); }
            } else {
                if !body.is_empty() { body.push('\n'); }
                body.push_str(t);
            }
        }
        if let Some(d) = j.push(&mut body) { out.push(d); }
        if let Some(d) = j.finish() { out.push(d); }
        out
    }

    #[test]
    fn pics_separes_par_des_lignes_vides_rattaches_au_spectre() {
        // Cas réel du MSP MoNA (MoNA_0015754) : une ligne vide entre chaque pic
        let t = "Name: A\nDB#: MoNA_0015754\nNum Peaks: 3\n147.0706 1.409837\n\n148.0828 3.317018\n\n172.0348 1.410634\n\nName: B\nNum Peaks: 1\n100 1\n";
        assert_eq!(spectra(t), vec![
            "Name: A\nDB#: MoNA_0015754\nNum Peaks: 3\n147.0706 1.409837\n148.0828 3.317018\n172.0348 1.410634".to_string(),
            "Name: B\nNum Peaks: 1\n100 1".to_string(),
        ]);
        // Métadonnées puis pics après une ligne vide
        assert_eq!(spectra("Name: C\nNum Peaks: 2\n\n85 999; 86 23;\n87\t5 \"ann\"\n"),
                   vec!["Name: C\nNum Peaks: 2\n85 999; 86 23;\n87\t5 \"ann\"".to_string()]);
    }

    #[test]
    fn fichier_sans_ce_defaut_inchange() {
        let t = "Name: A\nNum Peaks: 1\n100 1\n\n\nName: B\nComments: \"x=1\"\nNum Peaks: 2\n100 1\n101 2\n\nName: C\n50 2\n";
        assert_eq!(spectra(t), vec!["Name: A\nNum Peaks: 1\n100 1".to_string(),
            "Name: B\nComments: \"x=1\"\nNum Peaks: 2\n100 1\n101 2".to_string(), "Name: C\n50 2".to_string()]);
        // Un bloc de pics en tête de fichier (rien avant) reste un bloc à part
        assert_eq!(spectra("100 1\n\nName: A\n100 1\n"), vec!["100 1".to_string(), "Name: A\n100 1".to_string()]);
    }

    /// Copie exacte de l'ancien assemblage (buffer de lignes + join + regex + trim_end).
    fn ancien(filename: &str, file_hash: &str, lines: &[&str]) -> String {
        let re = Regex::new(r"(?i)FILENAME: .*\n").unwrap();
        let replacement = format!("FILENAME: {}\nFILEHASH: {}\n", filename, file_hash);
        let mut buffer = vec![format!("FILENAME: {}", filename)];
        for l in lines {
            buffer.push(l.to_string());
        }
        let spectrum = buffer.join("\n") + "\n";
        let replaced = re.replace(&spectrum, replacement.as_str()).to_string();
        replaced.trim_end().to_string()
    }

    #[test]
    fn assemblage_identique_a_l_ancien() {
        let re = Regex::new(r"(?i)FILENAME: .*\n").unwrap();
        let noms = ["MoNA-export-LC-MS-MS_Spectra.msp", "", "nom avec espaces.msp", "accentué_é.msp",
                    "prix$1.msp", "a$$b.msp", "${0}.msp", "ligne\nnouvelle.msp", "filename: piège.msp"];
        let spectres: [&[&str]; 4] = [
            &["Name: Caffeine", "Num Peaks: 2", "138.0662 100", "195.0877 45.2"],
            &["seule ligne"],
            &["FILENAME: autre.msp", "Comments: \"filename: x\"", "1 2"],
            &["Name: $1 et \\n", "Num Peaks: 1", "10,5 20,25"],
        ];
        for nom in noms {
            let hash = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
            let replacement = format!("FILENAME: {}\nFILEHASH: {}\n", nom, hash);
            for lignes in spectres {
                let body = lignes.join("\n");
                assert_eq!(assemble_msp_spectrum(nom, &replacement, &re, &body), ancien(nom, hash, lignes),
                           "nom {nom:?}, lignes {lignes:?}");
            }
        }
    }
}
