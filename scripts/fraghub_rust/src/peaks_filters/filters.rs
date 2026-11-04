// src/peaks_filters/filters.rs
use std::collections::HashMap;

/// Filtre obligatoire pour enlever les pics d'intensité négative ou nulle.
///
/// Pour un développeur Python : En Rust, passer `mut peaks: Vec` transfère la "propriété" (Ownership)
/// de la variable à cette fonction. La fonction a le droit de la modifier et doit la renvoyer.
///
/// # Arguments
/// * `peaks` (Vec<(f64, f64)>) : La liste des tuples `(m/z, intensité)` à filtrer.
///
/// # Returns
/// * `Vec<(f64, f64)>` : La liste de pics nettoyée.
pub fn remove_non_positive_peaks(mut peaks: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    peaks.retain(|&(_, int)| int > 0.0);
    peaks
}

/// Point d'entrée principal pour filtrer un spectre selon les préférences de l'utilisateur.
///
/// Pour un développeur Python : On reçoit un dictionnaire `parameters`.
/// Le type `&HashMap<String, f64>` indique un emprunt "en lecture seule" (Borrowing).
/// On ne peut pas modifier ce dictionnaire, mais on peut le lire librement.
///
/// L'utilisation de `Option<f64>` pour `precursormz` est le pendant Rust de `Optional[float]` 
/// (peut être un float ou None).
///
/// # Arguments
/// * `peaks` (Vec<(f64, f64)>) : La liste de pics du spectre à traiter.
/// * `precursormz` (Option<f64>) : La masse m/z du précurseur (optionnelle).
/// * `parameters` (&HashMap<String, f64>) : Le dictionnaire des paramètres définis par l'utilisateur.
/// * `deletion_reason` (&mut Option<String>) : Pointeur mutable pour inscrire la raison si le spectre est rejeté.
///
/// # Returns
/// * `Vec<(f64, f64)>` : La nouvelle liste de pics après tous les filtres (peut être vide si le spectre est rejeté).
pub fn apply_filters(
    mut peaks: Vec<(f64, f64)>,
    precursormz: Option<f64>,
    parameters: &HashMap<String, f64>,
    deletion_reason: &mut Option<String>
) -> Vec<(f64, f64)> {

    // Ordre des filtres :
    //   1. ceux qui RETIRENT des pics (m/z > précurseur, hors plage m/z, réduction au top N),
    //   2. puis ceux qui VÉRIFIENT le spectre final (nombre minimal de pics, pics intenses),
    //   3. la normalisation d'intensité en dernier parmi les transformations,
    // afin que les seuils et le maximum (= 1) s'appliquent au spectre réellement écrit.

    // --- Step 1: Mandatory filter ---
    peaks = remove_non_positive_peaks(peaks);

    // --- Filter 1: Remove peaks above precursor m/z ---
    if parameters.get("remove_peak_above_precursormz").copied().unwrap_or(0.0) == 1.0 {
        // `if let Some(pmz) = precursormz` vérifie que precursormz n'est pas None,
        // et déballe sa valeur dans `pmz` si elle existe. (Pattern matching)
        if let Some(pmz) = precursormz {
            peaks = super::remove_peak_above_precursormz::remove_peak_above_precursormz(peaks, pmz, deletion_reason);
            if peaks.is_empty() { return peaks; }
        }
    }

    // --- Filter 2: Keep peaks within a user-defined m/z range ---
    if parameters.get("keep_mz_in_range").copied().unwrap_or(0.0) == 1.0 {
        let mz_from = parameters.get("keep_mz_in_range_from_mz").copied().unwrap_or(0.0);
        let mz_to = parameters.get("keep_mz_in_range_to_mz").copied().unwrap_or(0.0);
        peaks = super::keep_mz_in_range::keep_mz_in_range(peaks, mz_from, mz_to, deletion_reason);
        if peaks.is_empty() { return peaks; }
    }

    // --- Filter 3: Reduce peak list to a maximum number of peaks ---
    if parameters.get("reduce_peak_list").copied().unwrap_or(0.0) == 1.0 {
        let max_peaks = parameters.get("reduce_peak_list_max_peaks").copied().unwrap_or(0.0) as usize;
        peaks = super::reduce_peak_list::reduce_peak_list(peaks, max_peaks);
    }

    // --- Filter 4: Check minimum required peak count (sur la liste finale) ---
    // En Rust, pour lire une valeur d'un HashMap, on utilise `.get()`. Cela retourne un `Option<&f64>`.
    // `.copied()` transforme le `Option<&f64>` en `Option<f64>`.
    // `.unwrap_or(0.0)` extrait la valeur, ou retourne `0.0` si la clé n'existe pas. (Remplace le `.get(key, 0.0)` en Python).
    if parameters.get("check_minimum_peak_requiered").copied().unwrap_or(0.0) == 1.0 {
        // `as usize` effectue un cast de `f64` en `usize` (entier positif).
        let n_peaks = parameters.get("check_minimum_peak_requiered_n_peaks").copied().unwrap_or(0.0) as usize;
        peaks = super::check_minimum_peak_requiered::check_minimum_peak_requiered(peaks, n_peaks, deletion_reason);
        if peaks.is_empty() { return peaks; }
    }

    // --- Filter 5: Normalize peak intensity (après toute suppression de pics, le max vaut donc 1) ---
    if parameters.get("normalize_intensity").copied().unwrap_or(0.0) == 1.0 {
        peaks = super::normalize_intensity::normalize_intensity(peaks);
        if peaks.is_empty() { return peaks; }
    }

    // --- Filter 6: Check minimum number of high-intensity peaks ---
    if parameters.get("check_minimum_of_high_peaks_requiered").copied().unwrap_or(0.0) == 1.0 {
        let intensity_percent = parameters.get("check_minimum_of_high_peaks_requiered_intensity_percent").copied().unwrap_or(0.0);
        let no_peaks = parameters.get("check_minimum_of_high_peaks_requiered_no_peaks").copied().unwrap_or(0.0) as usize;
        peaks = super::check_minimum_of_high_peaks_requiered::check_minimum_of_high_peaks_requiered(peaks, intensity_percent, no_peaks, deletion_reason);
        if peaks.is_empty() { return peaks; }
    }

    peaks
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(pairs: &[(&str, f64)]) -> HashMap<String, f64> {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    #[test]
    fn min_peak_count_is_checked_after_peak_removal() {
        // 4 pics, dont 2 au-dessus du précurseur (100 + 5 Da) : il en reste 2 < 3 requis -> spectre rejeté.
        let peaks = vec![(50.0, 10.0), (80.0, 20.0), (150.0, 30.0), (200.0, 40.0)];
        let p = params(&[
            ("check_minimum_peak_requiered", 1.0), ("check_minimum_peak_requiered_n_peaks", 3.0),
            ("remove_peak_above_precursormz", 1.0),
        ]);
        let mut reason = None;
        let out = apply_filters(peaks, Some(100.0), &p, &mut reason);
        assert!(out.is_empty());
        assert!(reason.unwrap().contains("number of peaks is below"));
    }

    #[test]
    fn min_peak_count_is_checked_after_mz_range() {
        let peaks = vec![(30.0, 10.0), (40.0, 20.0), (60.0, 30.0), (70.0, 40.0)];
        let p = params(&[
            ("check_minimum_peak_requiered", 1.0), ("check_minimum_peak_requiered_n_peaks", 3.0),
            ("keep_mz_in_range", 1.0), ("keep_mz_in_range_from_mz", 50.0), ("keep_mz_in_range_to_mz", 2000.0),
        ]);
        let mut reason = None;
        assert!(apply_filters(peaks, None, &p, &mut reason).is_empty());
    }

    #[test]
    fn max_intensity_is_one_after_all_filters() {
        // Le pic le plus intense (m/z 300) est hors plage : après filtrage, le maximum doit quand même valoir 1.
        let peaks = vec![(60.0, 10.0), (70.0, 40.0), (300.0, 100.0)];
        let p = params(&[
            ("normalize_intensity", 1.0),
            ("keep_mz_in_range", 1.0), ("keep_mz_in_range_from_mz", 50.0), ("keep_mz_in_range_to_mz", 200.0),
        ]);
        let mut reason = None;
        let out = apply_filters(peaks, None, &p, &mut reason);
        assert_eq!(out, vec![(60.0, 0.25), (70.0, 1.0)]);
    }

    #[test]
    fn valid_spectrum_passes_all_filters() {
        let peaks = vec![(55.0, 5.0), (60.0, 50.0), (75.0, 100.0), (90.0, 80.0)];
        let p = params(&[
            ("normalize_intensity", 1.0), ("remove_peak_above_precursormz", 1.0),
            ("check_minimum_peak_requiered", 1.0), ("check_minimum_peak_requiered_n_peaks", 3.0),
            ("reduce_peak_list", 1.0), ("reduce_peak_list_max_peaks", 500.0),
            ("keep_mz_in_range", 1.0), ("keep_mz_in_range_from_mz", 50.0), ("keep_mz_in_range_to_mz", 2000.0),
            ("check_minimum_of_high_peaks_requiered", 1.0),
            ("check_minimum_of_high_peaks_requiered_intensity_percent", 5.0),
            ("check_minimum_of_high_peaks_requiered_no_peaks", 2.0),
        ]);
        let mut reason = None;
        let out = apply_filters(peaks, Some(120.0), &p, &mut reason);
        assert_eq!(out.len(), 4);
        assert!(reason.is_none());
        assert_eq!(out.iter().map(|p| p.1).fold(0.0_f64, f64::max), 1.0);
    }
}
