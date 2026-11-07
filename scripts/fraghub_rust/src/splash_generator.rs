// src/splash_generator.rs
use pyo3::prelude::*;
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::cmp::Ordering;

const EPS_CORRECTION: f64 = 1.0e-7;
const MZ_PRECISION_FACTOR: f64 = 1_000_000.0;
const INTENSITY_PRECISION_FACTOR: f64 = 1.0;
const INTENSITY_MAP: &[u8] = b"0123456789abcdefghijklmnopqrstuvwxyz";

#[derive(Clone)]
struct SplashSpectrum {
    peaks: Vec<(f64, f64)>, // mz, intensity
}

impl SplashSpectrum {
    fn new(mut raw_peaks: Vec<(f64, f64)>) -> Self {
        let max_int = raw_peaks.iter().map(|&(_, i)| i).fold(0.0_f64, f64::max);
        if max_int > 0.0 {
            for (_, int) in &mut raw_peaks {
                *int = (*int / max_int) * 100.0;
            }
        }
        Self { peaks: raw_peaks }
    }

    fn filter(&self, top_ions: Option<usize>, base_peak_percentage: Option<f64>) -> Self {
        let mut filtered = self.peaks.clone();
        if let Some(pct) = base_peak_percentage {
            let base_int = filtered.iter().map(|&(_, i)| i).fold(0.0_f64, f64::max);
            filtered.retain(|&(_, i)| i + EPS_CORRECTION >= pct * base_int);
        }
        if let Some(top) = top_ions {
            filtered.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal).then(a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal)));
            filtered.truncate(top);
        }
        Self { peaks: filtered }
    }
}

/// Calcule le SPLASH (l'empreinte digitale unique) d'un spectre de masse.
///
/// Pour un développeur Python : Cette fonction effectue un hachage très complexe sur les pics.
/// En Python, de lourdes manipulations de chaînes et des tris (`sort()`) répétés sur chaque spectre
/// finiraient par bloquer le processeur. Ici, tout s'exécute à un niveau extrêmement bas (sur les
/// octets directement via la crate `sha2`) avec un tri très optimisé, ce qui est foudroyant de vitesse.
fn calculate_splash(spectrum: &SplashSpectrum) -> Option<String> {
    if spectrum.peaks.is_empty() { return None; }
    let initial_block = "splash10";

    let calc_hist = |spec: &SplashSpectrum, base: f64, length: usize, bin_size: f64| -> String {
        let mut hist = vec![0.0; length];
        for &(mz, int) in &spec.peaks {
            let idx = (mz / bin_size) as usize % length;
            hist[idx] += int;
        }
        let max_hist = hist.iter().cloned().fold(0.0_f64, f64::max);
        let mut out = String::with_capacity(length);
        for x in hist {
            let val = (EPS_CORRECTION + (base - 1.0) * x / max_hist) as usize;
            out.push(INTENSITY_MAP[val.min(INTENSITY_MAP.len() - 1)] as char);
        }
        out
    };

    let filtered_spec = spectrum.filter(Some(10), Some(0.1));
    let hist_str = calc_hist(&filtered_spec, 3.0, 10, 5.0);

    let mut n = i64::from_str_radix(&hist_str, 3).unwrap_or(0);
    let mut digits = Vec::new();
    if n == 0 { digits.push(0); }
    while n > 0 {
        digits.push((n % 36) as usize);
        n /= 36;
    }
    let mut prefilter_block = String::new();
    for &d in digits.iter().rev() {
        prefilter_block.push(INTENSITY_MAP[d] as char);
    }
    let prefilter_block = format!("{:0>4}", prefilter_block);
    let similarity_block = calc_hist(spectrum, 10.0, 10, 100.0);

    let mut formatted_peaks: Vec<(i64, i64)> = spectrum.peaks.iter().map(|&(mz, int)| {
        (((mz + EPS_CORRECTION) * MZ_PRECISION_FACTOR) as i64, ((int + EPS_CORRECTION) * INTENSITY_PRECISION_FACTOR) as i64)
    }).collect();

    formatted_peaks.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
    let mut spectrum_str = String::new();
    for (i, &(mz, int)) in formatted_peaks.iter().enumerate() {
        if i > 0 { spectrum_str.push(' '); }
        spectrum_str.push_str(&format!("{}:{}", mz, int));
    }

    let mut hasher = Sha256::new();
    hasher.update(spectrum_str.as_bytes());
    let hash_result = hex::encode(hasher.finalize());
    let exact_hash = &hash_result[..20];
    Some(format!("{}-{}-{}-{}", initial_block, prefilter_block, similarity_block, exact_hash))
}

/// Retire les pics inutilisables (intensité ou m/z nuls, négatifs ou non numériques) puis calcule le SPLASH.
///
/// Ces pics étaient retirés plus tard (filtre obligatoire du nettoyage), donc APRÈS le SPLASH : un même spectre
/// avait deux SPLASH selon qu'il était exporté avec ou sans ses pics d'intensité 0 (GNPS : le JSON les garde,
/// le MSP et le MGF non), et le dédoublonnage (SPLASH + InChIKey) ne reconnaissait pas le doublon. Un pic
/// d'intensité 0 n'apporte rien : le spectre gardé ne change pas, seul son SPLASH devient le même partout.
fn clean_peaks_and_splash(spec: &mut crate::spectrum::Spectrum) -> Option<String> {
    let peaks = std::mem::take(&mut spec.peaks);
    spec.peaks = crate::peaks_filters::filters::remove_non_positive_peaks(peaks);
    calculate_splash(&SplashSpectrum::new(spec.peaks.clone()))
}

/// Point d'entrée pour la génération parallèle des identifiants SPLASH.
///
/// Pour un développeur Python : Observez `par_iter()` de Rayon, qui découpe intelligemment 
/// les `chunks` de spectres entre les cœurs CPU. La fonction gère aussi des "Callbacks" (appels) 
/// vers Python/Vue.js pour animer la barre de progression graphique sans bloquer l'interface.
pub fn generate_splash_processing(
    py: Python,
    mut spectrum_list: Vec<crate::spectrum::Spectrum>,
    filename: String,
    progress_callback: Option<PyObject>,
    total_items_callback: Option<PyObject>,
    prefix_callback: Option<PyObject>,
    item_type_callback: Option<PyObject>,
) -> PyResult<Vec<crate::spectrum::Spectrum>> {

    let total = spectrum_list.len();

    // ⚠️ ORDRE CRITIQUE POUR L'INTERFACE VUE.JS
    if let Some(cb) = &total_items_callback { cb.call1(py, (total,))?; }
    if let Some(cb) = &prefix_callback { cb.call1(py, (format!("Generating SPLASH IDs for [{}]:", filename),))?; }
    if let Some(cb) = &item_type_callback { cb.call1(py, ("spectra",))?; }

    let chunk_size = 2000;
    let mut processed = 0;

    for chunk in spectrum_list.chunks_mut(chunk_size) {
        let splash_results: Vec<Option<String>> = chunk
            .par_iter_mut()
            .map(clean_peaks_and_splash)
            .collect();

        for (spec, splash_opt) in chunk.iter_mut().zip(splash_results.into_iter()) {
            if let Some(splash_str) = splash_opt {
                spec.metadata.insert("SPLASH".to_string(), splash_str);
            }
        }

        processed += chunk.len();
        if let Some(cb) = &progress_callback { cb.call1(py, (processed,))?; }
    }

    // ⚠️ GARANTIE DU 100% POUR CLÔTURER L'INTERFACE PROPREMENT
    if let Some(cb) = &progress_callback { cb.call1(py, (total,))?; }

    Ok(spectrum_list)
}


#[cfg(test)]
mod tests_pics_inutilisables {
    use super::*;
    use crate::spectrum::Spectrum;

    fn spec(peaks: &[(f64, f64)]) -> Spectrum {
        let mut s = Spectrum::default();
        s.peaks = peaks.to_vec();
        s
    }

    #[test]
    fn meme_splash_avec_ou_sans_pics_d_intensite_nulle() {
        // GNPS CCMSLIB00000001747 : le JSON a un pic à 0 que le MSP n'a pas
        let msp = [(100.1, 50.0), (150.2, 100.0), (200.3, 20.0)];
        let json = [(100.1, 50.0), (123.488281, 0.0), (150.2, 100.0), (200.3, 20.0)];
        let mut a = spec(&msp);
        let mut b = spec(&json);
        let sa = clean_peaks_and_splash(&mut a);
        let sb = clean_peaks_and_splash(&mut b);
        assert!(sa.is_some());
        assert_eq!(sa, sb);
        assert_eq!(b.peaks, msp.to_vec());
        // Spectre sans pic utilisable : pas de SPLASH, pics vides (supprimé ensuite, motif « peaks list is empty »)
        let mut c = spec(&[(100.0, 0.0), (f64::NAN, 5.0)]);
        assert_eq!(clean_peaks_and_splash(&mut c), None);
        assert!(c.peaks.is_empty());
    }
}
