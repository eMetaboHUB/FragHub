// src/chemistry_engine.rs
use std::collections::HashMap;
use lazy_static::lazy_static;
use regex::Regex;

// Masses monoisotopiques des éléments (en Daltons)
lazy_static! {
    static ref MONOISOTOPIC: HashMap<&'static str, f64> = {
        let mut m = HashMap::new();
        m.insert("H", 1.00782503207);
        m.insert("C", 12.0);
        m.insert("N", 14.0030740048);
        m.insert("O", 15.99491461956);
        m.insert("P", 30.97376163);
        m.insert("S", 31.97207100);
        m.insert("F", 18.99840322);
        m.insert("Cl", 34.96885268);
        m.insert("Br", 78.91833710);
        m.insert("I", 126.904473);
        m.insert("Si", 27.9769265325);
        m.insert("B", 11.0093054);
        m.insert("Se", 79.9165218);
        m.insert("Li", 7.01600455);
        m.insert("Na", 22.9897692809);
        m.insert("K", 38.96370668);
        m.insert("Ca", 39.96259098);
        m.insert("Mg", 23.98504170);
        m.insert("Fe", 55.9349375);
        m.insert("Zn", 63.92914201);
        m.insert("Cu", 62.92959772);
        m.insert("Mn", 54.9380451);
        m.insert("Co", 58.9331950);
        // Ajoutés le 6 novembre 2026 : oligo-éléments, métaux lourds, médicaments et agents de contraste.
        // Masse de l'isotope le plus abondant, valeurs NIST (« Atomic Weights and Isotopic Compositions »,
        // mêmes valeurs que les 23 éléments ci-dessus, écart < 1e-6 Da), recoupées avec l'évaluation AME2020
        // (écart < 1e-5 Da). Tc est absent : il n'a pas d'isotope stable.
        m.insert("Al", 26.98153863);   // Al-27
        m.insert("V", 50.9439595);     // V-51
        m.insert("Cr", 51.9405075);    // Cr-52
        m.insert("Ni", 57.9353429);    // Ni-58
        m.insert("Ga", 68.9255736);    // Ga-69
        m.insert("Ge", 73.9211778);    // Ge-74
        m.insert("As", 74.9215965);    // As-75
        m.insert("Rb", 84.911789738);  // Rb-85
        m.insert("Sr", 87.9056121);    // Sr-88
        m.insert("Mo", 97.9054082);    // Mo-98
        m.insert("Ag", 106.905097);    // Ag-107
        m.insert("Cd", 113.9033585);   // Cd-114
        m.insert("Sn", 119.9021947);   // Sn-120
        m.insert("Sb", 120.9038157);   // Sb-121
        m.insert("Te", 129.9062244);   // Te-130
        m.insert("Cs", 132.905451933); // Cs-133
        m.insert("Ba", 137.9052472);   // Ba-138
        m.insert("Gd", 157.9241039);   // Gd-158
        m.insert("W", 183.9509312);    // W-184
        m.insert("Pt", 194.9647911);   // Pt-195
        m.insert("Au", 196.9665687);   // Au-197
        m.insert("Hg", 201.970643);    // Hg-202
        m.insert("Tl", 204.9744275);   // Tl-205
        m.insert("Pb", 207.9766521);   // Pb-208
        m.insert("Bi", 208.9803987);   // Bi-209
        m
    };

    static ref ELEM_RE: Regex = Regex::new(r"([A-Z][a-z]?)(\d*)").unwrap();
    static ref CHARGED_RE: Regex = Regex::new(r"^\[([^\]]+)\](\d*)([+-])$").unwrap();
    static ref ADDUCT_RE: Regex = Regex::new(r"^\[(\d*)M((?:[+-][^]]+)*)\](\d*)(\*?[+-])$").unwrap();
    static ref MOD_SPLIT_RE: Regex = Regex::new(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)").unwrap();

    // Nombre de molécules M en tête d'adduit : "[2M+H]+", "[(2M-2H+Ca)-H]-", "[M+H]+" (vide = 1).
    static ref ADDUCT_NMOL_RE: Regex = Regex::new(r"^\s*\[\(?(\d*)M").unwrap();
    // Charge après le crochet fermant : "]2-", "]+", "]*+" (radical), suffixe isotopique toléré ("]2-_37Cl1").
    static ref ADDUCT_CHARGE_RE: Regex = Regex::new(r"\](\d*)\*?[+-]").unwrap();
}

/// Lit le nombre de molécules (n) et la valeur absolue de la charge (z) d'un adduit au format FragHub
/// (colonne `fraghub_default` de `adduct_to_convert.csv`).
///
/// Exemples : "[M+H]+" -> (1, 1) ; "[2M+H]+" -> (2, 1) ; "[M+2H]2+" -> (1, 2) ; "[(2M-2H+Ca)-H]-" -> (2, 1).
/// Retourne `None` si le format n'est pas reconnu.
pub fn adduct_nmol_and_charge(adduct: &str) -> Option<(u32, u32)> {
    let n_caps = ADDUCT_NMOL_RE.captures(adduct)?;
    let z_caps = ADDUCT_CHARGE_RE.captures(adduct)?;
    let n_str = n_caps.get(1).map_or("", |m| m.as_str());
    let z_str = z_caps.get(1).map_or("", |m| m.as_str());
    let n = if n_str.is_empty() { 1 } else { n_str.parse::<u32>().ok()? };
    let z = if z_str.is_empty() { 1 } else { z_str.parse::<u32>().ok()? };
    if n == 0 || z == 0 { return None; }
    Some((n, z))
}

/// m/z du précurseur à partir de la masse exacte neutre M de la molécule et de la différence de masse
/// de l'adduit (`massdiff`, masse TOTALE ajoutée/retirée, électrons compris, telle que dans
/// `adduct_to_convert.csv`) : m/z = (n·M + massdiff) / z.
///
/// Si le format de l'adduit n'est pas reconnu, on retombe sur n = 1, z = 1 (comportement historique).
pub fn precursor_mz_from_exact_mass(exact_mass: f64, mass_diff: f64, adduct: &str) -> f64 {
    let (n, z) = adduct_nmol_and_charge(adduct).unwrap_or((1, 1));
    (n as f64 * exact_mass + mass_diff) / z as f64
}

pub const ELECTRON_MASS: f64 = 0.000548579909;

/// Convertit une formule brute en dictionnaire d'atomes.
pub fn parse_formula(formula: &str) -> HashMap<String, i32> {
    if formula.is_empty() || formula == "nan" || formula == "NOT FOUND" {
        return HashMap::new();
    }
    
    let raw = if let Some(caps) = CHARGED_RE.captures(formula) {
        caps.get(1).map_or(formula, |m| m.as_str())
    } else {
        formula
    };

    let mut counts = HashMap::new();
    for cap in ELEM_RE.captures_iter(raw) {
        let elem = cap[1].to_string();
        let cnt_str = &cap[2];
        let cnt = if cnt_str.is_empty() { 1 } else { cnt_str.parse::<i32>().unwrap_or(1) };
        *counts.entry(elem).or_insert(0) += cnt;
    }
    counts
}

/// Convertit un dictionnaire d'atomes en chaîne de formule brute (convention Hill)
pub fn formula_to_str(counts: &HashMap<String, i32>) -> String {
    if counts.is_empty() {
        return String::new();
    }
    
    let mut parts = Vec::new();
    for elem in ["C", "H"] {
        if let Some(&cnt) = counts.get(elem) {
            if cnt > 0 {
                if cnt > 1 {
                    parts.push(format!("{}{}", elem, cnt));
                } else {
                    parts.push(elem.to_string());
                }
            }
        }
    }
    
    let mut other_elems: Vec<&String> = counts.keys().filter(|&k| k != "C" && k != "H").collect();
    other_elems.sort();
    
    for elem in other_elems {
        if let Some(&cnt) = counts.get(elem) {
            if cnt > 0 {
                if cnt > 1 {
                    parts.push(format!("{}{}", elem, cnt));
                } else {
                    parts.push(elem.to_string());
                }
            }
        }
    }
    
    parts.join("")
}

/// Parse une formule au format crochet et retourne sa charge nette
pub fn formula_net_charge(formula: &str) -> i32 {
    if let Some(caps) = CHARGED_RE.captures(formula) {
        let charge_str = caps.get(2).map_or("", |m| m.as_str());
        let charge = if charge_str.is_empty() { 1 } else { charge_str.parse::<i32>().unwrap_or(1) };
        let sign = caps.get(3).map_or("+", |m| m.as_str());
        if sign == "+" { charge } else { -charge }
    } else {
        0
    }
}

pub fn format_charged_formula(formula_str: &str, charge_abs: i32, sign: &str) -> String {
    if formula_str.is_empty() {
        return String::new();
    }
    let suffix = if charge_abs == 1 {
        sign.to_string()
    } else {
        format!("{}{}", charge_abs, sign)
    };
    format!("[{}]{}", formula_str, suffix)
}

/// Parse un adduct et retourne (n_mol, mods_str, charge_abs, sign)
pub fn parse_precursor_type(pt: &str) -> Option<(i32, String, i32, String)> {
    let pt = pt.trim();
    if let Some(caps) = ADDUCT_RE.captures(pt) {
        let n_mol_str = caps.get(1).map_or("", |m| m.as_str());
        let n_mol = if n_mol_str.is_empty() { 1 } else { n_mol_str.parse::<i32>().unwrap_or(1) };
        
        let mods_str = caps.get(2).map_or("", |m| m.as_str()).to_string();
        
        let charge_str = caps.get(3).map_or("", |m| m.as_str());
        let sign_str = caps.get(4).map_or("+", |m| m.as_str());
        
        let sign_clean = sign_str.replace("*", "");
        let mut charge = if charge_str.is_empty() { 1 } else { charge_str.parse::<i32>().unwrap_or(1) };
        
        let pos = sign_clean.matches('+').count() as i32;
        let neg = sign_clean.matches('-').count() as i32;
        charge *= pos - neg;
        
        let sign = if charge > 0 { "+" } else { "-" }.to_string();
        Some((n_mol, mods_str, charge.abs(), sign))
    } else {
        None
    }
}

pub fn exact_mass_from_counts(counts: &HashMap<String, i32>) -> f64 {
    let mut mass = 0.0;
    for (e, &n) in counts {
        if let Some(&w) = MONOISOTOPIC.get(e.as_str()) {
            mass += w * (n as f64);
        }
    }
    mass
}

pub fn apply_adduct_to_formula(base_formula: &str, precursor_type: &str) -> String {
    if base_formula.is_empty() || precursor_type.is_empty() {
        return String::new();
    }
    
    let mol_charge = formula_net_charge(base_formula);
    
    if precursor_type == "[M]" || precursor_type == "[M]*+" || precursor_type == "[M]+" {
        let total_charge = if mol_charge > 0 { mol_charge } else { 1 };
        let f_str = formula_to_str(&parse_formula(base_formula));
        return format_charged_formula(&f_str, total_charge, "+");
    }
    
    if let Some((n_mol, mods_str, adduct_charge, sign)) = parse_precursor_type(precursor_type) {
        let mut counts = HashMap::new();
        for (elem, &cnt) in parse_formula(base_formula).iter() {
            *counts.entry(elem.clone()).or_insert(0) += cnt * n_mol;
        }
        
        let mut effective_mods = mods_str;
        if mol_charge == -1 {
            effective_mods = format!("+H{}", effective_mods);
        }
        
        for cap in MOD_SPLIT_RE.captures_iter(&effective_mods) {
            let op = &cap[1];
            let coeff_str = &cap[2];
            let coeff = if coeff_str.is_empty() { 1 } else { coeff_str.parse::<i32>().unwrap_or(1) };
            let fragment = &cap[3];
            
            for (elem, &cnt) in parse_formula(fragment).iter() {
                let change = if op == "+" { cnt * coeff } else { -cnt * coeff };
                *counts.entry(elem.clone()).or_insert(0) += change;
            }
        }
        
        counts.retain(|_, &mut v| v > 0);
        let f_str = formula_to_str(&counts);
        
        let h_compensation = if mol_charge == -1 { 1 } else { 0 };
        let total_charge = mol_charge * n_mol + adduct_charge + h_compensation;
        
        if total_charge <= 0 {
            return String::new();
        }
        
        return format_charged_formula(&f_str, total_charge, &sign);
    }
    
    String::new()
}

pub fn compute_mz(precursor_formula: &str) -> Option<f64> {
    if precursor_formula.is_empty() {
        return None;
    }
    
    if let Some(caps) = CHARGED_RE.captures(precursor_formula) {
        let formula_str = caps.get(1).map_or("", |m| m.as_str());
        let charge_str = caps.get(2).map_or("", |m| m.as_str());
        let sign = caps.get(3).map_or("+", |m| m.as_str());
        
        let charge_abs = if charge_str.is_empty() { 1 } else { charge_str.parse::<i32>().unwrap_or(1) };
        if charge_abs <= 0 {
            return None;
        }
        
        let counts = parse_formula(formula_str);
        if counts.is_empty() {
            return None;
        }
        
        let neutral_mass = exact_mass_from_counts(&counts);
        let mz = if sign == "+" {
            (neutral_mass - (charge_abs as f64) * ELECTRON_MASS) / (charge_abs as f64)
        } else {
            (neutral_mass + (charge_abs as f64) * ELECTRON_MASS) / (charge_abs as f64)
        };
        
        Some((mz * 1_000_000.0).round() / 1_000_000.0) // Arrondi à 6 décimales
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn masses_des_elements_ajoutes() {
        // Valeurs de référence calculées avec pyteomics (masses NIST), indépendamment de cette table
        let close = |f: &str, v: f64| (exact_mass_from_counts(&parse_formula(f)) - v).abs() < 1e-6;
        assert!(close("PtCl2H6N2", 298.95559466));      // cisplatine
        assert!(close("C5H11AsO2", 177.99750109));      // arsénobétaïne
        assert!(close("C16H28GdN6O9", 606.11588040));   // gadotéridol
        assert!(close("C2H6Hg", 232.01759319));         // diméthylmercure
    }

    use super::*;

    const CAFFEINE: f64 = 194.080376; // masse monoisotopique neutre

    fn close(a: f64, b: f64) -> bool { (a - b).abs() < 1e-5 }

    #[test]
    fn reads_nmol_and_charge() {
        assert_eq!(adduct_nmol_and_charge("[M+H]+"), Some((1, 1)));
        assert_eq!(adduct_nmol_and_charge("[2M+H]+"), Some((2, 1)));
        assert_eq!(adduct_nmol_and_charge("[M+2H]2+"), Some((1, 2)));
        assert_eq!(adduct_nmol_and_charge("[3M+H]+"), Some((3, 1)));
        assert_eq!(adduct_nmol_and_charge("[2M+2Cl]2-"), Some((2, 2)));
        assert_eq!(adduct_nmol_and_charge("[2M+3Cl]3-_37Cl1"), Some((2, 3)));
        assert_eq!(adduct_nmol_and_charge("[(2M-2H+Ca)-H]-"), Some((2, 1)));
        assert_eq!(adduct_nmol_and_charge("[M]*+"), Some((1, 1)));
        assert_eq!(adduct_nmol_and_charge("M+H"), None);
    }

    #[test]
    fn single_charge_monomer_is_unchanged() {
        // [M+H]+ : 194.080376 + 1.007276
        assert!(close(precursor_mz_from_exact_mass(CAFFEINE, 1.007276, "[M+H]+"), 195.087652));
    }

    #[test]
    fn dimer_uses_twice_the_mass() {
        // [2M+H]+ : 2 x 194.080376 + 1.007276
        assert!(close(precursor_mz_from_exact_mass(CAFFEINE, 1.007276, "[2M+H]+"), 389.168028));
    }

    #[test]
    fn multiply_charged_is_divided_by_z() {
        // [M+2H]2+ : (194.080376 + 2.014552) / 2
        assert!(close(precursor_mz_from_exact_mass(CAFFEINE, 2.014552, "[M+2H]2+"), 98.047464));
        // [2M+2Cl]2- : (2 x 194.080376 + 69.938802) / 2
        assert!(close(precursor_mz_from_exact_mass(CAFFEINE, 69.938802, "[2M+2Cl]2-"), 229.049777));
    }

    #[test]
    fn unknown_format_falls_back_to_old_behaviour() {
        assert!(close(precursor_mz_from_exact_mass(CAFFEINE, 1.007276, "M+H"), 195.087652));
    }
}
