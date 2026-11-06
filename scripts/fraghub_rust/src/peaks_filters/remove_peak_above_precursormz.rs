// src/peaks_filters/remove_peak_above_precursormz.rs

/// Supprime les pics dont la valeur m/z est supérieure à la masse du précurseur (+ 5 Da).
///
/// Précurseur multichargé (ex. `[M+2H]2+`) : un fragment qui garde une seule charge peut avoir un m/z
/// bien plus grand que celui du précurseur (jusqu'à environ charge x m/z du précurseur). La limite
/// est donc `charge x m/z du précurseur + 5 Da`. Pour une charge 1, la limite est inchangée.
/// 
/// Pour un développeur Python : On utilise la méthode `.retain()`.
/// C'est l'équivalent ultra-optimisé d'une list comprehension `[p for p in peaks if p[0] < limit]`.
/// En Rust, `.retain()` modifie le vecteur sur place (in-place) sans allouer de nouvelle mémoire.
///
/// # Arguments
/// * `peaks` - Le vecteur `mut`able de pics (la fonction s'approprie la liste et a le droit de la modifier).
/// * `precursormz` - La masse m/z du précurseur (type `f64`, un float 64-bits).
/// * `charge` - Valeur absolue de la charge du précurseur, lue dans l'adduit (1 si inconnue).
/// * `deletion_reason` - Référence mutable pour écrire la raison en cas de suppression totale.
///
/// # Returns
/// * `Vec<(f64, f64)>` : La liste des pics sans ceux au-dessus du seuil.
pub fn remove_peak_above_precursormz(
    mut peaks: Vec<(f64, f64)>,
    precursormz: f64,
    charge: u32,
    deletion_reason: &mut Option<String>
) -> Vec<(f64, f64)> {
    // En Rust, le typage est strict, on écrit 5.0 (float) et non 5 (integer).
    // Pour une charge 1, `precursormz * 1.0` vaut exactement `precursormz` : limite identique à l'ancienne.
    let limit = precursormz * (charge.max(1) as f64) + 5.0;
    
    // `.retain()` garde uniquement les éléments pour lesquels la closure (fonction anonyme) renvoie `true`.
    // Le `&` dans `|&(mz, _)|` sert à déstructurer la référence vers le tuple sans la consommer.
    // L'underscore `_` ignore la seconde valeur (intensité) dont on n'a pas besoin ici.
    peaks.retain(|&(mz, _)| mz < limit);

    if peaks.is_empty() {
        *deletion_reason = Some("spectrum deleted because peaks list is empty after removing peaks above precursor m/z".to_string());
    }

    peaks // Retourne la liste filtrée.
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn charge_un_limite_inchangee() {
        let mut reason = None;
        let out = remove_peak_above_precursormz(vec![(100.0, 1.0), (104.9, 1.0), (105.0, 1.0)], 100.0, 1, &mut reason);
        assert_eq!(out, vec![(100.0, 1.0), (104.9, 1.0)]);
    }

    #[test]
    fn multicharge_garde_les_fragments_au_dessus_du_precurseur() {
        // [M+2H]2+ à m/z 500 : un fragment simplement chargé à m/z 800 est réel, gardé ; 1006 est au-delà de 2 x 500 + 5.
        let mut reason = None;
        let out = remove_peak_above_precursormz(vec![(200.0, 1.0), (800.0, 1.0), (1004.0, 1.0), (1006.0, 1.0)], 500.0, 2, &mut reason);
        assert_eq!(out, vec![(200.0, 1.0), (800.0, 1.0), (1004.0, 1.0)]);
        // Charge 3 : limite 3 x 400 + 5 = 1205.
        let out = remove_peak_above_precursormz(vec![(900.0, 1.0), (1204.0, 1.0), (1210.0, 1.0)], 400.0, 3, &mut reason);
        assert_eq!(out, vec![(900.0, 1.0), (1204.0, 1.0)]);
    }

    #[test]
    fn charge_zero_traitee_comme_un() {
        let mut reason = None;
        let out = remove_peak_above_precursormz(vec![(100.0, 1.0), (300.0, 1.0)], 100.0, 0, &mut reason);
        assert_eq!(out, vec![(100.0, 1.0)]);
    }
}
