// src/globals_vars.rs

use once_cell::sync::Lazy;
use regex::Regex;

/// Expressions régulières compilées une seule fois pour tout le programme.
///
/// Pour un développeur Python : Si vous utilisez `re.match()` ou compilez une regex à l'intérieur
/// d'une boucle Python, cela ralentit tout. Ici, on utilise `Lazy` pour compiler chaque Regex 
/// au tout premier appel, puis on les garde en RAM. C'est l'un des secrets des 
/// performances exceptionnelles de ce parseur natif.

// =================================================== REGEX PATTERN ====================================================

// ============ Parsors regex ============

pub static IS_ADDUCT_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)m\]?(\-|\+)").unwrap());

pub static METADATA_STRIP_VALUE_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r#"^"|"$"#).unwrap());

pub static METADATA_FIELDS_NAME_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"^[\W_]+|[\W_]+$").unwrap());

pub static METADATA_PATTERN_MGF: Lazy<Regex> = Lazy::new(|| Regex::new(r"([^:\n]*?)=\s*([^\n]*)(?:\n|$)").unwrap());

pub static METADATA_PATTERN_MSP: Lazy<Regex> = Lazy::new(|| Regex::new(r"([^:]*):(?: )?([^\n]*)(?:\n|$)").unwrap());

pub static COMPUTED_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)computed").unwrap());

pub static COMMENT_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)comment.*").unwrap());

pub static PEAK_LIST_SPLIT_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?:^|\n)(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)").unwrap());

pub static PEAK_LIST_JSON_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:|,|, )(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)").unwrap());

pub static SUB_FIELDS_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(\S+?)="([^"]*)"|"(\w+?)=([^"]*)"|"([^"]*?)=([^"]*)"|(\S+?)=(\d+(?:[.,]\d*)?)|(\S+?)=(.*?)(?:;|\n|$)"#).unwrap());

// Version corrigée (Greedy, identique à Python) : Pour le MGF (Même chose, on remet [\s\S]* et le =)
pub static METADATA_PEAK_LIST_SPLIT_PATTERN_MGF: Lazy<Regex> = Lazy::new(|| Regex::new(r"([\s\S]*=.*[0-9]*\n)((?:(?:-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(?:-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?).*(?:\n|$))*)").unwrap());

// Version corrigée (Greedy, identique à Python) : Pour le MSP (On remet [\s\S]* au début et on supprime (?s))
pub static METADATA_PEAK_LIST_SPLIT_PATTERN_MSP: Lazy<Regex> = Lazy::new(|| Regex::new(r"([\s\S]*:.*[0-9]*\n)((?:(?:-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(?:-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?).*(?:\n|$))*)").unwrap());

// ======================================

// ===== normalizers regex pattern ======

pub static INDIGO_SMILES_CORRECTION_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"\|[\s\S]*").unwrap());

pub static SUB_SIGNE_END_ADDUCT_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(\-|\+)$").unwrap());

pub static SUB_ADDUCT_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"\(|\)|(.*\[)|(\]([\d\+\-\*]*)?)").unwrap());

pub static FLOAT_CHECK_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)").unwrap());

pub static MS_LEVEL_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)(?:ms)?(\d)").unwrap());

pub static IONMODE_POS_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)^p|^\+|^pos").unwrap());

pub static IONMODE_NEG_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)^n|^\-|^neg").unwrap());

// Préfixe "InChI=" éventuel (y compris doublé, ex. "InChI=InChI=1S/...", vu dans MoNA), remplacé par un seul "InChI=".
pub static REPAIR_INCHI_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)^\s*(?:inchi=)*").unwrap());

// Match inchi : la valeur doit COMMENCER par "InChI=1S/..." ou "1S/..." (préfixe absent).
// Ancré au début : un SMILES ne peut pas commencer par un chiffre, et un SMILES stéréo
// contenant "/C/" n'est plus confondu avec un InChI.
pub static INCHI_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)^\s*(?:InChI=)*1S?/\S+").unwrap());

// Match smiles
pub static SMILES_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)[^J][a-z0-9@+\-\[\]\(\)\\/%=#$]{6,}").unwrap());

// Match inchikey complet (XXXXXXXXXXXXXX-YYYYYYYYYY-Z) ou court (14 majuscules), valeur ENTIÈRE uniquement.
// Ancré des deux côtés : un SMILES contenant 14 lettres consécutives (chaîne alkyle C14+) n'est plus reconnu.
// À n'utiliser que sur le champ INCHIKEY lui-même (une forme courte seule reste ambiguë, ex. tétradécane "CCCCCCCCCCCCCC").
pub static INCHIKEY_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"^\s*(?:(?i:[A-Z]{14}-[A-Z]{10}-[A-Z])|[A-Z]{14})\s*$").unwrap());

// Match inchikey complet uniquement. Sert à détecter un InChIKey rangé dans un AUTRE champ (SMILES, INCHI),
// où la forme courte serait indiscernable d'un SMILES purement carboné.
pub static INCHIKEY_FULL_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)^\s*[A-Z]{14}-[A-Z]{10}-[A-Z]\s*$").unwrap());

// Termes qui signalent un spectre prédit (in silico). Testés sur COMMENT (après retrait des formulations
// de IN_SILICO_EXCLUDE_PATTERN), et sur FILENAME + NAME.
// (Ancien terme "Annotation\.level\.3" retiré : il ne reconnaissait rien, et chez MassBank "Annotation level-3"
//  désigne la confiance de l'identification d'un spectre expérimental.)
//
// Outils qui GÉNÈRENT des spectres (CFM-ID, QCEIMS, NEIMS, MassFormer, ICEBERG, SCARF, 3DMolMS, GRAFF-MS, FraGNNet) et
// bibliothèques in silico (ISDB, LipidBlast) : reconnus en mots entiers ("_" et "-" comptent comme séparateurs, pour
// "UNPD_ISDB_R_p01.mgf"), sinon "ISDB" se trouverait au milieu d'InChIKeys (ex. "...-GOSISDBHSA-N").
// ICEBERG et SCARF : en majuscules seulement (sinon "iceberg lettuce" dans un nom de composé).
// Volontairement ABSENTS : SIRIUS, MS-FINDER, MetFrag, MAGMa, CSI:FingerID. Ce sont des outils d'IDENTIFICATION, cités
// dans les commentaires de spectres expérimentaux (MassBank "PROPOSED_COMPOUND SIRIUS", MoNA "MSFinder used to
// identify", GNPS "sirius atomic difference").
pub static IN_SILICO_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(&format!(
    r"(?i)in.silico|insilico|predicted|theoretical|{}", IN_SILICO_TOOLS
)).unwrap());

const IN_SILICO_TOOLS: &str =
    r"(?:^|[^a-z0-9])(?:cfm.?id|qce?ims|qcxms|neims|massformer|3dmolms|graff.?ms|fragnnet|isdb|lipidblast|(?-i:ICEBERG|SCARF))(?:[^a-z0-9]|$)";

// Formulations qui contiennent "predicted" ou "theoretical" sans dire que le SPECTRE est prédit :
// retirées du texte avant le test (le moteur de regex de Rust n'a pas de "sauf si suivi de").
// - "predicted rt rp=…", "predicted rt hilic=…", "predicted retention time", "predicted RI", "predicted CCS" (MoNA)
// - "predicted molecular formula SIRIUS" (GNPS)
// - "Annotation Predicted", "origin=MassBank; Predicted" : annotation des pics prédite (MassBank)
// - "theoretical m/z …, mass diff … (… ppm)" : m/z théorique d'une annotation (MoNA)
// - clauses d'IDENTIFICATION d'un spectre mesuré, qui citent les outils utilisés pour proposer un composé :
//   "PROPOSED_COMPOUND CFM-ID:MsfinderDB 1", "PROPOSED_MF SIRIUS5_Zodiac", "Proposed by SIRIUS5_Confidence",
//   "Supported by CFM-ID:MsfinderDB (Ranking 1), MS-FINDER:…", "Similarity by SIMILARITY:GNPS" (MassBank),
//   "annotation=MS-FINDER ver 3" (MoNA), "identified trough in silico generated spectrum" (MoNA, dans NAME).
//   Retirées jusqu'au prochain ";" (ou guillemet / parenthèse) : CFM-ID reste reconnu quand il décrit l'origine du spectre.
pub static IN_SILICO_EXCLUDE_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(
    r#"(?i)predicted\s+(?:rt|retention|ri\b|ccs|collision.cross.section|molecular\s+formula|formula)|annotation\s+predicted|origin=massbank;\s*predicted|theoretical\s+(?:m/z|mz|mass|exact\s+mass|precursor)|proposed_(?:compound|mf)\s+[^;"]*|(?:proposed|supported|similarity)\s+by\s+[^;")]*|annotation\s*=\s*[^;"]*|identified\s+(?:through|trough|thru|via|by|using|with)\s+(?:an?\s+)?in.silico"#
).unwrap());

// Marqueurs sans ambiguïté, cherchés aussi dans INSTRUMENTTYPE (ex. LipidBlast : "in-silico QTOF").
// Pas "predicted" ni "theoretical" dans ce champ : trop ambigus.
// Les noms d'outils de prédiction y sont aussi cherchés (MoNA : "Instrument_type: QCEIMS").
pub static IN_SILICO_STRONG_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(&format!(
    r"(?i)in.silico|insilico|{}", IN_SILICO_TOOLS
)).unwrap());

pub static RETENTION_TIME_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\W)?(m|min|minute|minutes|s|sec|second|seconds|ms|millisecond|milliseconds)(?:\W)?").unwrap());

// Regex optimisée pour capturer les modes d'ionisation avec détection des limites de mots (\b)
pub static IONIZATION_MODE_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)\b(APCI|ACPI|APPI|EI|ESI|FAB|MALDI)\b").unwrap());

pub static EMPTY_PATTERN: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"(?i)^(?:(?:CCS:( .*)?)|(?:\$:00in-source( .*)?)|(?:0( .*)?)|(?:0\.0( .*)?)|(?:)|(?:na( .*)?)|(?:n/a( .*)?)|(?:nan( .*)?)|(?:unknown( .*)?)|(?:unknow( .*)?)|(?:none( .*)?)|(?:\?( .*)?)|(?:unk( .*)?)|(?:x( .*)?))$").unwrap()
});

// =====================================

// ======================================================================================================================
// VARIABLES AJOUTÉES CÔTÉ RUST (N'existaient pas dans globals_vars.py)
// ======================================================================================================================

// Spectre GC (ou EI, ionisation électronique, en masse nominale et sans adduit) : "GC" ou "EI" en MOTS ENTIERS
// dans le type d'instrument ("GC-EI-TOF", "EI-B", "EI-EBEB", "GC-APCI-QTOF"…).
pub static GC_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)\b(?:GC|EI)\b").unwrap());

/// Règle UNIQUE « ce spectre est-il GC ? », utilisée partout (normalisation de l'adduit, contrôle de l'adduit,
/// test du m/z précurseur, correction de masse, tri final) :
/// - le nom du fichier contient "_GC" (convention du README pour un fichier qui ne contient que du GC), ou
/// - le type d'instrument contient "GC" ou "EI" en mots entiers.
/// Avant, chaque étape avait sa règle : le tri final rangeait "EI-B" en GC, mais le contrôle de l'adduit ne le
/// reconnaissait pas et supprimait ces spectres (sans adduit) ; le suffixe "_GC" n'était pas vu partout.
pub fn is_gc_spectrum(filename: &str, instrument_type: &str) -> bool {
    filename.contains("_GC") || GC_PATTERN.is_match(instrument_type)
}

// Regex pour extraire Collision Energy du Compound_Name (ex: GNPS, supporte 20-30, 102040/13334, 10 20 30)
pub static GNPS_COLLISION_ENERGY_EXTRACT_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?i)(?:\W|_)*collision[-\s_]?energy[\s:=]*(\d+(?:[\s.,\-/]+\d+)*\s*(?:eV|V)?)(?:\W|_)*").unwrap());

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_reconnu_par_le_type_d_instrument_ou_le_nom_de_fichier() {
        for it in ["GC-EI-TOF", "GC-APCI-QTOF", "EI-B", "EI-EBEB", "EI-QTOF", "GC-HRMS", "gc-ei-q", "GC"] {
            assert!(is_gc_spectrum("a.msp", it), "{it}");
        }
        assert!(is_gc_spectrum("GMD_20111121_MDN35_ALK_MSP_GC.msp", "")); // suffixe _GC du README
        for it in ["LC-ESI-QTOF", "Q-TOF", "Orbitrap", "ESI-ITFT", "HEIDELBERG", "GCMS", "LC-ESI-QQ", ""] {
            assert!(!is_gc_spectrum("a.msp", it), "{it}");
        }
    }
}
