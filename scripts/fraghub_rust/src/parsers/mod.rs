// src/parsers/mod.rs

/// Module `parsers` : analyse les spectres bruts (lus par le module `loaders`) en structures Rust `Spectrum`,
/// un parseur par format. `dispatcher` envoie chaque fichier d'entrée vers le bon parseur ; `raw_metadata`
/// renomme les champs bruts vers les clés FragHub (règles communes à tous les formats).
pub mod csv_parser;
pub mod json_parser;
pub mod mgf_parser;
pub mod msp_parser;
pub mod dispatcher;
pub mod raw_metadata;