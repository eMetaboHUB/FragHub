// src/convertors/mod.rs

/// Module `convertors` : analyse les spectres bruts (lus par le module `loaders`) en structures Rust `Spectrum`.
pub mod csv_to_dict;
pub mod json_to_dict;
pub mod mgf_to_dict;
pub mod msp_to_dict;
pub mod parsing_to_dict;
pub mod raw_metadata;
