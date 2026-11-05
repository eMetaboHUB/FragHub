//! build_dataset : fabrique un jeu de données de test plus petit à partir d'un dossier d'entrée FragHub.
//!
//! Pour chaque fichier .msp, .mgf, .json et .csv du dossier, on garde une partie des spectres,
//! tirés au hasard sur tout le fichier (donc représentatifs), et on l'écrit dans un fichier de
//! même nom et de même format dans le dossier de sortie.
//!
//! Deux modes de tirage :
//! - `--percent P` (défaut : 20) : chaque spectre est gardé avec une probabilité de P %.
//!   Le tri se fait au fil de la lecture : rien n'est gardé en mémoire, même pour 20 % de 11 Go.
//!   Le nombre gardé est proche de P % du total (tirage indépendant pour chaque spectre).
//! - `--per-file N` : exactement N spectres par fichier (ou tous s'il y en a moins),
//!   tirage "réservoir" uniforme. Les N spectres sont gardés en mémoire avant écriture.
//!
//! Les spectres sont recopiés OCTET POUR OCTET (aucune réécriture), dans leur ordre d'origine, et
//! découpés avec les mêmes règles que les lecteurs de FragHub (scripts/fraghub_rust/src/convertors/loaders.rs) :
//! - MSP  : un spectre = un bloc de lignes non vides, séparé par une ligne vide ;
//! - MGF  : un spectre = toutes les lignes jusqu'à "END IONS" inclus ;
//! - JSON : un spectre = un objet {...} de premier niveau (format tableau GNPS ou JSONL) ;
//! - CSV  : la ligne d'en-tête + un enregistrement par spectre (les champs entre guillemets
//!          sur plusieurs lignes sont gérés).
//! Le tirage est reproductible : même graine (--seed) = même échantillon.
//!
//! Usage :
//!   cargo run --release -- <dossier_entrée> [dossier_sortie] [--percent P | --per-file N] [--seed S]
//! Par défaut : 20 % des spectres de chaque fichier, graine 42, sortie = "<dossier_entrée>-sample".

use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Tirage aléatoire reproductible (sans dépendance externe)
// ---------------------------------------------------------------------------

/// Petit générateur pseudo-aléatoire (SplitMix64) : rapide, reproductible, suffisant pour un tirage.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Entier uniforme dans 0..=max.
    fn below_or_equal(&mut self, max: u64) -> u64 {
        // Rejet pour éviter le biais du modulo
        let n = max.wrapping_add(1);
        if n == 0 {
            return self.next_u64();
        }
        let zone = u64::MAX - (u64::MAX % n);
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % n;
            }
        }
    }

    /// Nombre réel uniforme dans [0, 1).
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ---------------------------------------------------------------------------
// Écriture du fichier de sortie, selon le format
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum Format {
    Msp,
    Mgf,
    Json,
    Csv,
}

/// Écrit les spectres retenus dans le fichier de sortie, avec les séparateurs du format.
struct OutWriter {
    w: BufWriter<File>,
    format: Format,
    /// Fin de ligne du fichier d'origine (CRLF ou LF), pour écrire les séparateurs dans le même style
    nl: &'static [u8],
    written: usize,
    /// En-tête CSV (écrit avant le premier enregistrement)
    header: Vec<u8>,
    started: bool,
}

impl OutWriter {
    fn create(path: &Path, format: Format) -> io::Result<Self> {
        Ok(OutWriter {
            w: BufWriter::with_capacity(8 * 1024 * 1024, File::create(path)?),
            format,
            nl: b"\n",
            written: 0,
            header: Vec::new(),
            started: false,
        })
    }

    fn start(&mut self) -> io::Result<()> {
        if self.started {
            return Ok(());
        }
        self.started = true;
        match self.format {
            Format::Json => self.w.write_all(b"[\n")?,
            Format::Csv => {
                self.w.write_all(&self.header)?;
                if !self.header.is_empty() && !self.header.ends_with(b"\n") {
                    self.w.write_all(self.nl)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn block(&mut self, b: &[u8]) -> io::Result<()> {
        self.start()?;
        match self.format {
            Format::Msp => {
                self.w.write_all(b)?;
                if !b.ends_with(b"\n") {
                    self.w.write_all(self.nl)?; // dernière ligne du fichier sans retour à la ligne
                }
                self.w.write_all(self.nl)?; // ligne vide de séparation
            }
            Format::Mgf | Format::Csv => {
                self.w.write_all(b)?;
                if !b.ends_with(b"\n") {
                    self.w.write_all(self.nl)?;
                }
            }
            Format::Json => {
                if self.written > 0 {
                    self.w.write_all(b",\n")?;
                }
                self.w.write_all(b)?;
            }
        }
        self.written += 1;
        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        self.start()?;
        if self.format == Format::Json {
            self.w.write_all(b"\n]\n")?;
        }
        self.w.flush()
    }
}

// ---------------------------------------------------------------------------
// Sélection des spectres
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
enum Mode {
    /// Pourcentage de spectres (0 < p <= 100), tirage indépendant pour chacun
    Percent(f64),
    /// Nombre fixe de spectres par fichier (tirage réservoir)
    Count(usize),
}

/// Décide quels spectres garder. Les lecteurs de format l'appellent ainsi, pour chaque spectre :
/// `begin()` (au début), `push()` (ses octets), puis `commit()` (spectre complet) ou `abort()`
/// (faux spectre, par exemple une fin de fichier incomplète).
/// La décision est prise dans `begin()` : un spectre non retenu n'est jamais copié en mémoire.
struct Selector {
    mode: Mode,
    rng: Rng,
    seen: u64,
    kept: usize,
    /// Mode pourcentage : octets du spectre en cours, s'il est retenu
    current: Option<Vec<u8>>,
    /// Mode nombre fixe : (numéro du spectre, octets) des spectres retenus
    slots: Vec<(u64, Vec<u8>)>,
    current_slot: Option<usize>,
    /// Spectre remplacé par le spectre en cours, gardé de côté tant que celui-ci n'est pas complet
    evicted: Option<(u64, Vec<u8>)>,
}

impl Selector {
    fn new(mode: Mode, seed: u64) -> Self {
        Selector { mode, rng: Rng(seed), seen: 0, kept: 0, current: None, slots: Vec::new(), current_slot: None, evicted: None }
    }

    /// Début d'un nouveau spectre.
    fn begin(&mut self) {
        let i = self.seen;
        self.seen += 1;
        match self.mode {
            Mode::Percent(p) => {
                self.current = if self.rng.unit() * 100.0 < p { Some(Vec::new()) } else { None };
            }
            Mode::Count(k) => {
                self.current_slot = if (i as usize) < k {
                    self.slots.push((i, Vec::new()));
                    Some(self.slots.len() - 1)
                } else {
                    let j = self.rng.below_or_equal(i) as usize;
                    if j < k {
                        self.evicted = Some(std::mem::replace(&mut self.slots[j], (i, Vec::new())));
                        Some(j)
                    } else {
                        None
                    }
                };
            }
        }
    }

    fn push(&mut self, bytes: &[u8]) {
        if let Some(buf) = self.current.as_mut() {
            buf.extend_from_slice(bytes);
        } else if let Some(slot) = self.current_slot {
            self.slots[slot].1.extend_from_slice(bytes);
        }
    }

    /// Le spectre en cours est complet : il est validé (et écrit tout de suite en mode pourcentage).
    fn commit(&mut self, out: &mut OutWriter) -> io::Result<()> {
        if let Some(buf) = self.current.take() {
            out.block(&buf)?;
            self.kept += 1;
        }
        self.current_slot = None;
        self.evicted = None;
        Ok(())
    }

    /// Le spectre en cours n'était pas un vrai spectre : on annule comme s'il n'avait pas existé.
    fn abort(&mut self) {
        self.current = None;
        if let Some(slot) = self.current_slot.take() {
            match self.evicted.take() {
                Some(old) => self.slots[slot] = old, // on remet le spectre remplacé
                None => {
                    self.slots.pop(); // il venait d'être ajouté à la fin
                }
            }
        }
        self.seen -= 1;
    }

    /// Fin du fichier : en mode nombre fixe, écrit les spectres retenus dans leur ordre d'origine.
    fn finish(mut self, mut out: OutWriter) -> io::Result<Sample> {
        if let Mode::Count(_) = self.mode {
            self.slots.sort_by_key(|(i, _)| *i);
            for (_, b) in &self.slots {
                out.block(b)?;
            }
            self.kept = self.slots.len();
        }
        out.finish()?;
        Ok(Sample { total: self.seen, kept: self.kept })
    }
}

/// Résultat de l'échantillonnage d'un fichier.
struct Sample {
    total: u64,
    kept: usize,
}

// ---------------------------------------------------------------------------
// Outils de lecture
// ---------------------------------------------------------------------------

/// Même test que les lecteurs FragHub : `line.trim()` (espaces Unicode compris).
fn trimmed(line: &[u8]) -> std::borrow::Cow<'_, str> {
    match std::str::from_utf8(line) {
        Ok(s) => std::borrow::Cow::Borrowed(s.trim()),
        Err(_) => std::borrow::Cow::Owned(String::from_utf8_lossy(line).trim().to_string()),
    }
}

/// Fin de ligne utilisée par une ligne (CRLF ou LF).
fn newline_of(line: &[u8]) -> &'static [u8] {
    if line.ends_with(b"\r\n") { b"\r\n" } else { b"\n" }
}

fn open_reader(path: &Path) -> io::Result<BufReader<File>> {
    Ok(BufReader::with_capacity(8 * 1024 * 1024, File::open(path)?))
}

// ---------------------------------------------------------------------------
// MSP : blocs de lignes non vides séparés par une ligne vide
// ---------------------------------------------------------------------------
fn sample_msp(input: &Path, output: &Path, mode: Mode, seed: u64) -> io::Result<Sample> {
    let mut reader = open_reader(input)?;
    let mut out = OutWriter::create(output, Format::Msp)?;
    let mut sel = Selector::new(mode, seed);
    let mut nl_known = false;
    let mut in_block = false;
    let mut line = Vec::with_capacity(4096);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if !nl_known {
            out.nl = newline_of(&line);
            nl_known = true;
        }
        if trimmed(&line).is_empty() {
            if in_block {
                sel.commit(&mut out)?; // fin du spectre en cours
            }
            in_block = false;
        } else {
            if !in_block {
                sel.begin();
                in_block = true;
            }
            sel.push(&line);
        }
    }
    if in_block {
        sel.commit(&mut out)?; // dernier spectre, sans ligne vide après
    }
    sel.finish(out)
}

// ---------------------------------------------------------------------------
// MGF : toutes les lignes jusqu'à "END IONS" inclus
// ---------------------------------------------------------------------------
fn sample_mgf(input: &Path, output: &Path, mode: Mode, seed: u64) -> io::Result<Sample> {
    let mut reader = open_reader(input)?;
    let mut out = OutWriter::create(output, Format::Mgf)?;
    let mut sel = Selector::new(mode, seed);
    let mut nl_known = false;
    let mut started = false;
    let mut line = Vec::with_capacity(4096);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if !nl_known {
            out.nl = newline_of(&line);
            nl_known = true;
        }
        if !started {
            sel.begin();
            started = true;
        }
        sel.push(&line);
        if trimmed(&line) == "END IONS" {
            sel.commit(&mut out)?;
            started = false; // le spectre suivant commence à la ligne d'après
        }
    }
    // Comme FragHub : ce qui suit le dernier "END IONS" n'est pas un spectre.
    if started {
        sel.abort();
    }
    sel.finish(out)
}

// ---------------------------------------------------------------------------
// JSON : objets {...} de premier niveau (même automate que load_spectrum_list_json)
// ---------------------------------------------------------------------------
fn sample_json(input: &Path, output: &Path, mode: Mode, seed: u64) -> io::Result<Sample> {
    let mut reader = open_reader(input)?;
    let mut out = OutWriter::create(output, Format::Json)?;
    let mut sel = Selector::new(mode, seed);
    let mut depth: u64 = 0;
    let mut in_string = false;
    let mut escaped = false;
    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        let n = buf.len();
        // Début de l'objet en cours dans ce morceau (pour recopier ses octets d'un coup)
        let mut start_in_chunk = if depth > 0 { Some(0) } else { None };
        for (pos, &b) in buf.iter().enumerate() {
            if depth == 0 && b == b'{' && !in_string {
                sel.begin();
                start_in_chunk = Some(pos);
            }
            if in_string {
                if escaped {
                    escaped = false;
                } else if b == b'\\' {
                    escaped = true;
                } else if b == b'"' {
                    in_string = false;
                }
            } else {
                match b {
                    b'"' => in_string = true,
                    b'{' => depth += 1,
                    b'}' => {
                        if depth > 0 {
                            depth -= 1;
                            if depth == 0 {
                                sel.push(&buf[start_in_chunk.unwrap_or(0)..=pos]);
                                sel.commit(&mut out)?;
                                start_in_chunk = None;
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        // Objet commencé mais pas fini dans ce morceau : on garde la partie déjà lue
        if depth > 0 {
            sel.push(&buf[start_in_chunk.unwrap_or(0)..n]);
        }
        reader.consume(n);
    }
    // Objet non terminé en fin de fichier : FragHub l'ignore, on fait pareil.
    if depth > 0 {
        sel.abort();
    }
    sel.finish(out)
}

// ---------------------------------------------------------------------------
// CSV : en-tête + enregistrements (un champ entre guillemets peut contenir des retours à la ligne)
// ---------------------------------------------------------------------------
fn sample_csv(input: &Path, output: &Path, mode: Mode, seed: u64) -> io::Result<Sample> {
    let mut reader = open_reader(input)?;
    let mut out = OutWriter::create(output, Format::Csv)?;
    let mut sel = Selector::new(mode, seed);
    let mut header_done = false;
    let mut record_open = false; // enregistrement en cours
    let mut quotes_odd = false; // guillemets non refermés
    let mut line = Vec::with_capacity(4096);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let odd_in_line = line.iter().filter(|&&c| c == b'"').count() % 2 == 1;
        if !header_done {
            out.header.extend_from_slice(&line);
            quotes_odd ^= odd_in_line;
            if !quotes_odd {
                header_done = true;
                out.nl = newline_of(&line);
            }
            continue;
        }
        if !record_open {
            if trimmed(&line).is_empty() {
                continue; // ligne vide entre deux enregistrements
            }
            sel.begin();
            record_open = true;
            quotes_odd = false;
        }
        sel.push(&line);
        quotes_odd ^= odd_in_line;
        if !quotes_odd {
            sel.commit(&mut out)?;
            record_open = false;
        }
    }
    // Enregistrement entre guillemets jamais refermé : on le garde tel quel (le lecteur CSV décidera)
    if record_open {
        sel.commit(&mut out)?;
    }
    sel.finish(out)
}

// ---------------------------------------------------------------------------
// Programme principal
// ---------------------------------------------------------------------------

struct Args {
    input: PathBuf,
    output: PathBuf,
    mode: Mode,
    seed: u64,
}

const USAGE: &str = "Usage : build_dataset <dossier_entrée> [dossier_sortie] [--percent P | --per-file N] [--seed S]
  --percent P  : pourcentage de spectres gardés dans chaque fichier, tirés au hasard (défaut : 20)
  --per-file N : nombre fixe de spectres gardés par fichier, tirés au hasard
  --seed S     : graine du tirage aléatoire (défaut : 42)
  Sortie par défaut : <dossier_entrée>-sample";

fn parse_args() -> Result<Args, String> {
    let mut positional = Vec::new();
    let mut mode: Option<Mode> = None;
    let mut seed = 42u64;
    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "-h" | "--help" => return Err(String::new()),
            "--percent" => {
                if mode.is_some() {
                    return Err("--percent et --per-file ne peuvent pas être utilisés ensemble".into());
                }
                let p: f64 = it.next().map(|v| v.trim_end_matches('%').replace(',', "."))
                    .and_then(|v| v.parse().ok()).filter(|&v: &f64| v > 0.0 && v <= 100.0)
                    .ok_or("--percent attend un nombre entre 0 (exclu) et 100, par exemple 20")?;
                mode = Some(Mode::Percent(p));
            }
            "--per-file" => {
                if mode.is_some() {
                    return Err("--percent et --per-file ne peuvent pas être utilisés ensemble".into());
                }
                let k = it.next().and_then(|v| v.parse().ok()).filter(|&v| v > 0)
                    .ok_or("--per-file attend un nombre entier > 0")?;
                mode = Some(Mode::Count(k));
            }
            "--seed" => {
                seed = it.next().and_then(|v| v.parse().ok()).ok_or("--seed attend un nombre entier")?;
            }
            _ if a.starts_with("--") => return Err(format!("option inconnue : {a}")),
            _ => positional.push(PathBuf::from(a)),
        }
    }
    if positional.len() > 2 {
        return Err("trop d'arguments".into());
    }
    let input = positional.first().cloned().ok_or("dossier d'entrée manquant")?;
    let output = match positional.get(1) {
        Some(p) => p.clone(),
        None => {
            let name = input.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_else(|| "input".into());
            input.with_file_name(format!("{name}-sample"))
        }
    };
    Ok(Args { input, output, mode: mode.unwrap_or(Mode::Percent(20.0)), seed })
}

fn fmt_size(bytes: u64) -> String {
    let mb = bytes as f64 / (1024.0 * 1024.0);
    if mb >= 1024.0 { format!("{:.2} Go", mb / 1024.0) } else { format!("{:.1} Mo", mb) }
}

fn describe(mode: Mode) -> String {
    match mode {
        Mode::Percent(p) => format!("{p} % des spectres de chaque fichier"),
        Mode::Count(k) => format!("{k} spectres par fichier"),
    }
}

fn run(args: &Args) -> Result<(), String> {
    if !args.input.is_dir() {
        return Err(format!("dossier d'entrée introuvable : {}", args.input.display()));
    }
    // Sécurité : ne jamais écrire dans le dossier d'entrée (on écraserait les originaux)
    let input_abs = fs::canonicalize(&args.input).map_err(|e| e.to_string())?;
    fs::create_dir_all(&args.output).map_err(|e| format!("création de {} : {e}", args.output.display()))?;
    let output_abs = fs::canonicalize(&args.output).map_err(|e| e.to_string())?;
    if output_abs == input_abs || output_abs.starts_with(&input_abs) {
        return Err("le dossier de sortie doit être en dehors du dossier d'entrée".into());
    }

    let mut files: Vec<PathBuf> = fs::read_dir(&args.input).map_err(|e| e.to_string())?
        .flatten().map(|e| e.path()).filter(|p| p.is_file()).collect();
    files.sort();

    println!("Entrée  : {}", args.input.display());
    println!("Sortie  : {}", args.output.display());
    println!("Tirage  : {}, au hasard, graine {}\n", describe(args.mode), args.seed);

    let mut summary = format!("Échantillon créé par scripts/build_dataset\nSource : {}\nTirage : {}\nGraine : {}\n\n",
                              args.input.display(), describe(args.mode), args.seed);
    let start_all = Instant::now();
    let (mut total_in, mut total_out) = (0u64, 0u64);
    for path in &files {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        // Même règle que FragHub (parsing_to_dict.rs) : extension en minuscules, sensible à la casse
        let sampler: fn(&Path, &Path, Mode, u64) -> io::Result<Sample> =
            if name.ends_with(".msp") { sample_msp }
            else if name.ends_with(".mgf") { sample_mgf }
            else if name.ends_with(".json") { sample_json }
            else if name.ends_with(".csv") { sample_csv }
            else {
                println!("  (ignoré, format non lu par FragHub) {name}");
                continue;
            };
        let out_path = args.output.join(&name);
        let size_in = fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        print!("  {name} ({}) ... ", fmt_size(size_in));
        io::stdout().flush().ok();
        let t = Instant::now();
        match sampler(path, &out_path, args.mode, args.seed) {
            Ok(s) => {
                let size_out = fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0);
                total_in += size_in;
                total_out += size_out;
                let pct = if s.total > 0 { 100.0 * s.kept as f64 / s.total as f64 } else { 0.0 };
                let line = format!("{} / {} spectres ({:.1} %), {} ({:.1} s)",
                                   s.kept, s.total, pct, fmt_size(size_out), t.elapsed().as_secs_f64());
                println!("{line}");
                summary.push_str(&format!("{name} : {line}\n"));
            }
            Err(e) => {
                println!("ERREUR : {e}");
                let _ = fs::remove_file(&out_path);
                summary.push_str(&format!("{name} : ERREUR {e}\n"));
            }
        }
    }
    // Résumé à côté des fichiers (extension .txt : ignorée par FragHub)
    let _ = fs::write(args.output.join("_sample_info.txt"), &summary);
    println!("\nTerminé en {:.1} s : {} → {}", start_all.elapsed().as_secs_f64(), fmt_size(total_in), fmt_size(total_out));
    Ok(())
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(msg) => {
            if !msg.is_empty() {
                eprintln!("Erreur : {msg}\n");
            }
            eprintln!("{USAGE}");
            std::process::exit(2);
        }
    };
    if let Err(e) = run(&args) {
        eprintln!("Erreur : {e}");
        std::process::exit(1);
    }
}
