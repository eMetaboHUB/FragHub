# Change Log
All notable changes to this project will be documented in this file.

- **19_06_2026** (Version 2.0.0 Major Update):
  - **Rust Core Engine**: Completely rewrote the data processing pipeline (parsers, molecule derivation, filters, writers) in Rust via PyO3, replacing the pure Python implementation for massively increased high-throughput processing speed.
  - **Headless CLI Mode**: Introduced a new `--cli` mode to run FragHub without the graphical interface, ideal for servers and automated data pipelines. Features advanced real-time metrics (ETA, speed) and a `--quiet` flag.
  - **Recursive Directory Parsing**: Added input resolution in CLI mode to automatically scan entire folder structures and batch-process valid MS files.
  - **Backend Modernization**: Updated the FastAPI Python backend to use modern `lifespan` handlers, resolving deprecation warnings, and upgraded Pydantic schemas to v2 standards (`ConfigDict`).
  - **Agnostic Callback System**: Built a robust event bridge allowing the Rust engine to safely stream asynchronous progress logs to both the CLI terminal and the Electron/Nuxt GUI via WebSockets.
  - **Toolchain Migration**: Transitioned the entire Python environment and dependency management to `uv`, significantly speeding up setup and resolving packaging constraints with PyInstaller across OS environments.
  - **Interactive Reporting**: Developed a comprehensive HTML execution report incorporating interactive Sunburst and UpSetPlot charts that summarize dataset evolution after processing.
  - **UI/UX & Logs Polishing**: Disabled the automatic launch of Electron Developer Tools at startup, unified the application icons, and drastically cleaned/standardized all the real-time progress steps displayed during processing.
  - **Packaging Fixes**: Resolved Python multiprocessing freeze issues when packaging the Rust-RDKit bridge with PyInstaller, and updated Electron metadata.
  - **RDKit C++**: molecule derivation and exact masses are computed by RDKit C++ called directly from Rust, on all cores; RDKit Python is no longer used nor bundled.
  - **License**: FragHub is now distributed under the CeCILL 2.1 license (`LICENSE`), replacing CC BY-NC 4.0.

- **31_05_2026**:
  - modernizing GUI and desktop integration with Nuxt 4 + electron.js 

- **01_10_2025**:
  - fixing missing headers issue when writing a previously not existing csv file.
  - optimizing peaks filter with Numba just-in-time (jit)


- **26_09_2025**:
  - remove peaks with intensity <= 0.
  - bugfix with json update

- **15_09_2025**:
  - replace UNKNOWN by NOT FOUND for missing value
  - improve JSON output with pretty humane-readable JSON.
  - adding De Novo fragments formula calculations.
  - adding safety checks to the GUI start button.
  - adding file hash to spectrums (SHA-256 of the file size).
  - disabling automatic instrument deduction (e.g., inferring LC if ESI is found).

- **05_06_2025**:
  - replace RETENTIONTIME default field by RT.
  - Return to the main window when FINISH is pressed.
  - Correction of spectra moved to in-silico if FragHub output is used as input.
  - Deduplicate now on SPLASH **and** INCHIKEY. spectra without INCHIKEY are not deleted at the beginning of the process.

- **14_05_2025**:
  - fix adduct in silico correction

- **18_04_2025**:
  - fix adduct regex pattern
  - improve spectrum deletion callback by writing deleted spectrum in DELETION_REASONS sub folder with a detailed reason.
  - extend precursormz and adduct checks exception to all GC spectrums.
  - modifiy adduct ionmode check with "pos", "neg" in adduct dico.
  - direct integration of spectra-hash (https://github.com/berlinguyinca/spectra-hash) into fraghub.
  - Correcting and adduct some adducts to adduct dict
  - Adding loading screen to GUI
  - adding new gc instruments to instrument dict.
  - auto add [M+H]+ or [M-H]- in In-Silico if adduct is missing.
  - creating FragHub executable for Windows, Linux, and macOS with Python fully integrated.


- **03_04_2025**:
  - fixing missing last msp spectrum in msp files in some cases.


- **25_03_2025**:
  - fixing peaks list regex parsing error when formulas in peaks comment.
  

- **29_01_2025**:
  - changing licence from MIT to CC-BY-NC 4.0


- **27_01_2025**:
  - Refactoring updates and project logic.


- **24_01_2025**: 
  - deleting spectrum if neg adduct in pos spectrum, or pos adduct in neg spectrum
  - bug(fix): let thread safely terminate in no new spectrum to process
  - removing no or bad adduct spectrum
  - refactoring adduct normalization
  - adding report.txt in output directory


- **20_12_2024**: 
  - New modern GUI (Graphical User Interface).
  - Adding dependencies to requirements.txt
  - Removing MacOS and Linux support 😢 (this is just good-bye).
  - Refactoring duplicatas removal, now by sames SPLASH key.
  - Refactoring filtering logic with GC case.
  - Adding logs to GUI for monitoring suppressed spectrum.
  - Now completing NAMES and descriptor from RDkit and PubChem datas (offline).
  - Now completing Classyfire and NPclassifier from local datas.
  - Auto calculating chunk size for multi threading pool.
  - Now deleting spectrum with no SMILES no InChI **AND no inchikey**.
  - Refactoring .json reader for standard ISO/IEC 20802-2:2016 .json **and non-standard formats**.
  - Moving all globals variables to a single file.
