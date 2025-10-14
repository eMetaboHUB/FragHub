# Change Log
All notable changes to this project will be documented in this file.

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
