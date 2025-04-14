# Change Log
All notable changes to this project will be documented in this file.

- **14_04_2025**:
  - fix adduct regex pattern
  - improve spectrum deletion callback by writing deleted spectrum in DELETION_REASONS sub folder with a detailed reason.
  - extend precursormz and adduct checks exception to all GC spectrums.
  - modifiy adduct ionmode check with "pos", "neg" in adduct dico.


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
