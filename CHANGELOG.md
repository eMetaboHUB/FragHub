# Change Log
All notable changes to this project will be documented in this file.
- **27_01_2025**:
  - Refactoring updates and project logic.


- **24_01_2025**: 
  - deleting spectrum if neg adduct in pos spectrum, or pos adduct in neg spectrum
  - bug(fix): let thread safely terminate in no new spectrum to process
  - removing no or bad adduct spectrum
  - refactoring adduct normalization


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
