# Change Log
All notable changes to this project will be documented in this file.


- **20_12_2024**: 
  - New modern GUI (Graphical User Interface).
  - Adding dependencies to requirements.txt
  - Removing MacOS and Linux support 😢 (this is just good-bye).
  - Refactoring duplicatas removal, now by sames SPLASH key.
  - Refactoring filtering logic with GC case.
  - Adding logs to GUI for monitoring suppressed spectrum.
  - Now completing NAMES and descriptor from RDkit and PubChem datas (offline).
  - Now completing Classyfire and NPclassifier from ...
  - Auto calculating chunk size for multi threading pool.
  - Now deleting spectrum with no SMILES no InChI **AND no inchikey**.
  - Refactoring .json reader for standard ISO/IEC 20802-2:2016 .json **and non-standard formats**.
  - Moving all globals variables to a single file.
