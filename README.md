<div align="center">
  <img src="scripts/GUI/assets/FragHub_icon.png" alt="FragHub Icon" /> 2.0.0
</div>

# FragHub  (2.0.0)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white&style=for-the-badge)
![Rust](https://img.shields.io/badge/Rust-1.96.0-000000?logo=rust&logoColor=white&style=for-the-badge)
![Electron](https://img.shields.io/badge/Electron-40.4.0-47848F?logo=electron&logoColor=white&style=for-the-badge)
![Nuxt](https://img.shields.io/badge/Nuxt-4.4.7-00C58E?logo=nuxt&logoColor=white&style=for-the-badge)
![Vite](https://img.shields.io/badge/Vite-8.0.16-646CFF?logo=vite&logoColor=white&style=for-the-badge)
![npm](https://img.shields.io/badge/npm-11.13.0-CB3837?logo=npm&logoColor=white&style=for-the-badge)
![License: CeCILL](https://img.shields.io/badge/License-CeCILL-yellow.svg?style=for-the-badge)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17235587-1f73b7?style=for-the-badge)](https://doi.org/10.5281/zenodo.17235587)

## INSTALL INSTRUCTIONS BELLOW
[See install instructions](#installation)

## DESCRIPTION
![img_10.png](img_10.png)
FragHub is a powerful tool designed to standardize and organize mass spectrometry (MS) data from OMSLs (Open Mass Spectra Libraries). The main objective of FragHub is to simplify and improve the process of MS data analysis by providing standardized, consistent, and easily accessible data.<br>
![img_11.png](img_11.png)
Key features:<br>
- **High-throughput processing**: Powered by a Rust-based core engine (via PyO3), FragHub 2.0 efficiently processes large-scale mass spectrometry libraries with significantly reduced computation times.
- **Data standardization**: FragHub standardizes field names and values of MS spectra from various databases, ensuring data consistency and compatibility.![img_12.png](img_12.png)![img_7.png](img_7.png)
- **Optional peak list filtering**: FragHub applies filters to peak lists to streamline them by retaining only essential data, significantly reducing data size without compromising relevance for analysis.
- **Recalculation and normalization of chemical identifiers**: The program recalculates and normalizes chemical identifiers such as SMILES, InChI, and InChIKey, ensuring a uniform representation of molecular structures.
  - Salts are removed (desalting): for a structure with several components (e.g. "quinine sulfate"), the component matching the precursor m/z is kept, since only that molecule is measured.
  - Charges are neutralized when possible.
  - **InChI and InChIKey** are the standard IUPAC identifiers of this desalted, neutralized molecule, so they can be compared with PubChem and other databases.
  - **Output SMILES are normalized with RDKit**: they are written in the RDKit canonical tautomer form (stereochemistry kept), so that the same molecule is always written the same way. When the canonical tautomer would describe a different molecule (different InChI, e.g. keto/enol), the SMILES of the molecule as given is written instead. Output SMILES may therefore differ from the SMILES of the input files.
- **Ontologies**: FragHub complete ontologic information from PubChem datas bases on Inchikey.
- **Spectra organization**: FragHub separates spectra based on different experimental parameters such as polarity (positive/negative), chromatographic mode (LC/GC), and acquisition type (experimental/in silico), facilitating their use and subsequent analysis.![img_8.png](img_8.png)
- **Compatibility with various analysis software**: Standardized spectra produced by FragHub are compatible with multiple analysis software, including MSdial, MZmine, and Flash Entropy Search, providing users with maximum flexibility in choosing analysis tools.

**Warning**: All spectra deemed inconsistent, i.e., those lacking SMILES and InChI, precursor m/z, and adduct information, are removed during the processing by FragHub.![img_9.png](img_9.png)


## CONTACT

* **Developer :** <axel.dablanc@utoulouse.fr>
* **Project Manager :** <guillaume.marti@utoulouse.fr>


## INSTALLATION
1) Download FragHub executable corresponding to your OS HERE:<br>
   ~~**Windows** x64: [Download Windows x64 Installer](https://zenodo.org/records/17235741/files/FragHub_2.0.0_setup_windows_x64.exe?download=1)<br>
   **Linux** x64: [Download Linux x64 Installer](https://zenodo.org/records/17235741/files/FragHub_2.0.0_setup_Linux_x64?download=1)<br>
   **macOS** arm64: [Download macOS arm64 app](https://zenodo.org/records/17235741/files/FragHub_macOS_2.0.0_arm64.dmg?download=1)<br>~~
2) Execute setup installer (macOS: drag to App)
3) select installation folder
4) click install

## USE

To use this programme:

1) Go to install directory then run FragHub executable (feel free to create shortcut).
>NB: If you have a file that contains **only** In-Silico spectra AND this is not specified within the filename or the spectrum, you can simply suffix the filename with "_insilico", like this: "UNPD_ISDB_R_p01_insilico.mgf".<br>
<br>
>NB: If you have a file that contains **only** GC spectra AND this is not specified within the filename or the spectrum, you can simply suffix the filename with "_GC", like this: "GMD_20111121_MDN35_ALK_MSP_GC.msp".<br>

2) FragHub GUI start: **at any time you can move the mouse over 🛈 to view a detailed description**<br>
   1) First tab: Select single or multiple .json, .csv, .msp, or .mgf files.<br>![img.png](img.png)
   2) Second tab: Create a new empty directory or Select an existing directory where FragHub has already written files<br>![img_1.png](img_1.png)
   3) Third tab: This area allows users to select specific filters for inclusion during the processing stage. Moreover, it provides the option to adjust the respective parameters of each filters.<br>![img_2.png](img_2.png)
   4) Fourth tab: This area allows the user to choose whether or not to perform the calculation of de novo fragment formulas.<br>![img_13.png](img_13.png)
   4) Fifth tab: Select the output file format of your preference. By default, all formats are selected.<br>![img_3.png](img_3.png)
   5) Sixth tab: FragHub save the splash keys of previous spectra processed.<br>So that at the next update, only new spectra from the database are processed, and added to previous FragHub processes.<br>Reseting project allows you to delete splash keys and output files from the selected project, in order to start a new project from scratch.<br>The "Write perf logs" switch (off by default) writes a performance log of the run (duration, memory and CPU usage of each step) in the PERF_LOGS sub-folder of the output directory.<br>![img_4.png](img_4.png)
   6) Click on START button to run the programme.
3) Progress start. At any time you can click on STOP button to close the programme.<br>![img_5.png](img_5.png)
4) When Progress is finish, click on FINISH button to close the programme.<br>![img_6.png](img_6.png)
5) **DO NOT DELETE FILES INTO YOUR SELECTED OUTPUT DIRECTORY AFTER COPY CLEANED VERSIONS.**
---

## TECHNICAL DOCUMENTATION

Please refer to the [technical documentation](documentation/) for comprehensive details regarding the architecture, inner workings, and development of FragHub.

---

## HEADLESS CLI MODE

FragHub 2.0 introduces a headless Command Line Interface (CLI). This feature enables seamless integration into automated data pipelines and server environments, utilizing the Rust core engine without requiring a graphical user interface.

Please refer to the [CLI Documentation](documentation/Documentation_CLI.md) for comprehensive usage instructions, available arguments, and pipeline integration examples.

---

## FILTERS

**check_minimum_peak_requiered(peak_array, n_peaks)**<br>
This function checks whether a given mass spectrum contains a minimum number of peaks. If the spectrum contains fewer peaks than the minimum requirement, it ignores the spectrum.<br>
<br>
**remove_peak_above_precursormz(peak_array, precursormz)**<br>
This function removes all peaks from the spectrum whose m/z value is greater than the precursor's m/z value plus 5 Da.<br>
<br>
**reduce_peak_list(peak_array, max_peaks)**<br>
This function reduces the peak list to a specified maximum number of peaks. The peaks to retain are chosen based on their intensity, with peaks of greater intensity being selected.<br>
<br>
**normalize_intensity(peak_array)**<br>
This function normalizes the intensity of all the peaks in a given spectrum to the maximum intensity.<br>
<br>
**keep_mz_in_range(peak_array, mz_from, mz_to)**<br>
This function takes an array of peak data (representing mass-to-charge ratio, or m/z) and returns a new array containing only those peaks whose m/z value falls between mz_from and mz_to.<br>

**check_minimum_of_high_peaks_requiered(peak_array, intensity_percent, no_peaks)**<br>
This function is used to check whether a given array containing peak data has a required minimum number of "high peaks". A "high peak" is defined as a peak whose intensity is above a certain percentage (intensity_percent) of the maximum intensity. If the array does not contain a sufficient number of "high peaks", the function ignore the spectrum.<br>

**remove_spectrum_under_entropy_score(score)**<br>
The Shannon entropy score of the spectrum is calculated during processing. If a spectrum has an entropy score lower than the minimum required, it is deleted.<br>

## required csv file
1) CSV files need to be separated by '**;**' or '**\t**' with quotechar '**"**'.<br>
2) peaks columns need to be named '**peaks**'.<br>
3) '**peaks**' column need to be formatted with one of the following format, in string:
   1) >"[[79.054840, 12486.074219], [79.629868, 854.089905]]"
   2) > "<br>
   57.07042529 0.7697591662<br>
   71.08607535 1.507457981<br>
   97.06533991 0.4893302623<br>
   99.08098997 0.4737337839<br>
   137.09664 0.498920401<br>
   165.0915547 0.4243093978<br>
   "<br>
   3) > "314.00000000 1.00000000;315.00000000 0.15715716;316.00000000 0.35735736;317.00000000 0.05105105"
    