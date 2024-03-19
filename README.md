# FragHub  (1.0.0)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Required_Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10837523.svg)(https://doi.org/10.5281/zenodo.10837523)

## INSTALLATION

To install all dependencies, double-click on the install script corresponding to your OS.<br>
NB: **Make sure that python is in the path variables and that you run Python >= 3.9**<br>

## USE

To use this programme:

1) Put your msp, mgf, json, csv or xml files into 'INPUT/\<dedicated folder\>'.<br>
>NB: If you have a file that contains **only** In-Silico spectra AND this is not specified within the filename or the spectrum, you can simply suffix the filename with "_insilico", like this: "UNPD_ISDB_R_p01_insilico.mgf".<br>
2) Double-click on your corresponding OS run script into scripts folder.<br>
   1) FragHub GUI start<br>
   2) First tab: This area allows users to select specific functions for inclusion during the processing stage. Moreover, it provides the option to adjust the respective parameters of each function.<br>![img.png](img.png)
   3) Second tab: Select the output file format of your preference. By default, all formats are selected.<br>![img_1.png](img_1.png)
   4) Third tab: This tab facilitates the management of distinct profiles (like 'internal lab standards' or 'In-Silico DB', etc). Either select from a previously created profile or create a new one by just entering your desired profile name.<br>![img_2.png](img_2.png)
   4) Fourth tab: The 'Reset updates' option, when checked, allows for a reset of everything related to previously encountered spectra from the **current selected profil**. This will also delete all existing files located in the OUTPUT/{**current selected profil**} folder.<br>![img_3.png](img_3.png)
3) When the execution is complete, please remember to take a copy of your cleaned files from the OUTPUT folder and place them in a different location.
4) **DO NOT DELETE FILES INTO 'OUTPUT' AFTER COPY CLEANED VERSIONS.**

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

## required csv file
1) CSV files need to be separated by '**;**' with quotechar '**"**'.<br>
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
    