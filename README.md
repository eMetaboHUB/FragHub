# FragHub  (4.0.0)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Required_Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)


## INSTALLATION

To install all dependencies, double-click on the install script corresponding to your OS.<br>
NB: **Make sure that python is in the path variables and that you run Python >= 3.9**<br>

## USE

To use this programme:

1) Put your msp, mgf, json, csv or xml files into 'INPUT/\<dedicated folder\>'.
2) Double-click on your corresponding OS run script into RUN folder.<br>
   1) FragHub GUI start<br>
   2) First tab: The 'Reset updates' option, when checked, allows for a reset of everything related to previously encountered spectra. This will also delete all existing files located in the OUTPUT folder.<br>![img.png](img.png)
   3) Second tab: This area allows users to select specific functions for inclusion during the processing stage. Moreover, it provides the option to adjust the respective parameters of each function.<br>![img_1.png](img_1.png)
   4) Third tab: Select the output file format of your preference. By default, all formats are selected.<br>![img_2.png](img_2.png)
3) When the execution is complete, please remember to take a copy of your cleaned files from the OUTPUT folder and place them in a different location.
4) **DO NOT DELETE FILES INTO 'OUTPUT' AFTER COPY CLEANED VERSIONS.**


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
    