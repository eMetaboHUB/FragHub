# FragHub  (4.0.0)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Required_Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)


## INSTALLATION

To install all dependencies, double-click on "Install.bat".<br>
NB: Make sure that python is in the path variables and that you run Python >= 3.9<br>

## USE

To use this programme:

1) Put your msp, json, csv or xml files into "INPUT/\<dedicated folder\>".
2) Double-click on your corresponding OS run script into RUN folder.<br>
3) when the execution is finished, copy your cleaned msp files <br>"OUTPUT/MSP/(POS or NEG)"<br>and<br>csv version of the cleaned msp (POS and NEG) into<br>"OUTPUT/CSV/(POS or NEG)"
4) DO NOT DELETE FILES INTO "OUTPUT" AFTER COPY CLEANED VERSIONS.


## required csv file
1) CSV files need to be separated by ";" with quotechar '"'.<br>
2) peaks columns need to be named "peaks".<br>
3) "peaks" column need to be formatted with one of the following format, in string:
   1) >"[[79.054840, 12486.074219], [79.629868, 854.089905]]"
   2) > "<br>
   57.07042529 0.7697591662<br>
   71.08607535 1.507457981<br>
   97.06533991 0.4893302623<br>
   99.08098997 0.4737337839<br>
   137.09664 0.498920401<br>
   165.0915547 0.4243093978<br>
   "<br>
    