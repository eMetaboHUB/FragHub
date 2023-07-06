# MSP_V3 (MultiThreaded) (3.6.0)



## INSTALLATION

To install all dependencies, double click on "Install.bat".<br>
NB: Make sure that python is in the path variables.<br>
Copie the 'known_key_conversions.csv' into the folder 'data' of matchms installation. Replace it.

## USE

To use this programme:

1) Put your msp, json or xml files into "INPUT/\<dedicated folder\>".
2.1) If you have HMDB xml spectrums, run HMDB_completion.bat first.
2.2) Else: Double click on "Use.bat".
3) when the execution is finished, retrieve your cleaned msp files into <br>"OUTPUT/CLEAN_MSP/(FINAL_POS or FINAL_NEG)"<br>and<br>csv version of the cleaned msp (POS and NEG) into<br>"OUTPUT/CSV/(FINAL_POS or FINAL_NEG)"
4) DO NOT FORGET TO DELETE FILES INTO "INPUT" and "OUTPUT" AFTER RETRIEVING CLEANED VERSIONS.