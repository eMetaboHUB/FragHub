# FragHub  (4.0.0)



## INSTALLATION

To install all dependencies, double click on "Install.bat".<br>
NB: Make sure that python is in the path variables.<br>
Copie the 'known_key_conversions.csv' into the folder 'data' of matchms installation. Replace it.<br>
Install the CDK python wrapper at https://github.com/hcji/pycdk

## USE

To use this programme:

1) Put your msp, json, csv or xml files into "INPUT/\<dedicated folder\>".
2) CSV files need to be separated by ";" and peaks columns need to be named "PEAK_LIST"
3.1) If you have HMDB xml spectrums, run HMDB_completion.bat first.
3.2) Else: Double click on "Run.bat".
4) when the execution is finished, retrieve your cleaned msp files into <br>"OUTPUT/CLEAN_MSP/(FINAL_POS or FINAL_NEG)"<br>and<br>csv version of the cleaned msp (POS and NEG) into<br>"OUTPUT/CSV/(FINAL_POS or FINAL_NEG)"
5) DO NOT FORGET TO DELETE FILES INTO "INPUT" and "OUTPUT" AFTER RETRIEVING CLEANED VERSIONS.