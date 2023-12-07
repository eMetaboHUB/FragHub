import pandas as pd
import os
import re

def in_silico(spectrums):
    if re.search("in-silico|insilico| predicted |predicted: true|theoretical|Annotation level-3",spectrums, flags=re.I):
        return True

def is_GC_MS(spectrums):
    if re.search(" GC-",spectrums) or re.search(" EI-",spectrums):
        return True

input_path = r"C:\Users\Axel\Documents\PYTHON\MSP_V3\ORIGINALS_msp_DB\MSP\Sept2023"

dictionary = {"FILENAME": [], "INCHIKEY": []}
df = pd.DataFrame()

for files in os.listdir(input_path):
    if files.endswith(".msp"):
        if files != "MoNA-export-In-Silico_Spectra.msp":

            with open(os.path.join(input_path,files),"r",encoding="UTF-8") as buffer:
                content = buffer.read()

            content = content.split("\n\n")

            for spectrums in content:
                if not in_silico(spectrums):
                    if not is_GC_MS(spectrums):
                        if re.search("([A-Z]{14}-[A-Z]{10}-N)",spectrums):
                            dictionary["FILENAME"].append(files.replace(".msp",""))
                            dictionary["INCHIKEY"].append(re.search("([A-Z]{14}-[A-Z]{10}-N)",spectrums).group(1))

final_df = df.from_dict(dictionary)

final_df.to_csv(r"C:\Users\Axel\Documents\KNIME\MSP3_check\original_db.csv", sep = ";", encoding="UTF-8", index=False)




