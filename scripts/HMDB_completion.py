from tqdm import tqdm
import pandas as pd
import os
import re

HMDB_df = pd.read_csv("../datas/HMDB.csv",sep=";",encoding="UTF-8")

dirs = os.listdir("../INPUT/XML")

print("-- COMPLETING HMDB SPECTRUMS --")

for files in tqdm(dirs, total=len(dirs), unit="files", colour="green"):
    if files.endswith(".xml"):
        with open(os.path.join("../INPUT/XML", files), "r", encoding="UTF-8") as buffer:
            xml_content = buffer.read()

        HMDB_ID = re.search("<database-id>(.*)</database-id>",xml_content).group(1) # hmdb id retrieval for csv matching
        if HMDB_ID != None:
            line = HMDB_df.loc[HMDB_df['EXTERNAL_ID'] == HMDB_ID]

            if not line.empty:
                INCHIKEY = line["INCHIKEY"].values[0]
                SMILES = line["SMILES"].values[0]
                INCHI = line["INCHI"].values[0]
                MOLECULAR_FORMULA = line["MOLECULAR_FORMULA"].values[0]
                ACC_MASS = line["ACC_MASS"].values[0]
            else:
                continue

            if re.search("<sample-mass>(.*)</sample-mass>",xml_content):
                old_tag = "</sample-mass>\n"
                new_tag = f"</sample-mass>\n  <inchikey>{INCHIKEY}</inchikey>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "</sample-mass>\n"
                new_tag = f"</sample-mass>\n  <SMILES>{SMILES}</SMILES>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "</sample-mass>\n"
                new_tag = f"</sample-mass>\n  <INCHI>{INCHI}</INCHI>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "</sample-mass>\n"
                new_tag = f"</sample-mass>\n  <MOLECULAR_FORMULA>{MOLECULAR_FORMULA}</MOLECULAR_FORMULA>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "</sample-mass>\n"
                new_tag = f"</sample-mass>\n  <ACC_MASS>{ACC_MASS}</ACC_MASS>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

            elif re.search("<sample-mass nil=\"true\"/>",xml_content):
                old_tag = "<sample-mass nil=\"true\"/>\n"
                new_tag = f"<sample-mass nil=\"true\"/>\n  <inchikey>{INCHIKEY}</inchikey>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "<sample-mass nil=\"true\"/>\n"
                new_tag = f"<sample-mass nil=\"true\"/>\n  <SMILES>{SMILES}</SMILES>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "<sample-mass nil=\"true\"/>\n"
                new_tag = f"<sample-mass nil=\"true\"/>\n  <INCHI>{INCHI}</INCHI>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "<sample-mass nil=\"true\"/>\n"
                new_tag = f"<sample-mass nil=\"true\"/>\n  <MOLECULAR_FORMULA>{MOLECULAR_FORMULA}</MOLECULAR_FORMULA>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

                old_tag = "<sample-mass nil=\"true\"/>\n"
                new_tag = f"<sample-mass nil=\"true\"/>\n  <ACC_MASS>{ACC_MASS}</ACC_MASS>\n"
                xml_content = xml_content.replace(old_tag, new_tag)

            with open(os.path.join("../INPUT/XML",files),"w",encoding="UTF-8") as buffer:
                buffer.write(xml_content)
        else:
            continue

