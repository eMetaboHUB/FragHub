import concurrent.futures
import pandas as pd
import numpy as np
import re
import os

def concatenate_clean_msp(clean_msp_path):
    CONCATENATE_LIST = []
    # append all spectrums of all cleaned files into CONCATENATE_LIST
    for files in os.listdir(clean_msp_path):
        if files.endswith("_clean.msp"):
            with open(os.path.join(clean_msp_path,files),"r",encoding="UTF-8") as buffer:
                temp = [element for element in buffer.read().split("\n\n") if element != "\n"]

            CONCATENATE_LIST.extend(temp)

    return CONCATENATE_LIST

def split_pos_neg(CONCATENATE_LIST):
    POS = []
    NEG = []
    for spectrum in CONCATENATE_LIST:
        SEARCH = re.search("(CHARGE): (.*)", spectrum)
        if SEARCH != None:
            if int(SEARCH.group(2)) > 0:
                POS.append(spectrum+"\n")
            else:
                NEG.append(spectrum+"\n")
    return POS, NEG

def harmonize_fields_names(spectrum):
    if spectrum is not None:
        expected_fields = ["SYNON","INCHIKEY","INSTRUMENT","FORMULA","SMILES","INCHI","COMMENT","IONIZATION","RESOLUTION","FRAGMENTATIONMODE","COMPOUNDNAME","SPECTRUMID","ADDUCT","MSLEVEL",
                           "INSTRUMENTTYPE","IONMODE","COLLISIONENERGY","PARENTMASS","PRECURSORMZ","CHARGE","NUM PEAKS","PREDICTED","RETENTIONTIME","FILENAME"]

        spectrum = re.sub("COMPOUND_NAME","COMPOUNDNAME",spectrum,flags=re.I)
        spectrum = re.sub("PRECURSOR_MZ", "PRECURSORMZ", spectrum, flags=re.I)
        spectrum = re.sub("INCHIKEY: \n", "INCHIKEY: None\n", spectrum)
        spectrum = re.sub("((^|\n)(.*?):) \n", "\n",spectrum)
        spectrum = re.sub("\n{2,}", "\n", spectrum)
        fields = re.finditer("(^|\n)(.*?):",spectrum)
        fields_names = [matche.group(2) for matche in fields]


        # Remove undesired punctuation
        for field in fields_names:
            spectrum = spectrum.replace(field,re.sub("\!|\(|\)|\[|\]|\{|\}|\;|\:|\'|\\|\,|\<|\>|\.|\/|\?|\@|\#|\$|\%|\^|\&|\*|\~|\+","",field))

        fields_names = [re.sub("\!|\(|\)|\[|\]|\{|\}|\;|\:|\'|\\|\,|\<|\>|\.|\/|\?|\@|\#|\$|\%|\^|\&|\*|\~|\+","",field) for field in fields_names]

        for field in fields_names:
            if field not in expected_fields: # Si champ dans le spectre pas voulu, on le supprime
                spectrum = re.sub(rf"{field}:.*\n","",spectrum)
        for field in expected_fields:
            if field not in fields_names: # Si un champ voulu est manquant, on le rajoute
                spectrum = field+": None\n"+spectrum

        return spectrum

def harmonize_adduct(spectrum):
    if spectrum != None:
        if re.search("((^|\n)(ADDUCT:)) (.*)\n",spectrum):
            adduct = re.search("((^|\n)(ADDUCT:)) (.*)\n",spectrum).group(4)
            if adduct != "None":
                if "[" not in adduct or "]" not in adduct:
                    if re.search("((^|\n)(IONMODE:)) (.*)\n",spectrum):
                        ionmode = re.search("((^|\n)(IONMODE:)) (.*)\n",spectrum).group(4)
                        if ionmode != "None":
                            if ionmode.lower().startswith("p"):
                                adduct = "["+adduct+"]"+"+\n"
                                spectrum = re.sub("(^|\n)(ADDUCT:) (.*)\n",f"\nADDUCT: {adduct}",spectrum)
                                return spectrum
                            elif ionmode.lower().startswith("n"):
                                adduct = "[" + adduct + "]" + "-\n"
                                spectrum = re.sub("((^|\n)(ADDUCT:)) (.*)\n", f"\nADDUCT: {adduct}", spectrum)
                                return spectrum
                        else:
                            return spectrum
                    else:
                        return spectrum
                elif re.search("((ADDUCT:)) ((.*\+)\*)\n",spectrum):
                    adduct = re.search("((ADDUCT:)) ((.*\+)\*)\n", spectrum).group(4)
                    spectrum = re.sub("((ADDUCT:)) ((.*\+)\*)\n",f"ADDUCT: {adduct}\n",spectrum)
                    return spectrum
                elif re.search("((ADDUCT:)) ((.*\-)\*)\n",spectrum):
                    adduct = re.search("((ADDUCT:)) ((.*\-)\*)\n", spectrum).group(4)
                    spectrum = re.sub("((ADDUCT:)) ((.*\-)\*)\n",f"ADDUCT: {adduct}\n",spectrum)
                    return spectrum
                else:
                    return spectrum
            else:
                return spectrum
        else:
            return spectrum
    else:
        return spectrum

def harmonize_retention_time(spectrum):
    if spectrum != None:
        if re.search("RETENTIONTIME: 0\n",spectrum):
            spectrum = re.sub("RETENTIONTIME: 0\n", "RETENTIONTIME: None\n", spectrum)
            return spectrum
        elif re.search("(^|\n)(RETENTIONTIME:) (.*)\n",spectrum):
            RT = re.search("((^|\n)(RETENTIONTIME:)) (.*)\n", spectrum).group(4)
            try:
                RT_test = float(RT)
                return spectrum
            except:
                spectrum = re.sub("((^|\n)(RETENTIONTIME:)) (.*)\n", "\nRETENTIONTIME: None\n", spectrum)
                return spectrum

def harmonize_ms_level(spectrum):
    if spectrum != None:
        spectrum = re.sub("MSLEVEL: MS\n", "MSLEVEL: 1\n", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: MS1","MSLEVEL: 1", spectrum,flags=re.I)
        spectrum = re.sub("MSLEVEL: MS2", "MSLEVEL: 2", spectrum,flags=re.I)
        spectrum = re.sub("MSLEVEL: MS3", "MSLEVEL: 3", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: MS4", "MSLEVEL: 4", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: 2-MS4 Composite", "MSLEVEL: 4", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: 2-MS5 Composite", "MSLEVEL: 5", spectrum, flags=re.I)

    return spectrum

def harmonize_collisionenergy(spectrum):
    if spectrum != None:
        if re.search("(^|\n)(COLLISIONENERGY:) ([0-9]*)((?: )?)((?:\()?)([a-zA-Z]*)((?:\))?)\n", spectrum):
            collisionenergy =  re.search("(^|\n)(COLLISIONENERGY:) ([0-9]*)((?: )?)((?:\()?)([a-zA-Z]*)((?:\))?)\n", spectrum)
            if collisionenergy.group(3) == '':
                return spectrum
            elif collisionenergy.group(3).isnumeric():
                if collisionenergy.group(6) == '':
                    return spectrum
                elif collisionenergy.group(6).isalpha():
                    fragmentation = collisionenergy.group(6)
                    # deleting alpha caracters
                    spectrum = re.sub("(^|\n)(COLLISIONENERGY:) ([0-9]*)((?: )?)((?:\()?)([a-zA-Z]*)((?:\))?)\n",f"\nCOLLISIONENERGY: {collisionenergy.group(3)}\n",spectrum)
                    # check if fragmentation mode already exist
                    if re.search("(^|\n)(FRAGMENTATIONMODE:) (.*)\n", spectrum):
                        fragmentationmode = re.search("(^|\n)(FRAGMENTATIONMODE:) (.*)\n", spectrum).group(3)
                        if fragmentationmode == "None":
                            spectrum = re.sub("FRAGMENTATIONMODE: None",f"FRAGMENTATIONMODE: {fragmentation}",spectrum)

    return spectrum

def harmonize_syns(spectrum):
    if spectrum != None:
        if re.search("\$:00in-source",spectrum):
            spectrum = re.sub("\$:00in-source","None",spectrum)

    return spectrum

def harmonize_formula(spectrum):
    if spectrum != None:
        if re.search("FORMULA: \[(.*)\](\+|-)\n",spectrum):
            FORMULA = re.search("FORMULA: \[(.*)\](\+|-)\n",spectrum).group(1)
            spectrum = re.sub("FORMULA: (.*)\n",f"FORMULA: {FORMULA}\n",spectrum)
        elif re.search("FORMULA: (.*)(\+|-)\n",spectrum):
            FORMULA = re.search("FORMULA: (.*)(\+|-)\n",spectrum).group(1)
            spectrum = re.sub("FORMULA: (.*)\n", f"FORMULA: {FORMULA}\n",spectrum)
        elif re.search(r"FORMULA: N\\A",spectrum,flags=re.I):
            spectrum = re.sub(r"FORMULA: N\\A", "FORMULA: None", spectrum, flags=re.I)

    return spectrum

def harmonize_empties(spectrum):
    if spectrum != None:
        if re.search(": \n",spectrum):
            spectrum = re.sub(": \n",": None\n",spectrum)
    return spectrum

def predicted_correction(spectrum):
    if spectrum != None:
        temp_spectrum = re.sub("FILENAME: (.*)\n", "", spectrum)
        temp_spectrum = re.sub("PREDICTED: (.*)\n", "", temp_spectrum)
        if re.search("in-silico|insilico|predicted",temp_spectrum,flags=re.I):
            spectrum = re.sub("PREDICTED: .*\n","PREDICTED: true\n",spectrum)
        else:
            spectrum = re.sub("PREDICTED: .*\n", "PREDICTED: false\n", spectrum)
    return spectrum

def remove_no_inchikey(spectrum):
    if spectrum != None:
        if re.search("INCHIKEY: \n|INCHIKEY: None\n",spectrum):
            return None
        else:
            return spectrum

def remove_no_mass(spectrum):
    if spectrum != None:
        if re.search("PRECURSORMZ: \n|PRECURSORMZ: None\n",spectrum):
            return None
        else:
            return spectrum

def harmonize_db_informations(spectrum):
    if spectrum != None:
        # GNPS CAS
        if re.search("(DB#=(.*)(;|,)((?: )?)(origin=(.*)))",spectrum):
            matche = re.search("(DB#=(.*)(;|,)((?: )?)(origin=(.*)))",spectrum)
            spectrum_id = matche.group(2)
            # adding spectrum id
            spectrum = re.sub("SPECTRUMID: None\n",f"SPECTRUMID: {spectrum_id}\n",spectrum)
        # MSMS CAS
        elif re.search("(spec_id=(.*)(;|,)((?: )?)(origin=(.*)))", spectrum):
            matche = re.search("(spec_id=(.*)(;|,)((?: )?)(origin=(.*)))", spectrum)
            spectrum_id = matche.group(2)
            # adding spectrum id
            spectrum = re.sub("SPECTRUMID: None\n", f"SPECTRUMID: {spectrum_id}\n", spectrum)

    return spectrum

def harmonize_fields_values(spectrum):
    spectrum = remove_no_inchikey(spectrum)
    spectrum = remove_no_mass(spectrum)
    spectrum = harmonize_adduct(spectrum)
    spectrum = harmonize_retention_time(spectrum)
    spectrum = harmonize_ms_level(spectrum)
    # spectrum = harmonize_collisionenergy(spectrum)
    spectrum = harmonize_syns(spectrum)
    spectrum = harmonize_formula(spectrum)
    spectrum = harmonize_empties(spectrum)
    spectrum = predicted_correction(spectrum)
    spectrum = harmonize_db_informations(spectrum)

    return spectrum

def correct_uncomplete_charge(msp_path):
    with open(msp_path, "r", encoding="UTF-8") as msp_buffer:
        content = msp_buffer.read()

    content = re.sub("charge: -\n","charge: -1\n",content,flags=re.I)

    with open(msp_path, "w", encoding="UTF-8") as msp_buffer:
        msp_buffer.write(content)

def msp_to_csv(clean_msp_path):
    # POS
    POS_DIR = os.path.join(clean_msp_path, "FINAL_POS")

    dictionary = {"PEAKS_LIST": []}

    first = True

    for files in os.listdir(POS_DIR):
        if files.endswith(".msp"):
            with open(os.path.join(POS_DIR, files), "r", encoding="UTF-8") as msp_file_buffer:
                msp_file = msp_file_buffer.read()

            spectrum_list = msp_file.split("\n\n")  # une liste de spectres

            for spectrum in spectrum_list:
                if spectrum != "\n":
                    fields = re.findall(r"(.+?):(.*)\n", spectrum)
                    if first == True:
                        for element in fields:
                            dictionary[element[0]] = []
                        first = False

                    if first == False:
                        for element in fields:
                            dictionary[element[0]].append(element[1])

                    dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2).split("\n"))

    POS_df = pd.DataFrame.from_dict(dictionary)
    POS_df.to_csv("../OUTPUT/CSV/FINAL_POS/POS_clean.csv", sep=";", encoding="UTF-8", index=False)

    # NEG
    NEG_DIR = os.path.join(clean_msp_path, "FINAL_NEG")

    dictionary = {"PEAKS_LIST": []}

    first = True

    for files in os.listdir(NEG_DIR):
        if files.endswith(".msp"):
            with open(os.path.join(NEG_DIR, files), "r", encoding="UTF-8") as msp_file_buffer:
                msp_file = msp_file_buffer.read()

            spectrum_list = msp_file.split("\n\n")  # une liste de spectres

            for spectrum in spectrum_list:
                if spectrum != "\n":
                    fields = re.findall(r"(.+?):(.*)", spectrum)
                    if first == True:
                        for element in fields:
                            dictionary[element[0]] = []
                        first = False

                    if first == False:
                        for element in fields:
                            dictionary[element[0]].append(element[1])

                    dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2).split("\n"))

    NEG_df = pd.DataFrame.from_dict(dictionary)
    NEG_df.to_csv("../OUTPUT/CSV/FINAL_NEG/NEG_clean.csv", sep=";", encoding="UTF-8", index=False)

def find_indices(l, value):
    duplicatas_index = [index for index, item in enumerate(l) if item == value]
    if duplicatas_index != None:
        if len(duplicatas_index) >= 2:
            duplicatas_index = duplicatas_index[1:]
            return duplicatas_index
        else:
            return []
    else:
        return []

def remove_duplicatas(POS, NEG):
    # POS
    POS_list = []
    POS_index_to_delete = []
    for spectrum in POS: # Générer une liste donnée réduite à inchikey + peak_list
        POS_list.append(str({"INCHIKEY": re.search("INCHIKEY: (.*)\n",spectrum).group(1), "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)}))

    for current_dict in POS_list:
        POS_index_to_delete.extend(find_indices(POS_list,current_dict))

    POS_FILTERED = []
    compteur = 0
    for spectrum in POS:
        if compteur not in POS_index_to_delete:
            POS_FILTERED.append(spectrum)
        compteur += 1

    # NEG
    NEG_list = []
    NEG_index_to_delete = []
    for spectrum in NEG:  # Générer une liste donnée réduite à inchikey + peak_list
        NEG_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                         "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2)})

    for current_dict in NEG_list:
        NEG_index_to_delete.extend(find_indices(NEG_list, current_dict))

    NEG_FILTERED = []
    compteur = 0
    for spectrum in NEG:
        if compteur not in NEG_index_to_delete:
            NEG_FILTERED.append(spectrum)
        compteur += 1

    return POS_FILTERED,NEG_FILTERED











