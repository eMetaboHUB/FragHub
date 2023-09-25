from tqdm import tqdm
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

def harmonize_fields_names(spectrum):
    if spectrum is not None:
        expected_fields = ["SYNON","INCHIKEY","INSTRUMENT","FORMULA","SMILES","INCHI","COMMENT","IONIZATION","RESOLUTION","FRAGMENTATIONMODE","NAME","SPECTRUMID","PRECURSORTYPE","MSLEVEL",
                           "INSTRUMENTTYPE","IONMODE","COLLISIONENERGY","PARENTMASS","PRECURSORMZ","CHARGE","NUM PEAKS","PREDICTED","RETENTIONTIME","FILENAME"]

        spectrum = re.sub("COMPOUND_NAME:","NAME:",spectrum,flags=re.I)
        spectrum = re.sub("PRECURSOR_MZ:", "PRECURSORMZ:", spectrum, flags=re.I)
        spectrum = re.sub("ADDUCT:", "PRECURSORTYPE:", spectrum, flags=re.I)
        spectrum = re.sub("INCHIKEY: \n", "INCHIKEY: None\n", spectrum)
        spectrum = re.sub("((^|\n)(.*?):) \n", "\n",spectrum)
        spectrum = re.sub("\n{2,}", "\n", spectrum)
        fields = re.finditer("(^|\n)(.*?):",spectrum)
        fields_names = [matche.group(2) for matche in fields]


        # Remove undesired punctuation
        for field in fields_names:
            spectrum = spectrum.replace(field,re.sub("\!|\(|\)|\[|\]|\{|\}|\;|\:|\'|\\|\,|\<|\>|\.|\/|\?|\@|\#|\$|\%|\^|\&|\*|\~|\+","",field))

        fields_names = [re.sub("\!|\(|\)|\[|\]|\{|\}|\;|\'|\\|\,|\<|\>|\.|\/|\?|\@|\#|\$|\%|\^|\&|\*|\~|\+","",field) for field in fields_names]

        for field in fields_names:
            if field not in expected_fields: # Si champ dans le spectre pas voulu, on le supprime
                spectrum = re.sub(rf"(^|\n){field}:.*\n","\n",spectrum)
        for field in expected_fields:
            if field not in fields_names: # Si un champ voulu est manquant, on le rajoute
                spectrum = field+": None\n"+spectrum

        # Sort fields
        SPECTRUM = ""

        SPECTRUM = SPECTRUM + re.search("((^|\n)(FILENAME: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PREDICTED: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(SPECTRUMID: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(RESOLUTION: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(SYNON: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(CHARGE: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PARENTMASS: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(IONIZATION: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(MSLEVEL: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(FRAGMENTATIONMODE: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(NAME: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PRECURSORMZ: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PRECURSORTYPE: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INSTRUMENTTYPE: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INSTRUMENT: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(SMILES: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INCHI: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INCHIKEY: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(COLLISIONENERGY: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(FORMULA: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(RETENTIONTIME: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(IONMODE: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(COMMENT: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(NUM PEAKS: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)

        return SPECTRUM

def harmonize_adduct(spectrum):
    if spectrum != None:
        if re.search("((^|\n)(PRECURSORTYPE:)) (.*)\n",spectrum):
            adduct = re.search("((^|\n)(PRECURSORTYPE:)) (.*)\n",spectrum).group(4)
            if adduct != "None":
                if "[" not in adduct or "]" not in adduct:
                    if re.search("((^|\n)(IONMODE:)) (.*)\n",spectrum):
                        ionmode = re.search("((^|\n)(IONMODE:)) (.*)\n",spectrum).group(4)
                        if ionmode != "None":
                            if ionmode.lower().startswith("p"):
                                adduct = "["+adduct+"]"+"+\n"
                                spectrum = re.sub("(^|\n)(PRECURSORTYPE:) (.*)\n",f"\nPRECURSORTYPE: {adduct}",spectrum)
                                return spectrum
                            elif ionmode.lower().startswith("n"):
                                adduct = "[" + adduct + "]" + "-\n"
                                spectrum = re.sub("((^|\n)(PRECURSORTYPE:)) (.*)\n", f"\nPRECURSORTYPE: {adduct}", spectrum)
                                return spectrum
                        else:
                            return spectrum
                    else:
                        return spectrum
                elif re.search("((PRECURSORTYPE:)) ((.*\+)\*)\n",spectrum):
                    adduct = re.search("((PRECURSORTYPE:)) ((.*\+)\*)\n", spectrum).group(4)
                    spectrum = re.sub("((PRECURSORTYPE:)) ((.*\+)\*)\n",f"PRECURSORTYPE: {adduct}\n",spectrum)
                    return spectrum
                elif re.search("((PRECURSORTYPE:)) ((.*\-)\*)\n",spectrum):
                    adduct = re.search("((PRECURSORTYPE:)) ((.*\-)\*)\n", spectrum).group(4)
                    spectrum = re.sub("((PRECURSORTYPE:)) ((.*\-)\*)\n",f"PRECURSORTYPE: {adduct}\n",spectrum)
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
        if not re.search("FILENAME: MSMS_Public.*",spectrum):
            temp_spectrum = re.sub("PREDICTED: (.*)\n", "", spectrum)
            if re.search("in-silico|insilico|predicted|theoretical|Annotation level-3",temp_spectrum,flags=re.I):
                spectrum = re.sub("PREDICTED: .*\n","PREDICTED: true\n",spectrum)
            else:
                spectrum = re.sub("PREDICTED: .*\n", "PREDICTED: false\n", spectrum)
        else:
            temp_spectrum = re.sub("FILENAME: (.*)\n", "", spectrum)
            temp_spectrum = re.sub("PREDICTED: (.*)\n", "", temp_spectrum)
            if re.search("in-silico|insilico|predicted|theoretical|Annotation level-3", temp_spectrum, flags=re.I):
                spectrum = re.sub("PREDICTED: .*\n", "PREDICTED: true\n", spectrum)
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

def correct_ionmode(spectrum):
    if spectrum != None:
        if re.search("IONMODE: n/a",spectrum):
            if re.search("PRECURSORTYPE: (.*)\-\n",spectrum):
                spectrum = re.sub("IONMODE: n/a","IONMODE: negative",spectrum)
                spectrum = re.sub("CHARGE: (.*)\n", "CHARGE: -1\n", spectrum)
                return spectrum
            elif re.search("PRECURSORTYPE: (.*)\+\n",spectrum):
                spectrum = re.sub("IONMODE: n/a","IONMODE: positive",spectrum)
                spectrum = re.sub("CHARGE: (.*)\n", "CHARGE: 1\n", spectrum)
                return spectrum
        else:
            return spectrum

def harmonize_fields_values(spectrum):
    spectrum = remove_no_inchikey(spectrum)
    spectrum = remove_no_mass(spectrum)
    spectrum = harmonize_adduct(spectrum)
    spectrum = correct_ionmode(spectrum)
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

def names_completion(CONCATENATE_LIST):
    inchikey_names = {}
    inchikey = "None"
    name = "None"

    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        if re.search("\nNAME: (.*)\n", spectrum):
            name = re.search("\nNAME: (.*)\n", spectrum).group(1)

        # If InChIKey and name are valid, add them to the dictionary
        if inchikey and name != "None":
            inchikey_names[inchikey] = name

    # Update missing names with corresponding names in dictionary list
    updated_spetcrum_list = []
    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t    updating"):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        if re.search("\nNAME: (.*)\n", spectrum):
            name = re.search("\nNAME: (.*)\n", spectrum).group(1)
        if name == "None":
            if inchikey in inchikey_names.keys():
                spectrum = re.sub("\nNAME: (.*)\n", f"\nNAME: {inchikey_names[inchikey]}\n", spectrum)

        updated_spetcrum_list.append(spectrum)

    return updated_spetcrum_list

def inchi_smiles_completion(CONCATENATE_LIST):
    inchikey_inchi = {}
    inchikey_smiles = {}

    inchikey = "None"
    inchi = "None"
    smiles = "None"

    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        if re.search("\nINCHI: (.*)\n", spectrum):
            inchi = re.search("\nINCHI: (.*)\n", spectrum).group(1)
        if re.search("\nSMILES: (.*)\n", spectrum):
            smiles = re.search("\nSMILES: (.*)\n", spectrum).group(1)

        if inchikey and inchi != "None":
            inchikey_inchi[inchikey] = inchi

        if inchikey and smiles != "None":
            inchikey_smiles[inchikey] = smiles

    # Update missing inchi/smiles with corresponding inchikey in dictionary list
    updated_spetcrum_list = []
    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t    updating"):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        # INCHI
        if re.search("\nINCHI: (.*)\n", spectrum):
            inchi = re.search("\nINCHI: (.*)\n", spectrum).group(1)
        if inchi == "None":
            if inchikey in inchikey_inchi.keys():
                spectrum = re.sub("\nINCHI: (.*)\n", rf"\nINCHI: {re.escape(inchikey_inchi[inchikey])}\n", spectrum)
                spectrum = re.sub(r"\\\\",r"\\",spectrum)
        # SMILES
        if re.search("\nSMILES: (.*)\n", spectrum):
            smiles = re.search("\nSMILES: (.*)\n", spectrum).group(1)
        if smiles == "None":
            if inchikey in inchikey_smiles.keys():
                spectrum = re.sub("\nSMILES: (.*)\n", rf"\nSMILES: {re.escape(inchikey_smiles[inchikey])}\n", spectrum)
                spectrum = re.sub(r"\\\\", r"\\", spectrum)

        updated_spetcrum_list.append(spectrum)

    return updated_spetcrum_list

def remove_no_smiles_inchi(CONCATENATE_LIST):
    CONCATENATE_LIST_temp = []
    for spectrums in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
        if not re.search("SMILES: None\n",spectrums) or not re.search("INCHI: None\n",spectrums):
            CONCATENATE_LIST_temp.append(spectrums)

    return CONCATENATE_LIST_temp
