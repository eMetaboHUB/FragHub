import concurrent.futures
import pubchempy as pcp
from rdkit import Chem
from tqdm import tqdm
import threading
import uuid
import re
import os

INCHI_DICT_LOCK = threading.Lock()

def concatenate_clean_msp(clean_msp_path):
    CONCATENATE_LIST = []
    # append all spectrums of all cleaned files into CONCATENATE_LIST
    for files in os.listdir(clean_msp_path):
        if files.endswith("_clean.msp"):
            with open(os.path.join(clean_msp_path,files),"r",encoding="UTF-8") as buffer:
                temp = [element for element in buffer.read().split("\n\n") if element != "\n"]

            CONCATENATE_LIST.extend(temp)

    return CONCATENATE_LIST

def correct_uncomplete_charge(msp_path):
    with open(msp_path, "r", encoding="UTF-8") as msp_buffer:
        content = msp_buffer.read()

    content = re.sub("charge: -\n","charge: -1\n",content,flags=re.I)

    with open(msp_path, "w", encoding="UTF-8") as msp_buffer:
        msp_buffer.write(content)

def names_completion(CONCATENATE_LIST):
    inchi_names = {}
    inchi = "None"
    name = "None"

    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
        if re.search("INCHI: (.*)\n", spectrum):
            inchi = re.search("INCHI: (.*)\n", spectrum).group(1)
        if re.search("\nNAME: (.*)\n", spectrum):
            name = re.search("\nNAME: (.*)\n", spectrum).group(1)

        # If InChIKey and name are valid, add them to the dictionary
        if inchi and name != "None":
            inchi_names[inchi] = name

    # Update missing names with corresponding names in dictionary list
    updated_spetcrum_list = []
    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t    updating"):
        if re.search("INCHI: (.*)\n", spectrum):
            inchi = re.search("INCHI: (.*)\n", spectrum).group(1)
        if re.search("\nNAME: (.*)\n", spectrum):
            name = re.search("\nNAME: (.*)\n", spectrum).group(1)
        if name == "None":
            if inchi in inchi_names.keys():
                spectrum = re.sub("\nNAME: (.*)\n", f"\nNAME: {inchi_names[inchi]}\n", spectrum)

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
        if not (re.search("SMILES: None\n",spectrums) or re.search("INCHI: None\n",spectrums)):
            CONCATENATE_LIST_temp.append(spectrums)

    return CONCATENATE_LIST_temp

def unique_id_generator(CONCATENATE_LIST):
    temp_list = []
    for spectrums in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
        if re.search("(PREDICTED:(.*)\n)",spectrums):
            predicted_line = re.search("(PREDICTED:(.*)\n)",spectrums).group(1)
            FRAGBANKID = "FRAGBANKID: "+str(uuid.uuid4())+"\n"
            temp_list.append(re.sub("(PREDICTED:(.*)\n)",rf"{predicted_line}{FRAGBANKID}",spectrums))

    return temp_list

def cas_1(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" in re.search("INCHI: (.*)\n", spectrum).group(1):
            if not re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" not in re.search("SMILES: (.*)\n",spectrum).group(1)):
                cas = True

    return cas

def cas_2(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" in re.search("INCHI: (.*)\n", spectrum).group(1):
            if not re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def cas_3(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" in re.search("INCHI: (.*)\n", spectrum).group(1):
            if re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" not in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def cas_4(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" not in re.search("INCHI: (.*)\n", spectrum).group(1):
            if not re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ( "C" in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def cas_5(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" not in re.search("INCHI: (.*)\n", spectrum).group(1):
            if re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" not in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def cas_6(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" not in re.search("INCHI: (.*)\n", spectrum).group(1):
            if re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def cas_7(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" in re.search("INCHI: (.*)\n", spectrum).group(1):
            if re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def cas_8(spectrum):
    cas = False
    if re.search("INCHI: (.*)\n", spectrum):
        if "=" not in re.search("INCHI: (.*)\n", spectrum).group(1):
            if not re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum) and ("C" not in re.search("SMILES: (.*)\n", spectrum).group(1)):
                cas = True

    return cas

def spectrum_have_a_name(spectrum):
    name = True
    if re.search("\nNAME: None\n",spectrum):
        name = False

    return name

def modif_INCHI_DICT(spectrum,INCHIKEY,INCHI,SMILES):
    with INCHI_DICT_LOCK:
        global INCHI_DICT

        if INCHI not in INCHI_DICT.keys():
            sub_dict = {"INCHIKEY": "", "SMILES": "", "NAME": ""}
            INCHI_DICT[INCHI] = sub_dict
            INCHI_DICT[INCHI]["INCHIKEY"] = INCHIKEY
            INCHI_DICT[INCHI]["SMILES"] = SMILES
            if spectrum_have_a_name(spectrum):
                INCHI_DICT[INCHI]["NAME"] = re.search("\nNAME: (.*)\n", spectrum).group(1)

def generate_dict_inchikey_smiles_inchi(spectrum):
    # généré un dictionnaire des inchikey,smiles,inchi des DB

    INCHIKEY = ""
    INCHI = ""
    SMILES = ""

    if cas_1(spectrum) or cas_2(spectrum) or cas_3(spectrum) or cas_7(spectrum):
        INCHI = re.search("INCHI: (.*)\n", spectrum).group(1)
        INCHI = re.sub("inchi=","InChI=",INCHI,flags=re.IGNORECASE)
        try:
            INCHIKEY = Chem.MolToInchiKey(Chem.MolFromInchi(INCHI)) # sanitize=True, removeHs=True
            SMILES = Chem.MolToSmiles(Chem.MolFromInchi(INCHI))
        except: # Si InChI non parsé par RDkit, on essaye avec le SMILE:  !!! PAS TOTALEMENT FIABLE !!!
            if cas_2(spectrum) or cas_7(spectrum):
                try:
                    SMILES = re.search("SMILES: (.*)\n", spectrum).group(1)

                    INCHI = Chem.MolToInchi(Chem.MolFromSmiles(SMILES))
                    INCHIKEY = Chem.MolToInchiKey(Chem.MolFromInchi(INCHI))
                except:
                    return None
            else:
                return None

        modif_INCHI_DICT(spectrum,INCHIKEY,INCHI,SMILES)

    elif cas_5(spectrum):
        try:
            if re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum):
                INCHIKEY = re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrum).group(1)

                sub_dict = {"INCHIKEY": "", "SMILES": "","NAME": ""}
                compound = pcp.get_compounds(INCHIKEY, 'inchikey')[0]
                INCHI = compound.inchi # ATTENTION !!! Deux InChiKey identiques peuvent avoir des InChi différents (Bien que cela soit rare).
                SMILES = compound.canonical_smiles

                modif_INCHI_DICT(spectrum,INCHIKEY,INCHI,SMILES)

        except IndexError:
            return None

    elif cas_4(spectrum) or cas_6(spectrum):
        SMILES = re.search("SMILES: (.*)\n", spectrum).group(1)
        INCHI = None
        try:
            INCHI = Chem.MolToInchi(Chem.MolFromSmiles(SMILES))
            INCHIKEY = Chem.MolToInchiKey(Chem.MolFromInchi(INCHI))
        except: # SMILE:  !!! PAS TOTALEMENT FIABLE !!!
            return None

        modif_INCHI_DICT(spectrum,INCHIKEY,INCHI,SMILES)

    elif cas_8(spectrum): # Aucune smiles, inchi ou inchikey ==> on essaye avec le nom sur pubchem
        if spectrum_have_a_name(spectrum):
            try:
                NAME = re.search("\nNAME: (.*)\n",spectrum).group(1)

                compound = pcp.get_compounds(NAME, 'name')[0]
                INCHI = compound.inchi  # ATTENTION !!! Deux InChiKey identiques peuvent avoir des InChi différents (Bien que cela soit rare).
                INCHIKEY = compound.inchikey
                SMILES = compound.canonical_smiles

                modif_INCHI_DICT(spectrum,INCHIKEY,INCHI,SMILES)

            except IndexError:
                return None


def generate_dict_inchikey_smiles_inchi_processing(CONCATENATE_LIST):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        tqdm(executor.map(generate_dict_inchikey_smiles_inchi, CONCATENATE_LIST), total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  generating")



def mols_derivator(spectrums):
    global INCHI_DICT

    # for key,value in INCHIKEY_DICT.items():
    #     print(key, " ",value)

    temp_list = []
    if cas_1(spectrums) or cas_2(spectrums) or cas_3(spectrums) or cas_7(spectrums):
        INCHI = re.search("INCHI: (.*)\n", spectrums).group(1)

        spectrums = re.sub("INCHIKEY: (.*)\n", rf"INCHIKEY: {INCHI_DICT[INCHI]['INCHIKEY']}\n",spectrums)
        spectrums = re.sub("SMILES: (.*)\n", rf"SMILES: {INCHI_DICT[INCHI]['SMILES']}\n", spectrums)

        return spectrums

    elif cas_5(spectrums): # ATTENTION !!! Deux InChiKey identiques peuvent avoir des InChi différents (Bien que cela soit rare).
        INCHIKEY = re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n", spectrums).group(1)
        # retrouver l'inchi
        INCHI = None
        for key in INCHI_DICT.keys():
            if INCHI_DICT[key]["INCHIKEY"] == INCHIKEY:
                INCHI = key
                break

        if INCHI != None:
            spectrums = re.sub("INCHI: (.*)\n", rf"INCHI: {INCHI}\n", spectrums)
            spectrums = re.sub("SMILES: (.*)\n", rf"SMILES: {INCHI_DICT[INCHI]['SMILES']}\n", spectrums)

            return spectrums

    elif cas_4(spectrums) or cas_6(spectrums): # SMILE:  !!! PAS TOTALEMENT FIABLE !!!
        SMILES = re.search("SMILES: (.*)\n", spectrums).group(1)
        # retrouver l'inchi
        INCHI = None
        for key in INCHI_DICT.keys():
            if INCHI_DICT[key]["SMILES"] == SMILES:
                INCHI = key
                break

        if INCHI != None:
            spectrums = re.sub("INCHI: (.*)\n", rf"INCHI: {INCHI}\n", spectrums)
            spectrums = re.sub("INCHIKEY: (.*)\n", rf"INCHIKEY: {INCHI_DICT[INCHI]['INCHIKEY']}\n", spectrums)

            return spectrums

    elif cas_8(spectrums):
        if spectrum_have_a_name(spectrums):
            NAME = re.search("\nNAME: (.*)\n",spectrums).group(1)
            # retrouver l'inchi
            INCHI = None
            for key in INCHI_DICT.keys():
                if INCHI_DICT[key]["NAME"] == NAME:
                    INCHI = key
                    break

            if INCHI != None:
                spectrums = re.sub("INCHI: (.*)\n", rf"INCHI: {INCHI}\n", spectrums)
                spectrums = re.sub("INCHIKEY: (.*)\n", rf"INCHIKEY: {INCHI_DICT[INCHI]['INCHIKEY']}\n", spectrums)
                spectrums = re.sub("SMILES: (.*)\n", rf"SMILES: {INCHI_DICT[INCHI]['SMILES']}\n", spectrums)

                return spectrums

def mols_derivator_processing(CONCATENATE_LIST):
    global INCHI_DICT
    INCHI_DICT = {}

    generate_dict_inchikey_smiles_inchi_processing(CONCATENATE_LIST)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(mols_derivator, CONCATENATE_LIST), total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


