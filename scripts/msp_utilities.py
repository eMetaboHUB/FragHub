import pubchempy as pcp
from rdkit import Chem
from tqdm import tqdm
import uuid
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

def generate_dict_inchikey_smiles_inchi(CONCATENATE_LIST):
    # généré un dictionnaire des inchikey,smiles,inchi des DB
    INCHIKEY = ""
    INCHI = ""
    SMILES = ""

    INCHIKEY_DICT = {}
    sub_dict = {"INCHI": "", "SMILES": ""}

    for spectrum in CONCATENATE_LIST:
        if re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n",spectrum): # on à l'inchi
            INCHIKEY = re.search("INCHIKEY: ([A-Z]{14}-[A-Z]{10}-N)\n",spectrum).group(1)
            # Recherchez la substance correspondant à l'InChIKey dans PubChem
            if INCHIKEY not in INCHIKEY_DICT.keys():
                try:
                    compound = pcp.get_compounds(INCHIKEY, 'inchikey')[0]
                    INCHI = compound.inchi
                    SMILES = compound.canonical_smiles

                    INCHIKEY_DICT[INCHIKEY] = sub_dict
                    INCHIKEY_DICT[INCHIKEY]["INCHI"] = INCHI
                    INCHIKEY_DICT[INCHIKEY]["SMILES"] = SMILES
                except IndexError:
                    continue

        elif re.search("INCHI: (.*)\n",spectrum): # on à l'inchi
            if "=" in re.search("INCHI: (.*)\n",spectrum).group(1):
                INCHI = re.search("INCHI: (.*)\n",spectrum).group(1)
                INCHIKEY = Chem.MolToInchiKey(Chem.MolFromInchi(INCHI))
                SMILES = Chem.MolToSmiles(Chem.MolFromInchi(INCHI))

                if INCHIKEY not in INCHIKEY_DICT.keys():
                    INCHIKEY_DICT[INCHIKEY] = sub_dict
                    INCHIKEY_DICT[INCHIKEY]["INCHI"] = INCHI
                    INCHIKEY_DICT[INCHIKEY]["SMILES"] = SMILES

        elif re.search("SMILES: (.*)\n",spectrum):
            if "C" in re.search("SMILES: (.*)\n", spectrum).group(1):
                SMILES = re.search("SMILES: (.*)\n", spectrum).group(1)
                INCHI = Chem.MolToInchi(Chem.MolFromSmiles(SMILES))
                INCHIKEY = Chem.MolToInchiKey(Chem.MolFromInchi(INCHI))

                if INCHIKEY not in INCHIKEY_DICT.keys():
                    INCHIKEY_DICT[INCHIKEY] = sub_dict
                    INCHIKEY_DICT[INCHIKEY]["INCHI"] = INCHI
                    INCHIKEY_DICT[INCHIKEY]["SMILES"] = SMILES

    return INCHIKEY_DICT




def mols_derivator(CONCATENATE_LIST):
    INCHIKEY_DICT = generate_dict_inchikey_smiles_inchi(CONCATENATE_LIST)

    # temp_list = []
    # for spectrums in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green",desc="\t  processing"):



