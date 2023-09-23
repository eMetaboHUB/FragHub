from tqdm import tqdm
import pandas as pd
import time
import re

def remove_dupli_POS_LC(POS_LC):
    dictionary = {}
    POS_LC_df = pd.DataFrame()
    first = True
    empty = False
    if len(POS_LC) == 0:
        empty = True

    for spectrum in POS_LC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        POS_LC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(POS_LC_df), desc="\t\t  POS_LC", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        POS_LC_df = POS_LC_df.loc[~POS_LC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(POS_LC_df))

        # Fermer la barre de progression
        t.close()

    return  POS_LC_df

def remove_dupli_POS_LC_In_Silico(POS_LC_In_Silico):
    dictionary = {}
    POS_LC_df = pd.DataFrame()
    first = True
    empty = False
    if len(POS_LC_In_Silico) == 0:
        empty = True

    for spectrum in POS_LC_In_Silico:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        POS_LC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(POS_LC_df), desc="POS_LC_In_Silico", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        POS_LC_df = POS_LC_df.loc[~POS_LC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(POS_LC_df))

        # Fermer la barre de progression
        t.close()

    return POS_LC_df

def remove_dupli_POS_GC(POS_GC):
    dictionary = {}
    POS_GC_df = pd.DataFrame()
    first = True
    empty = False
    if len(POS_GC) == 0:
        empty = True

    for spectrum in POS_GC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        POS_GC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(POS_GC_df), desc="\t\t  POS_GC", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        POS_GC_df = POS_GC_df.loc[~POS_GC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(POS_GC_df))

        # Fermer la barre de progression
        t.close()

    return POS_GC_df

def remove_dupli_POS_GC_In_Silico(POS_GC_In_Silico):
    dictionary = {}
    POS_GC_df = pd.DataFrame()
    first = True
    empty = False
    if len(POS_GC_In_Silico) == 0:
        empty = True

    for spectrum in POS_GC_In_Silico:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        POS_GC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(POS_GC_df), desc="POS_GC_In_Silico", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        POS_GC_df = POS_GC_df.loc[~POS_GC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(POS_GC_df))

        # Fermer la barre de progression
        t.close()

    return POS_GC_df

def remove_dupli_NEG_LC(NEG_LC):
    dictionary = {}
    NEG_LC_df = pd.DataFrame()
    first = True
    empty = False
    if len(NEG_LC) == 0:
        empty = True

    for spectrum in NEG_LC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        NEG_LC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(NEG_LC_df), desc="\t\t  NEG_LC", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        NEG_LC_df = NEG_LC_df.loc[~NEG_LC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(NEG_LC_df))

        # Fermer la barre de progression
        t.close()

    return NEG_LC_df

def remove_dupli_NEG_LC_In_Silico(NEG_LC_In_Silico):
    dictionary = {}
    NEG_LC_df = pd.DataFrame()
    first = True
    empty = False
    if len(NEG_LC_In_Silico) == 0:
        empty = True

    for spectrum in NEG_LC_In_Silico:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        NEG_LC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(NEG_LC_df), desc="NEG_LC_In_Silico", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        NEG_LC_df = NEG_LC_df.loc[~NEG_LC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(NEG_LC_df))

        # Fermer la barre de progression
        t.close()

    return NEG_LC_df


def remove_dupli_NEG_GC(NEG_GC):
    dictionary = {}
    NEG_GC_df = pd.DataFrame()
    first = True
    empty = False
    if len(NEG_GC) == 0:
        empty = True

    for spectrum in NEG_GC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        NEG_GC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(NEG_GC_df), desc="\t\t  NEG_GC", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        NEG_GC_df = NEG_GC_df.loc[~NEG_GC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(NEG_GC_df))

        # Fermer la barre de progression
        t.close()

    return NEG_GC_df

def remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico):
    dictionary = {}
    NEG_GC_df = pd.DataFrame()
    first = True
    empty = False
    if len(NEG_GC_In_Silico) == 0:
        empty = True

    for spectrum in NEG_GC_In_Silico:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        NEG_GC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        # Créer la barre de progression
        t = tqdm(total=len(NEG_GC_df), desc="NEG_GC_In_Silico", colour="green", unit=" row")

        # Supprimer les doublons et mettre à jour la barre de progression
        NEG_GC_df = NEG_GC_df.loc[~NEG_GC_df.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
        t.update(len(NEG_GC_df))

        # Fermer la barre de progression
        t.close()

    return NEG_GC_df


def re_write_MSP_POS_LC(POS_LC_df):
    POS_LC = []
    for index,row in tqdm(POS_LC_df.iterrows(), total=len(POS_LC_df), desc="\t\t  POS_LC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_LC.append(SPECTRUM)

    return POS_LC

def re_write_MSP_POS_LC_In_Silico(POS_LC_df_insilico):
    POS_LC = []
    for index, row in tqdm(POS_LC_df_insilico.iterrows(), total=len(POS_LC_df_insilico), desc="POS_LC_In_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_LC.append(SPECTRUM)

    return POS_LC


def re_write_MSP_POS_GC(POS_GC_df):
    POS_GC = []
    for index, row in tqdm(POS_GC_df.iterrows(), total=len(POS_GC_df), desc="\t\t  POS_GC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_GC.append(SPECTRUM)

    return POS_GC

def re_write_MSP_POS_GC_In_Silico(POS_GC_df_insilico):
    POS_GC = []
    for index, row in tqdm(POS_GC_df_insilico.iterrows(), total=len(POS_GC_df_insilico), desc="POS_GC_In_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_GC.append(SPECTRUM)

    return POS_GC

def re_write_MSP_NEG_LC(NEG_LC_df):
    NEG_LC = []
    for index, row in tqdm(NEG_LC_df.iterrows(), total=len(NEG_LC_df), desc="\t\t  NEG_LC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_LC.append(SPECTRUM)

    return NEG_LC

def re_write_MSP_NEG_LC_In_Silico(NEG_LC_df_insilico):
    NEG_LC = []
    for index, row in tqdm(NEG_LC_df_insilico.iterrows(), total=len(NEG_LC_df_insilico), desc="NEG_LC_In_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_LC.append(SPECTRUM)

    return NEG_LC


def re_write_MSP_NEG_GC(NEG_GC_df):
    NEG_GC = []
    for index, row in tqdm(NEG_GC_df.iterrows(), total=len(NEG_GC_df), desc="\t\t  NEG_GC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_GC.append(SPECTRUM)

    return NEG_GC

def re_write_MSP_NEG_GC_In_Silico(NEG_GC_df_insilico):
    NEG_GC = []
    for index, row in tqdm(NEG_GC_df_insilico.iterrows(), total=len(NEG_GC_df_insilico), desc="NEG_GC_in_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_GC.append(SPECTRUM)

    return NEG_GC

def remove_duplicatas(POS_LC,POS_LC_In_Silico,POS_GC,POS_GC_In_Silico,NEG_LC,NEG_LC_In_Silico,NEG_GC,NEG_GC_In_Silico):
    # ========================================================================= POS_LC =========================================================================
    POS_LC_df = remove_dupli_POS_LC(POS_LC)
    # Re convert to MSP
    POS_LC = re_write_MSP_POS_LC(POS_LC_df)

    # ========================================================================= POS_LC_In_Silico =========================================================================
    POS_LC_df_insilico = remove_dupli_POS_LC_In_Silico(POS_LC_In_Silico)
    # Re convert to MSP
    POS_LC_In_Silico = re_write_MSP_POS_LC_In_Silico(POS_LC_df_insilico)

    # ========================================================================= POS_GC =========================================================================
    POS_GC_df = remove_dupli_POS_GC(POS_GC)
    # Re convert to MSP
    POS_GC = re_write_MSP_POS_GC(POS_GC_df)

    # ========================================================================= POS_GC_In_Silico =========================================================================
    POS_GC_df_insilico = remove_dupli_POS_GC_In_Silico(POS_GC_In_Silico)
    # Re convert to MSP
    POS_GC_In_Silico = re_write_MSP_POS_GC_In_Silico(POS_GC_df_insilico)

    # ========================================================================= NEG_LC =========================================================================
    NEG_LC_df = remove_dupli_NEG_LC(NEG_LC)
    # Re convert to MSP
    NEG_LC = re_write_MSP_NEG_LC(NEG_LC_df)

    # ========================================================================= NEG_LC_In_Silico =========================================================================
    NEG_LC_df_insilico = remove_dupli_NEG_LC_In_Silico(NEG_LC_In_Silico)
    # Re convert to MSP
    NEG_LC_In_Silico = re_write_MSP_NEG_LC_In_Silico(NEG_LC_df_insilico)

    # ========================================================================= NEG_GC =========================================================================
    NEG_GC_df = remove_dupli_NEG_GC(NEG_GC)
    # Re convert to MSP
    NEG_GC = re_write_MSP_NEG_GC(NEG_GC_df)

    # ========================================================================= NEG_GC_In_Silico =========================================================================
    NEG_GC_df_insilico = remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico)
    # Re convert to MSP
    NEG_GC_In_Silico = re_write_MSP_NEG_GC_In_Silico(NEG_GC_df_insilico)


    return POS_LC,POS_LC_df,POS_LC_df_insilico,POS_LC_In_Silico,POS_GC,POS_GC_df,POS_GC_df_insilico,POS_GC_In_Silico,NEG_LC,NEG_LC_df,NEG_LC_df_insilico,NEG_LC_In_Silico,NEG_GC,NEG_GC_df,NEG_GC_df_insilico,NEG_GC_In_Silico
