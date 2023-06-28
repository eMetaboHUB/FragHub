from tqdm import tqdm
import pandas as pd
import time
import re

# def find_indices_public(l, value):
#     duplicatas_index = [index for index, item in enumerate(l) if item == value]
#     if duplicatas_index != None:
#         if len(duplicatas_index) >= 2:
#             duplicatas_index = duplicatas_index[1:]
#             return duplicatas_index
#         else:
#             return []
#     else:
#         return []
#
# def remove_duplicatas_public(POS_LC,POS_GC,NEG_LC,NEG_GC):
#     # ========================================================================= POS_LC =========================================================================
#     POS_LC_list = []
#     POS_index_to_delete = []
#     for spectrum in POS_LC:  # Generate a given list reduced to inchikey + peak_list
#         POS_LC_list.append(str(re.search("INCHIKEY: (.*)\n", spectrum).group(1) + re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)))
#
#     for current_dict in POS_LC_list:
#         POS_index_to_delete.extend(find_indices_public(POS_LC_list, current_dict))
#
#     POS_LC_FILTERED = []
#     compteur = 0
#
#     for spectrum in list(tqdm(POS_LC, total=len(POS_LC), unit="spectrums", colour="green", desc="\tPOS_LC")):
#         if compteur not in POS_index_to_delete:
#             POS_LC_FILTERED.append(spectrum)
#         compteur += 1
#     #  ========================================================================= POS_GC =========================================================================
#     POS_GC_list = []
#     POS_index_to_delete = []
#     for spectrum in POS_GC:  # Generate a given list reduced to inchikey + peak_list
#         POS_GC_list.append(str(re.search("INCHIKEY: (.*)\n", spectrum).group(1) + re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)))
#
#     for current_dict in POS_GC_list:
#         POS_index_to_delete.extend(find_indices_public(POS_GC_list, current_dict))
#
#     POS_GC_FILTERED = []
#     compteur = 0
#
#     for spectrum in list(tqdm(POS_GC, total=len(POS_GC), unit="spectrums", colour="green", desc="\tPOS_GC")):
#         if compteur not in POS_index_to_delete:
#             POS_GC_FILTERED.append(spectrum)
#         compteur += 1
#
#     #  ========================================================================= NEG_LC =========================================================================
#     NEG_LC_list = []
#     NEG_index_to_delete = []
#     for spectrum in NEG_LC:  # Generate a given list reduced to inchikey + peak_list
#         NEG_LC_list.append(str(re.search("INCHIKEY: (.*)\n", spectrum).group(1) + re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)))
#
#     for current_dict in NEG_LC_list:
#         NEG_index_to_delete.extend(find_indices_public(NEG_LC_list, current_dict))
#
#     NEG_LC_FILTERED = []
#     compteur = 0
#
#     for spectrum in list(tqdm(NEG_LC, total=len(NEG_LC), unit="spectrums", colour="green", desc="\tNEG_LC")):
#         if compteur not in NEG_index_to_delete:
#             NEG_LC_FILTERED.append(spectrum)
#         compteur += 1
#
#     #  ========================================================================= NEG_GC =========================================================================
#     NEG_GC_list = []
#     NEG_index_to_delete = []
#     for spectrum in NEG_GC:  # Generate a given list reduced to inchikey + peak_list
#         NEG_GC_list.append(str(re.search("INCHIKEY: (.*)\n", spectrum).group(1) + re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)))
#
#     for current_dict in NEG_GC_list:
#         NEG_index_to_delete.extend(find_indices_public(NEG_GC_list, current_dict))
#
#     NEG_GC_FILTERED = []
#     compteur = 0
#
#     for spectrum in list(tqdm(NEG_GC, total=len(NEG_GC), unit="spectrums", colour="green", desc="\tNEG_GC")):
#         if compteur not in NEG_index_to_delete:
#             NEG_GC_FILTERED.append(spectrum)
#         compteur += 1
#
#     return POS_LC_FILTERED, POS_GC_FILTERED, NEG_LC_FILTERED, NEG_GC_FILTERED

def remove_duplicatas_public(POS_LC,POS_GC,NEG_LC,NEG_GC):
    # ========================================================================= POS_LC =========================================================================
    dictionary = {}
    POS_LC_df = pd.DataFrame()
    first = True
    empty = False
    if len(POS_LC) == 0:
        empty = True

    for spectrum in POS_LC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first == True:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if first == False:
                for element in fields:
                    dictionary[element[0]].append(element[1])

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if empty == False:
        # Creating Dataframe
        POS_LC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        POS_LC_df = POS_LC_df.drop_duplicates(subset=['INCHIKEY', 'PEAKS_LIST'])
        # Re convert to MSP
        POS_LC = []
        for index,row in POS_LC_df.iterrows():
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
    # ========================================================================= POS_GC =========================================================================
    dictionary = {}
    POS_GC_df = pd.DataFrame()
    first = True
    empty = False
    if len(POS_GC) == 0:
        empty = True

    for spectrum in POS_GC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first == True:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if first == False:
                for element in fields:
                    dictionary[element[0]].append(element[1])

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if empty == False:
        # Creating Dataframe
        POS_GC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        POS_GC_df = POS_GC_df.drop_duplicates(subset=['INCHIKEY', 'PEAKS_LIST'])
        # Re convert to MSP
        POS_GC = []
        for index, row in POS_GC_df.iterrows():
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
    # ========================================================================= NEG_LC =========================================================================
    dictionary = {}
    NEG_LC_df = pd.DataFrame()
    first = True
    empty = False
    if len(NEG_LC) == 0:
        empty = True

    for spectrum in NEG_LC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first == True:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if first == False:
                for element in fields:
                    dictionary[element[0]].append(element[1])

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if empty == False:
        # Creating Dataframe
        NEG_LC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        NEG_LC_df = NEG_LC_df.drop_duplicates(subset=['INCHIKEY', 'PEAKS_LIST'])
        # Re convert to MSP
        NEG_LC = []
        for index, row in NEG_LC_df.iterrows():
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
    # ========================================================================= NEG_GC =========================================================================
    dictionary = {}
    NEG_GC_df = pd.DataFrame()
    first = True
    empty = False
    if len(NEG_GC) == 0:
        empty = True

    for spectrum in NEG_GC:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first == True:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if first == False:
                for element in fields:
                    dictionary[element[0]].append(element[1])

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if empty == False:
        # Creating Dataframe
        NEG_GC_df = pd.DataFrame.from_dict(dictionary)
        # Removing duplicatas
        NEG_GC_df = NEG_GC_df.drop_duplicates(subset=['INCHIKEY', 'PEAKS_LIST'])
        # Re convert to MSP
        NEG_GC = []
        for index, row in NEG_GC_df.iterrows():
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

    return POS_LC,POS_LC_df,POS_GC,POS_GC_df,NEG_LC,NEG_LC_df,NEG_GC,NEG_GC_df

def find_indices_lrsv(l, value):
    duplicatas = [[index,item] for index, item in enumerate(l) if item == value]

    duplicatas_index = [index[0] for index in duplicatas if int(index[1]["NPEAKS"]) != max(int(d[1]["NPEAKS"]) for d in duplicatas)]
    if duplicatas_index != None:
        if len(duplicatas_index) >= 2:
            return duplicatas_index
        else:
            return []
    else:
        return []
def remove_duplicatas_lrsv(POS_LC,POS_GC,NEG_LC,NEG_GC):
    # ========================================================================= POS_LC =========================================================================
    POS_LC_list = []
    POS_index_to_delete = []
    for spectrum in POS_LC: # Generate a given list reduced to inchikey + peak_list
        POS_LC_list.append(str({"INCHIKEY": re.search("INCHIKEY: (.*)\n",spectrum).group(1),
                         "PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1) ,
                         "NPEAKS": re.search("NUM PEAKS: (.*)\n",spectrum).group(1)}))

    for current_dict in POS_LC_list:
        POS_index_to_delete.extend(find_indices_lrsv(POS_LC_list,current_dict))

    POS_LC_FILTERED = []
    compteur = 0

    for spectrum in list(tqdm(POS_LC, total=len(POS_LC), unit="spectrums", colour="green", desc="\tPOS_LC")):
        if compteur not in POS_index_to_delete:
            POS_LC_FILTERED.append(spectrum)
        compteur += 1
    #  ========================================================================= POS_GC =========================================================================
    POS_GC_list = []
    POS_index_to_delete = []
    for spectrum in POS_GC:  # Generate a given list reduced to inchikey + peak_list
        POS_GC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                            "PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n", spectrum).group(1),
                            "NPEAKS": re.search("NUM PEAKS: (.*)\n", spectrum).group(1)})

    for current_dict in POS_GC_list:
        POS_index_to_delete.extend(find_indices_lrsv(POS_GC_list, current_dict))

    POS_GC_FILTERED = []
    compteur = 0

    for spectrum in list(tqdm(POS_GC, total=len(POS_GC), unit="spectrums", colour="green", desc="\tPOS_GC")):
        if compteur not in POS_index_to_delete:
            POS_GC_FILTERED.append(spectrum)
        compteur += 1

    #  ========================================================================= NEG_LC =========================================================================
    NEG_LC_list = []
    NEG_index_to_delete = []
    for spectrum in NEG_LC:  # Generate a given list reduced to inchikey + peak_list
        NEG_LC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                         "PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1),
                         "NPEAKS": re.search("NUM PEAKS: (.*)\n",spectrum).group(1)})

    for current_dict in NEG_LC_list:
        NEG_index_to_delete.extend(find_indices_lrsv(NEG_LC_list, current_dict))

    NEG_LC_FILTERED = []
    compteur = 0

    for spectrum in list(tqdm(NEG_LC, total=len(NEG_LC), unit="spectrums", colour="green", desc="\tNEG_LC")):
        if compteur not in NEG_index_to_delete:
            NEG_LC_FILTERED.append(spectrum)
        compteur += 1

    #  ========================================================================= NEG_GC =========================================================================
    NEG_GC_list = []
    NEG_index_to_delete = []
    for spectrum in NEG_GC:  # Generate a given list reduced to inchikey + peak_list
        NEG_GC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                         "PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n", spectrum).group(1),
                         "NPEAKS": re.search("NUM PEAKS: (.*)\n", spectrum).group(1)})

    for current_dict in NEG_GC_list:
        NEG_index_to_delete.extend(find_indices_lrsv(NEG_GC_list, current_dict))

    NEG_GC_FILTERED = []
    compteur = 0

    for spectrum in list(tqdm(NEG_GC, total=len(NEG_GC), unit="spectrums", colour="green", desc="\tNEG_GC")):
        if compteur not in NEG_index_to_delete:
            NEG_GC_FILTERED.append(spectrum)
        compteur += 1

    return POS_LC_FILTERED,POS_GC_FILTERED,NEG_LC_FILTERED,NEG_GC_FILTERED