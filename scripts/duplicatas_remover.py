from tqdm import tqdm
import time
import re

def find_indices_public(l, value):
    duplicatas_index = [index for index, item in enumerate(l) if item == value]
    if duplicatas_index != None:
        if len(duplicatas_index) >= 2:
            duplicatas_index = duplicatas_index[1:]
            return duplicatas_index
        else:
            return []
    else:
        return []

def remove_duplicatas_public(POS, NEG):
    # POS
    POS_list = []
    POS_index_to_delete = []
    for spectrum in POS: # Generate a given list reduced to inchikey + peak_list
        POS_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n",spectrum).group(1),"PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1) ,"PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)})

    for current_dict in POS_list:
        POS_index_to_delete.extend(find_indices_public(POS_list,current_dict))

    POS_FILTERED = []
    compteur = 0
    print("POS")
    time.sleep(0.01)
    for spectrum in list(tqdm(POS, total=len(POS), unit="spectrums", colour="green")):
        if compteur not in POS_index_to_delete:
            POS_FILTERED.append(spectrum)
        compteur += 1

    # NEG
    NEG_list = []
    NEG_index_to_delete = []
    for spectrum in NEG:  # Generate a given list reduced to inchikey + peak_list
        NEG_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),"PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1), "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2)})

    for current_dict in NEG_list:
        NEG_index_to_delete.extend(find_indices_public(NEG_list, current_dict))

    NEG_FILTERED = []
    compteur = 0
    print("NEG")
    time.sleep(0.01)
    for spectrum in list(tqdm(NEG, total=len(NEG), unit="spectrums", colour="green")):
        if compteur not in NEG_index_to_delete:
            NEG_FILTERED.append(spectrum)
        compteur += 1

    return POS_FILTERED,NEG_FILTERED

def find_indices_lrsv(l, value):
    duplicatas = [(index,item) for index, item in enumerate(l) if item == value]
    duplicatas_index = [index[0] for index in duplicatas if int(index[1]["NPEAKS"]) != max(int(d["NPEAKS"]) for d in duplicatas)]
    if duplicatas_index != None:
        if len(duplicatas_index) >= 2:
            duplicatas_index = duplicatas_index
            return duplicatas_index
        else:
            return []
    else:
        return []
def remove_duplicatas_lrsv(POS, NEG):
    # POS
    POS_list = []
    POS_index_to_delete = []
    for spectrum in POS: # Generate a given list reduced to inchikey + peak_list
        POS_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n",spectrum).group(1),"PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1) ,"NPEAKS": re.search("NUM PEAKS: (.*)\n",spectrum).group(1) ,"PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)})

    for current_dict in POS_list:
        POS_index_to_delete.extend(find_indices_lrsv(POS_list,current_dict))

    POS_FILTERED = []
    compteur = 0
    print("POS")
    time.sleep(0.01)
    for spectrum in list(tqdm(POS, total=len(POS), unit="spectrums", colour="green")):
        if compteur not in POS_index_to_delete:
            POS_FILTERED.append(spectrum)
        compteur += 1

    # NEG
    NEG_list = []
    NEG_index_to_delete = []
    for spectrum in NEG:  # Generate a given list reduced to inchikey + peak_list
        NEG_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),"PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1), "NPEAKS": re.search("NUM PEAKS: (.*)\n",spectrum).group(1), "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2)})

    for current_dict in NEG_list:
        NEG_index_to_delete.extend(find_indices_lrsv(NEG_list, current_dict))

    NEG_FILTERED = []
    compteur = 0
    print("NEG")
    time.sleep(0.01)
    for spectrum in list(tqdm(NEG, total=len(NEG), unit="spectrums", colour="green")):
        if compteur not in NEG_index_to_delete:
            NEG_FILTERED.append(spectrum)
        compteur += 1

    return POS_FILTERED,NEG_FILTERED