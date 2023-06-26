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

def remove_duplicatas_public(POS_LC,POS_GC,NEG_LC,NEG_GC):
    # ========================================================================= POS_LC =========================================================================
    POS_LC_list = []
    POS_index_to_delete = []
    for spectrum in POS_LC:  # Generate a given list reduced to inchikey + peak_list
        POS_LC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                            "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)})

    for current_dict in POS_LC_list:
        POS_index_to_delete.extend(find_indices_public(POS_LC_list, current_dict))

    POS_LC_FILTERED = []
    compteur = 0
    print("POS_LC")
    time.sleep(0.01)
    for spectrum in list(tqdm(POS_LC, total=len(POS_LC), unit="spectrums", colour="green")):
        if compteur not in POS_index_to_delete:
            POS_LC_FILTERED.append(spectrum)
        compteur += 1
    #  ========================================================================= POS_GC =========================================================================
    POS_GC_list = []
    POS_index_to_delete = []
    for spectrum in POS_GC:  # Generate a given list reduced to inchikey + peak_list
        POS_GC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                            "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)})

    for current_dict in POS_GC_list:
        POS_index_to_delete.extend(find_indices_public(POS_GC_list, current_dict))

    POS_GC_FILTERED = []
    compteur = 0
    print("POS_GC")
    time.sleep(0.01)
    for spectrum in list(tqdm(POS_GC, total=len(POS_GC), unit="spectrums", colour="green")):
        if compteur not in POS_index_to_delete:
            POS_GC_FILTERED.append(spectrum)
        compteur += 1

    #  ========================================================================= NEG_LC =========================================================================
    NEG_LC_list = []
    NEG_index_to_delete = []
    for spectrum in NEG_LC:  # Generate a given list reduced to inchikey + peak_list
        NEG_LC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                            "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)})

    for current_dict in NEG_LC_list:
        NEG_index_to_delete.extend(find_indices_public(NEG_LC_list, current_dict))

    NEG_LC_FILTERED = []
    compteur = 0
    print("NEG_LC")
    time.sleep(0.01)
    for spectrum in list(tqdm(NEG_LC, total=len(NEG_LC), unit="spectrums", colour="green")):
        if compteur not in NEG_index_to_delete:
            NEG_LC_FILTERED.append(spectrum)
        compteur += 1

    #  ========================================================================= NEG_GC =========================================================================
    NEG_GC_list = []
    NEG_index_to_delete = []
    for spectrum in NEG_GC:  # Generate a given list reduced to inchikey + peak_list
        NEG_GC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n", spectrum).group(1),
                            "PEAK_LIST": re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)})

    for current_dict in NEG_GC_list:
        NEG_index_to_delete.extend(find_indices_public(NEG_GC_list, current_dict))

    NEG_GC_FILTERED = []
    compteur = 0
    print("NEG_GC")
    time.sleep(0.01)
    for spectrum in list(tqdm(NEG_GC, total=len(NEG_GC), unit="spectrums", colour="green")):
        if compteur not in NEG_index_to_delete:
            NEG_GC_FILTERED.append(spectrum)
        compteur += 1

    return POS_LC_FILTERED, POS_GC_FILTERED, NEG_LC_FILTERED, NEG_GC_FILTERED

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
        POS_LC_list.append({"INCHIKEY": re.search("INCHIKEY: (.*)\n",spectrum).group(1),
                         "PRECURSORTYPE": re.search("PRECURSORTYPE: (.*)\n",spectrum).group(1) ,
                         "NPEAKS": re.search("NUM PEAKS: (.*)\n",spectrum).group(1)})

    for current_dict in POS_LC_list:
        POS_index_to_delete.extend(find_indices_lrsv(POS_LC_list,current_dict))

    POS_LC_FILTERED = []
    compteur = 0
    print("POS_LC")
    time.sleep(0.01)
    for spectrum in list(tqdm(POS_LC, total=len(POS_LC), unit="spectrums", colour="green")):
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
    print("POS_GC")
    time.sleep(0.01)
    for spectrum in list(tqdm(POS_GC, total=len(POS_GC), unit="spectrums", colour="green")):
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
    print("NEG_LC")
    time.sleep(0.01)
    for spectrum in list(tqdm(NEG_LC, total=len(NEG_LC), unit="spectrums", colour="green")):
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
    print("NEG_GC")
    time.sleep(0.01)
    for spectrum in list(tqdm(NEG_GC, total=len(NEG_GC), unit="spectrums", colour="green")):
        if compteur not in NEG_index_to_delete:
            NEG_GC_FILTERED.append(spectrum)
        compteur += 1

    return POS_LC_FILTERED,POS_GC_FILTERED,NEG_LC_FILTERED,NEG_GC_FILTERED