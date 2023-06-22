import re
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