import re
def split_pos_neg(CONCATENATE_LIST):
    POS = []
    NEG = []
    for spectrum in CONCATENATE_LIST:
        #POS
        if re.search("CHARGE: [0-9]\n",spectrum, flags=re.I) or re.search("PRECURSORTYPE: (.*)\+\n",spectrum, flags=re.I) or re.search("IONMODE: p(.*)\n",spectrum, flags=re.I):
            POS.append(spectrum)
        # NEG
        elif re.search("CHARGE: \-[0-9]\n",spectrum, flags=re.I) or re.search("PRECURSORTYPE: (.*)\-\n",spectrum, flags=re.I) or re.search("IONMODE: n(.*)\n",spectrum, flags=re.I):
            NEG.append(spectrum)

    return POS, NEG

# def split_LC_GC(POS,NEG):
