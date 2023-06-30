from tqdm import tqdm
import re
def split_pos_neg(CONCATENATE_LIST):
    POS = []
    NEG = []
    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit="spectrums", colour="green", desc="\tprocessing"):
        #POS
        if re.search("CHARGE: [0-9]\n",spectrum, flags=re.I) or re.search("PRECURSORTYPE: (.*)\+\n",spectrum, flags=re.I) or re.search("IONMODE: p(.*)\n",spectrum, flags=re.I):
            POS.append(spectrum)
        # NEG
        elif re.search("CHARGE: \-[0-9]\n",spectrum, flags=re.I) or re.search("PRECURSORTYPE: (.*)\-\n",spectrum, flags=re.I) or re.search("IONMODE: n(.*)\n",spectrum, flags=re.I):
            NEG.append(spectrum)

    return POS, NEG

def split_LC_GC(POS,NEG):
    # POS
    POS_LC = []
    POS_GC = []
    for spectrums in tqdm(POS, total=len(POS), unit="spectrums", colour="green", desc="\t\t   POS"):
        if re.search("(?<![a-zA-Z0-9])GC(?![a-zA-Z0-9])",spectrums,flags=re.I) or re.search("(?<![a-zA-Z0-9])EI(?![a-zA-Z0-9])",spectrums,flags=re.I):
            POS_GC.append(spectrums)
        else:
            POS_LC.append(spectrums)

    # NEG
    NEG_LC = []
    NEG_GC = []
    for spectrums in tqdm(NEG, total=len(NEG), unit="spectrums", colour="green", desc="\t\t   NEG"):
        if re.search("(?<![a-zA-Z0-9])GC(?![a-zA-Z0-9])", spectrums, flags=re.I) or re.search("(?<![a-zA-Z0-9])EI(?![a-zA-Z0-9])", spectrums, flags=re.I):
            NEG_GC.append(spectrums)
        else:
            NEG_LC.append(spectrums)

    return POS_LC,POS_GC,NEG_LC,NEG_GC
