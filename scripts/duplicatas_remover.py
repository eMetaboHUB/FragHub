from tqdm import tqdm
import pandas as pd
import time
import re

def remove_duplicatas(POS_LC,POS_GC,NEG_LC,NEG_GC):
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
        tqdm.pandas(desc="\t\t POS_LC",colour="green",unit=" row")
        POS_LC_df = POS_LC_df.groupby(['INCHIKEY', 'PEAKS_LIST']).progress_apply(lambda x: x.drop_duplicates())
        # Re convert to MSP
        POS_LC = []
        for index,row in POS_LC_df.iterrows():
            SPECTRUM = ""
            SPECTRUM = SPECTRUM + "FILENAME:" + row["FILENAME"] + "\n"
            SPECTRUM = SPECTRUM + "PREDICTED:" + row["PREDICTED"] + "\n"
            SPECTRUM = SPECTRUM + "SPECTRUMID:" + row["SPECTRUMID"] + "\n"
            SPECTRUM = SPECTRUM + "RESOLUTION:" + row["RESOLUTION"] + "\n"
            SPECTRUM = SPECTRUM + "SYNON:" + row["SYNON"] + "\n"
            SPECTRUM = SPECTRUM + "CHARGE:" + row["CHARGE"] + "\n"
            SPECTRUM = SPECTRUM + "PARENTMASS:" + row["PARENTMASS"] + "\n"
            SPECTRUM = SPECTRUM + "IONIZATION:" + row["IONIZATION"] + "\n"
            SPECTRUM = SPECTRUM + "MSLEVEL:" + row["MSLEVEL"] + "\n"
            SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE:" + row["FRAGMENTATIONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "NAME:" + row["NAME"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORMZ:" + row["PRECURSORMZ"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORTYPE:" + row["PRECURSORTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENTTYPE:" + row["INSTRUMENTTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENT:" + row["INSTRUMENT"] + "\n"
            SPECTRUM = SPECTRUM + "SMILES:" + row["SMILES"] + "\n"
            SPECTRUM = SPECTRUM + "INCHI:" + row["INCHI"] + "\n"
            SPECTRUM = SPECTRUM + "INCHIKEY:" + row["INCHIKEY"] + "\n"
            SPECTRUM = SPECTRUM + "COLLISIONENERGY:" + row["COLLISIONENERGY"] + "\n"
            SPECTRUM = SPECTRUM + "FORMULA:" + row["FORMULA"] + "\n"
            SPECTRUM = SPECTRUM + "RETENTIONTIME:" + row["RETENTIONTIME"] + "\n"
            SPECTRUM = SPECTRUM + "IONMODE:" + row["IONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "COMMENT:" + row["COMMENT"] + "\n"
            SPECTRUM = SPECTRUM + "NUM PEAKS:" + row["NUM PEAKS"] + "\n"
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
        tqdm.pandas(desc="\t\t POS_GC", colour="green",unit=" row")
        POS_GC_df = POS_GC_df.groupby(['INCHIKEY', 'PEAKS_LIST']).progress_apply(lambda x: x.drop_duplicates())
        # Re convert to MSP
        POS_GC = []
        for index, row in POS_GC_df.iterrows():
            SPECTRUM = ""
            SPECTRUM = SPECTRUM + "FILENAME:" + row["FILENAME"] + "\n"
            SPECTRUM = SPECTRUM + "PREDICTED:" + row["PREDICTED"] + "\n"
            SPECTRUM = SPECTRUM + "SPECTRUMID:" + row["SPECTRUMID"] + "\n"
            SPECTRUM = SPECTRUM + "RESOLUTION:" + row["RESOLUTION"] + "\n"
            SPECTRUM = SPECTRUM + "SYNON:" + row["SYNON"] + "\n"
            SPECTRUM = SPECTRUM + "CHARGE:" + row["CHARGE"] + "\n"
            SPECTRUM = SPECTRUM + "PARENTMASS:" + row["PARENTMASS"] + "\n"
            SPECTRUM = SPECTRUM + "IONIZATION:" + row["IONIZATION"] + "\n"
            SPECTRUM = SPECTRUM + "MSLEVEL:" + row["MSLEVEL"] + "\n"
            SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE:" + row["FRAGMENTATIONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "NAME:" + row["NAME"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORMZ:" + row["PRECURSORMZ"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORTYPE:" + row["PRECURSORTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENTTYPE:" + row["INSTRUMENTTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENT:" + row["INSTRUMENT"] + "\n"
            SPECTRUM = SPECTRUM + "SMILES:" + row["SMILES"] + "\n"
            SPECTRUM = SPECTRUM + "INCHI:" + row["INCHI"] + "\n"
            SPECTRUM = SPECTRUM + "INCHIKEY:" + row["INCHIKEY"] + "\n"
            SPECTRUM = SPECTRUM + "COLLISIONENERGY:" + row["COLLISIONENERGY"] + "\n"
            SPECTRUM = SPECTRUM + "FORMULA:" + row["FORMULA"] + "\n"
            SPECTRUM = SPECTRUM + "RETENTIONTIME:" + row["RETENTIONTIME"] + "\n"
            SPECTRUM = SPECTRUM + "IONMODE:" + row["IONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "COMMENT:" + row["COMMENT"] + "\n"
            SPECTRUM = SPECTRUM + "NUM PEAKS:" + row["NUM PEAKS"] + "\n"
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
        tqdm.pandas(desc="\t\t NEG_LC", colour="green",unit=" row")
        NEG_LC_df = NEG_LC_df.groupby(['INCHIKEY', 'PEAKS_LIST']).progress_apply(lambda x: x.drop_duplicates())
        # Re convert to MSP
        NEG_LC = []
        for index, row in NEG_LC_df.iterrows():
            SPECTRUM = ""
            SPECTRUM = SPECTRUM + "FILENAME:" + row["FILENAME"] + "\n"
            SPECTRUM = SPECTRUM + "PREDICTED:" + row["PREDICTED"] + "\n"
            SPECTRUM = SPECTRUM + "SPECTRUMID:" + row["SPECTRUMID"] + "\n"
            SPECTRUM = SPECTRUM + "RESOLUTION:" + row["RESOLUTION"] + "\n"
            SPECTRUM = SPECTRUM + "SYNON:" + row["SYNON"] + "\n"
            SPECTRUM = SPECTRUM + "CHARGE:" + row["CHARGE"] + "\n"
            SPECTRUM = SPECTRUM + "PARENTMASS:" + row["PARENTMASS"] + "\n"
            SPECTRUM = SPECTRUM + "IONIZATION:" + row["IONIZATION"] + "\n"
            SPECTRUM = SPECTRUM + "MSLEVEL:" + row["MSLEVEL"] + "\n"
            SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE:" + row["FRAGMENTATIONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "NAME:" + row["NAME"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORMZ:" + row["PRECURSORMZ"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORTYPE:" + row["PRECURSORTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENTTYPE:" + row["INSTRUMENTTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENT:" + row["INSTRUMENT"] + "\n"
            SPECTRUM = SPECTRUM + "SMILES:" + row["SMILES"] + "\n"
            SPECTRUM = SPECTRUM + "INCHI:" + row["INCHI"] + "\n"
            SPECTRUM = SPECTRUM + "INCHIKEY:" + row["INCHIKEY"] + "\n"
            SPECTRUM = SPECTRUM + "COLLISIONENERGY:" + row["COLLISIONENERGY"] + "\n"
            SPECTRUM = SPECTRUM + "FORMULA:" + row["FORMULA"] + "\n"
            SPECTRUM = SPECTRUM + "RETENTIONTIME:" + row["RETENTIONTIME"] + "\n"
            SPECTRUM = SPECTRUM + "IONMODE:" + row["IONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "COMMENT:" + row["COMMENT"] + "\n"
            SPECTRUM = SPECTRUM + "NUM PEAKS:" + row["NUM PEAKS"] + "\n"
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
        tqdm.pandas(desc="\t\t NEG_GC", colour="green",unit=" row")
        NEG_GC_df = NEG_GC_df.groupby(['INCHIKEY', 'PEAKS_LIST']).progress_apply(lambda x: x.drop_duplicates())
        # Re convert to MSP
        NEG_GC = []
        for index, row in NEG_GC_df.iterrows():
            SPECTRUM = ""
            SPECTRUM = SPECTRUM + "FILENAME:" + row["FILENAME"] + "\n"
            SPECTRUM = SPECTRUM + "PREDICTED:" + row["PREDICTED"] + "\n"
            SPECTRUM = SPECTRUM + "SPECTRUMID:" + row["SPECTRUMID"] + "\n"
            SPECTRUM = SPECTRUM + "RESOLUTION:" + row["RESOLUTION"] + "\n"
            SPECTRUM = SPECTRUM + "SYNON:" + row["SYNON"] + "\n"
            SPECTRUM = SPECTRUM + "CHARGE:" + row["CHARGE"] + "\n"
            SPECTRUM = SPECTRUM + "PARENTMASS:" + row["PARENTMASS"] + "\n"
            SPECTRUM = SPECTRUM + "IONIZATION:" + row["IONIZATION"] + "\n"
            SPECTRUM = SPECTRUM + "MSLEVEL:" + row["MSLEVEL"] + "\n"
            SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE:" + row["FRAGMENTATIONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "NAME:" + row["NAME"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORMZ:" + row["PRECURSORMZ"] + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORTYPE:" + row["PRECURSORTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENTTYPE:" + row["INSTRUMENTTYPE"] + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENT:" + row["INSTRUMENT"] + "\n"
            SPECTRUM = SPECTRUM + "SMILES:" + row["SMILES"] + "\n"
            SPECTRUM = SPECTRUM + "INCHI:" + row["INCHI"] + "\n"
            SPECTRUM = SPECTRUM + "INCHIKEY:" + row["INCHIKEY"] + "\n"
            SPECTRUM = SPECTRUM + "COLLISIONENERGY:" + row["COLLISIONENERGY"] + "\n"
            SPECTRUM = SPECTRUM + "FORMULA:" + row["FORMULA"] + "\n"
            SPECTRUM = SPECTRUM + "RETENTIONTIME:" + row["RETENTIONTIME"] + "\n"
            SPECTRUM = SPECTRUM + "IONMODE:" + row["IONMODE"] + "\n"
            SPECTRUM = SPECTRUM + "COMMENT:" + row["COMMENT"] + "\n"
            SPECTRUM = SPECTRUM + "NUM PEAKS:" + row["NUM PEAKS"] + "\n"
            SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
            NEG_GC.append(SPECTRUM)

    return POS_LC,POS_LC_df,POS_GC,POS_GC_df,NEG_LC,NEG_LC_df,NEG_GC,NEG_GC_df
