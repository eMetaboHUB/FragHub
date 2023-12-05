from tqdm.auto import tqdm
from tqdm import tqdm
import pandas as pd
import time
import re

def remove_dupli_POS_LC(POS_LC):
    total_rows = len(POS_LC)
    t = tqdm(total=len(POS_LC), desc="\t\t  POS_LC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_LC = POS_LC.loc[~POS_LC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return  POS_LC

def remove_dupli_POS_LC_In_Silico(POS_LC_In_Silico):
    total_rows = len(POS_LC_In_Silico)
    t = tqdm(total=len(POS_LC_In_Silico), desc="POS_LC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_LC_In_Silico = POS_LC_In_Silico.loc[~POS_LC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return POS_LC_In_Silico

def remove_dupli_POS_GC(POS_GC):
    total_rows = len(POS_GC)
    t = tqdm(total=len(POS_GC), desc="\t\t  POS_GC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_GC = POS_GC.loc[~POS_GC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return POS_GC

def remove_dupli_POS_GC_In_Silico(POS_GC_In_Silico):
    total_rows = len(POS_GC_In_Silico)
    t = tqdm(total=len(POS_GC_In_Silico), desc="POS_GC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_GC_In_Silico = POS_GC_In_Silico.loc[~POS_GC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return POS_GC_In_Silico

def remove_dupli_NEG_LC(NEG_LC):
    total_rows = len(NEG_LC)
    t = tqdm(total=len(NEG_LC), desc="\t\t  NEG_LC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_LC = NEG_LC.loc[~NEG_LC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_LC

def remove_dupli_NEG_LC_In_Silico(NEG_LC_In_Silico):
    total_rows = len(NEG_LC_In_Silico)
    t = tqdm(total=len(NEG_LC_In_Silico), desc="NEG_LC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_LC_In_Silico = NEG_LC_In_Silico.loc[~NEG_LC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_LC_In_Silico


def remove_dupli_NEG_GC(NEG_GC):
    total_rows = len(NEG_GC)
    t = tqdm(total=len(NEG_GC), desc="\t\t  NEG_GC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_GC = NEG_GC.loc[~NEG_GC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_GC

def remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico):
    total_rows = len(NEG_GC_In_Silico)
    t = tqdm(total=len(NEG_GC_In_Silico), desc="NEG_GC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_GC_In_Silico = NEG_GC_In_Silico.loc[~NEG_GC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_GC_In_Silico


def re_write_MSP_POS_LC(POS_LC_df):
    POS_LC = []
    for index,row in tqdm(POS_LC_df.iterrows(), total=len(POS_LC_df), desc="\t\t  POS_LC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
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
