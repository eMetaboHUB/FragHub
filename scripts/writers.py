import pandas as pd
import os
import re

def writting_msp(clean_msp_path,POS_LC,POS_GC,NEG_LC,NEG_GC,POS_LC_In_Silico,POS_GC_In_Silico,NEG_LC_In_Silico,NEG_GC_In_Silico):
    """
    :param clean_msp_path: The path where the cleaned MSP files will be saved.
    :param POS_LC: A list of strings representing the positive LC data.
    :param POS_GC: A list of strings representing the positive GC data.
    :param NEG_LC: A list of strings representing the negative LC data.
    :param NEG_GC: A list of strings representing the negative GC data.
    :param POS_LC_In_Silico: A list of strings representing the positive LC In Silico data.
    :param POS_GC_In_Silico: A list of strings representing the positive GC In Silico data.
    :param NEG_LC_In_Silico: A list of strings representing the negative LC In Silico data.
    :param NEG_GC_In_Silico: A list of strings representing the negative GC In Silico data.
    :return: None

    This method writes the cleaned MSP files to the specified `clean_msp_path` location. The cleaned data is obtained by joining the strings in the respective input lists and applying regular
    * expression substitutions.

    The cleaned MSP files are saved with the following names and paths:
    - 'POS_LC_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'POS_GC_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'NEG_LC_clean.msp' in the 'NEG' folder under `clean_msp_path`
    - 'NEG_GC_clean.msp' in the 'NEG' folder under `clean_msp_path`
    - 'POS_LC_In_Silico_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'POS_GC_In_Silico_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'NEG_LC_In_Silico_clean.msp' in the 'NEG' folder under `clean_msp_path`
    - 'NEG_GC_In_Silico_clean.msp' in the 'NEG' folder under `clean_msp_path`
    """
    POS_LC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(POS_LC))
    del POS_LC
    POS_GC_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_GC))
    del POS_GC
    NEG_LC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG_LC))
    del NEG_LC
    NEG_GC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG_GC))
    del NEG_GC

    POS_LC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_LC_In_Silico))
    del POS_LC_In_Silico
    POS_GC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_GC_In_Silico))
    del POS_GC_In_Silico
    NEG_LC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(NEG_LC_In_Silico))
    del NEG_LC_In_Silico
    NEG_GC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(NEG_GC_In_Silico))
    del NEG_GC_In_Silico

    with open(os.path.join(clean_msp_path,"POS/POS_LC_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_LC_FULL)
    del POS_LC_FULL
    with open(os.path.join(clean_msp_path,"POS/POS_GC_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_GC_FULL)
    del POS_GC_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_LC_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_LC_FULL)
    del NEG_LC_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_GC_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_GC_FULL)
    del NEG_GC_FULL

    with open(os.path.join(clean_msp_path,"POS/POS_LC_In_Silico_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_LC_In_Silico_FULL)
    del POS_LC_In_Silico_FULL
    with open(os.path.join(clean_msp_path,"POS/POS_GC_In_Silico_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_GC_In_Silico_FULL)
    del POS_GC_In_Silico_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_LC_In_Silico_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_LC_In_Silico_FULL)
    del NEG_LC_In_Silico_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_GC_In_Silico_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_GC_In_Silico_FULL)
    del NEG_GC_In_Silico_FULL

def writting_csv(POS_LC_df,POS_GC_df,NEG_LC_df,NEG_GC_df,POS_LC_df_insilico,POS_GC_df_insilico,NEG_LC_df_insilico,NEG_GC_df_insilico):
    """
    Method to write dataframes to CSV files.

    :param POS_LC_df: Dataframe containing positive LC data.
    :param POS_GC_df: Dataframe containing positive GC data.
    :param NEG_LC_df: Dataframe containing negative LC data.
    :param NEG_GC_df: Dataframe containing negative GC data.
    :param POS_LC_df_insilico: Dataframe containing positive LC in silico data.
    :param POS_GC_df_insilico: Dataframe containing positive GC in silico data.
    :param NEG_LC_df_insilico: Dataframe containing negative LC in silico data.
    :param NEG_GC_df_insilico: Dataframe containing negative GC in silico data.
    :return: None
    """
    POS_LC_df.to_csv("../OUTPUT/CSV/POS/POS_LC_clean.csv", sep=";", encoding="UTF-8", index=False)
    del POS_LC_df
    POS_GC_df.to_csv("../OUTPUT/CSV/POS/POS_GC_clean.csv", sep=";", encoding="UTF-8", index=False)
    del POS_GC_df
    NEG_LC_df.to_csv("../OUTPUT/CSV/NEG/NEG_LC_clean.csv", sep=";", encoding="UTF-8", index=False)
    del NEG_LC_df
    NEG_GC_df.to_csv("../OUTPUT/CSV/NEG/NEG_GC_clean.csv", sep=";", encoding="UTF-8", index=False)
    del NEG_GC_df

    POS_LC_df_insilico.to_csv("../OUTPUT/CSV/POS/POS_LC_In_Silico_clean.csv", sep=";", encoding="UTF-8", index=False)
    del POS_LC_df_insilico
    POS_GC_df_insilico.to_csv("../OUTPUT/CSV/POS/POS_GC_In_Silico_clean.csv", sep=";", encoding="UTF-8", index=False)
    del POS_GC_df_insilico
    NEG_LC_df_insilico.to_csv("../OUTPUT/CSV/NEG/NEG_LC_In_Silico_clean.csv", sep=";", encoding="UTF-8", index=False)
    del NEG_LC_df_insilico
    NEG_GC_df_insilico.to_csv("../OUTPUT/CSV/NEG/NEG_GC_In_Silico_clean.csv", sep=";", encoding="UTF-8", index=False)
    del NEG_GC_df_insilico