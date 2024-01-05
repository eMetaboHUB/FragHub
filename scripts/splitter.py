from tqdm import tqdm
import re

def split_pos_neg(CONCATENATE_DF):
    """
    :param CONCATENATE_DF: DataFrame containing data to be split based on the value of the 'IONMODE' column.
    :return: Tuple containing two DataFrames - 'POS' and 'NEG', representing the subsets of CONCATENATE_DF where the value of 'IONMODE' is 'positive' and 'negative' respectively.
    """
    # Créer une barre de progression pour la première étape
    with tqdm(total=len(CONCATENATE_DF), unit=" spectrums", colour="green", desc="{:>80}".format("POS")) as pbar:
        # Séparer les lignes en fonction de la valeur de la colonne "IONMODE"
        POS = CONCATENATE_DF[CONCATENATE_DF['IONMODE'] == 'positive']

        # Mettre à jour la barre de progression
        pbar.update(len(POS))

    # Créer une barre de progression pour la deuxième étape
    with tqdm(total=len(CONCATENATE_DF), unit=" spectrums", colour="green", desc="{:>80}".format("NEG")) as pbar:
        NEG = CONCATENATE_DF[CONCATENATE_DF['IONMODE'] == 'negative']

        # Mettre à jour la barre de progression
        pbar.update(len(NEG))

    return POS, NEG

def split_LC_GC(POS,NEG):
    """
    :param POS: DataFrame containing positive spectrums
    :param NEG: DataFrame containing negative spectrums
    :return: Four DataFrames containing positive LC spectrums, positive GC spectrums, negative LC spectrums, and negative GC spectrums
    """
    # Séparer les lignes en fonction de la colonne "INSTRUMENTTYPE" pour POS

    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>80}".format("POS_GC")) as pbar:
        POS_GC = POS[POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(POS_GC))


    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>80}".format("POS_LC")) as pbar:
        POS_LC = POS[~POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(POS_LC))

    # Séparer les lignes en fonction de la colonne "INSTRUMENTTYPE" pour NEG
    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>80}".format("NEG_GC")) as pbar:
        NEG_GC = NEG[NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(NEG_GC))

    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>80}".format("NEG_LC")) as pbar:
        NEG_LC = NEG[~NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(NEG_LC))

    return POS_LC, POS_GC, NEG_LC, NEG_GC

def exp_in_silico_splitter(POS_LC,POS_GC,NEG_LC,NEG_GC):
    """
    :param POS_LC: The data frame containing the positive LC data
    :param POS_GC: The data frame containing the positive GC data
    :param NEG_LC: The data frame containing the negative LC data
    :param NEG_GC: The data frame containing the negative GC data
    :return: A tuple containing separate data frames for positive LC in silico, positive GC in silico,
        negative LC in silico, negative GC in silico, positive LC experimental, positive GC experimental,
        negative LC experimental, and negative GC experimental.
    """
    # Barres de progression pour chaque étape de séparation
    with tqdm(total=len(POS_LC), unit=" spectrums", colour="green", desc="{:>80}".format("POS_LC_In_Silico")) as pbar:
        POS_LC_In_Silico_temp = POS_LC[POS_LC['PREDICTED'] == "true"]
        pbar.update(len(POS_LC_In_Silico_temp))

    with tqdm(total=len(POS_GC), unit=" spectrums", colour="green", desc="{:>80}".format("POS_GC_In_Silico")) as pbar:
        POS_GC_In_Silico_temp = POS_GC[POS_GC['PREDICTED'] == "true"]
        pbar.update(len(POS_GC_In_Silico_temp))

    with tqdm(total=len(NEG_LC), unit=" spectrums", colour="green", desc="{:>80}".format("NEG_LC_In_Silico")) as pbar:
        NEG_LC_In_Silico_temp = NEG_LC[NEG_LC['PREDICTED'] == "true"]
        pbar.update(len(NEG_LC_In_Silico_temp))

    with tqdm(total=len(NEG_GC), unit=" spectrums", colour="green", desc="{:>80}".format("NEG_GC_In_Silico")) as pbar:
        NEG_GC_In_Silico_temp = NEG_GC[NEG_GC['PREDICTED'] == "true"]
        pbar.update(len(NEG_GC_In_Silico_temp))

    # Séparation pour les lignes contenant "False"
    with tqdm(total=len(POS_LC), unit=" spectrums", colour="green", desc="{:>80}".format("POS_LC_Exp")) as pbar:
        POS_LC_temp = POS_LC[POS_LC['PREDICTED'] == "false"]
        pbar.update(len(POS_LC_temp))

    with tqdm(total=len(POS_GC), unit=" spectrums", colour="green", desc="{:>80}".format("POS_GC_Exp")) as pbar:
        POS_GC_temp = POS_GC[POS_GC['PREDICTED'] == "false"]
        pbar.update(len(POS_GC_temp))

    with tqdm(total=len(NEG_LC), unit=" spectrums", colour="green", desc="{:>80}".format("NEG_LC_Exp")) as pbar:
        NEG_LC_temp = NEG_LC[NEG_LC['PREDICTED'] == "false"]
        pbar.update(len(NEG_LC_temp))

    with tqdm(total=len(NEG_GC), unit=" spectrums", colour="green", desc="{:>80}".format("NEG_GC_Exp")) as pbar:
        NEG_GC_temp = NEG_GC[NEG_GC['PREDICTED'] == "false"]
        pbar.update(len(NEG_GC_temp))

    return POS_LC_temp,POS_LC_In_Silico_temp,POS_GC_temp,POS_GC_In_Silico_temp,NEG_LC_temp,NEG_LC_In_Silico_temp,NEG_GC_temp,NEG_GC_In_Silico_temp