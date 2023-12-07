from tqdm import tqdm
import re

def split_pos_neg(CONCATENATE_DF):
    """
    Split the given DataFrame, CONCATENATE_DF, into two separate DataFrames based on the value of the 'IONMODE' column.

    :param CONCATENATE_DF: The DataFrame to be split.
    :return: Two DataFrames, POS and NEG, containing the rows with 'IONMODE' values of 'positive' and 'negative' respectively.

    """
    # Créer une barre de progression pour la première étape
    with tqdm(total=len(CONCATENATE_DF), unit=" spectrums", colour="green", desc="\t  POS") as pbar:
        # Séparer les lignes en fonction de la valeur de la colonne "IONMODE"
        POS = CONCATENATE_DF[CONCATENATE_DF['IONMODE'] == 'positive']

        # Mettre à jour la barre de progression
        pbar.update(len(POS))

    # Créer une barre de progression pour la deuxième étape
    with tqdm(total=len(CONCATENATE_DF), unit=" spectrums", colour="green", desc="\t  NEG") as pbar:
        NEG = CONCATENATE_DF[CONCATENATE_DF['IONMODE'] == 'negative']

        # Mettre à jour la barre de progression
        pbar.update(len(NEG))

    return POS, NEG

def split_LC_GC(POS,NEG):
    """
    Split the input data into LC and GC based on the "INSTRUMENTTYPE" column.

    :param POS: Positive data.
    :param NEG: Negative data.
    :return: Tuple containing POS_LC, POS_GC, NEG_LC, and NEG_GC.
    """
    # Séparer les lignes en fonction de la colonne "INSTRUMENTTYPE" pour POS

    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="\t  POS_GC") as pbar:
        POS_GC = POS[POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(POS_GC))


    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="\t  POS_LC") as pbar:
        POS_LC = POS[~POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(POS_LC))

    # Séparer les lignes en fonction de la colonne "INSTRUMENTTYPE" pour NEG
    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="\t  NEG_GC") as pbar:
        NEG_GC = NEG[NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(NEG_GC))

    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="\t  NEG_LC") as pbar:
        NEG_LC = NEG[~NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]

        # Mettre à jour la barre de progression
        pbar.update(len(NEG_LC))

    return POS_LC, POS_GC, NEG_LC, NEG_GC

def exp_in_silico_splitter(POS_LC,POS_GC,NEG_LC,NEG_GC):
    """
    :param POS_LC: The input DataFrame for positive LC spectrums.
    :param POS_GC: The input DataFrame for positive GC spectrums.
    :param NEG_LC: The input DataFrame for negative LC spectrums.
    :param NEG_GC: The input DataFrame for negative GC spectrums.
    :return: A tuple containing the following DataFrames:
        - POS_LC_temp: The separated DataFrame for positive LC spectrums with 'PREDICTED' column value "false".
        - POS_LC_In_Silico_temp: The separated DataFrame for positive LC spectrums with 'PREDICTED' column value "true".
        - POS_GC_temp: The separated DataFrame for positive GC spectrums with 'PREDICTED' column value "false".
        - POS_GC_In_Silico_temp: The separated DataFrame for positive GC spectrums with 'PREDICTED' column value "true".
        - NEG_LC_temp: The separated DataFrame for negative LC spectrums with 'PREDICTED' column value "false".
        - NEG_LC_In_Silico_temp: The separated DataFrame for negative LC spectrums with 'PREDICTED' column value "true".
        - NEG_GC_temp: The separated DataFrame for negative GC spectrums with 'PREDICTED' column value "false".
        - NEG_GC_In_Silico_temp: The separated DataFrame for negative GC spectrums with 'PREDICTED' column value "true".
    """
    # Barres de progression pour chaque étape de séparation
    with tqdm(total=len(POS_LC), unit=" spectrums", colour="green", desc="\t  POS_LC_In_Silico") as pbar:
        POS_LC_In_Silico_temp = POS_LC[POS_LC['PREDICTED'] == "true"]
        pbar.update(len(POS_LC_In_Silico_temp))

    with tqdm(total=len(POS_GC), unit=" spectrums", colour="green", desc="\t  POS_GC_In_Silico") as pbar:
        POS_GC_In_Silico_temp = POS_GC[POS_GC['PREDICTED'] == "true"]
        pbar.update(len(POS_GC_In_Silico_temp))

    with tqdm(total=len(NEG_LC), unit=" spectrums", colour="green", desc="\t  NEG_LC_In_Silico") as pbar:
        NEG_LC_In_Silico_temp = NEG_LC[NEG_LC['PREDICTED'] == "true"]
        pbar.update(len(NEG_LC_In_Silico_temp))

    with tqdm(total=len(NEG_GC), unit=" spectrums", colour="green", desc="\t  NEG_GC_In_Silico") as pbar:
        NEG_GC_In_Silico_temp = NEG_GC[NEG_GC['PREDICTED'] == "true"]
        pbar.update(len(NEG_GC_In_Silico_temp))

    # Séparation pour les lignes contenant "False"
    with tqdm(total=len(POS_LC), unit=" spectrums", colour="green", desc="\t  POS_LC_Exp") as pbar:
        POS_LC_temp = POS_LC[POS_LC['PREDICTED'] == "false"]
        pbar.update(len(POS_LC_temp))

    with tqdm(total=len(POS_GC), unit=" spectrums", colour="green", desc="\t  POS_GC_Exp") as pbar:
        POS_GC_temp = POS_GC[POS_GC['PREDICTED'] == "false"]
        pbar.update(len(POS_GC_temp))

    with tqdm(total=len(NEG_LC), unit=" spectrums", colour="green", desc="\t  NEG_LC_Exp") as pbar:
        NEG_LC_temp = NEG_LC[NEG_LC['PREDICTED'] == "false"]
        pbar.update(len(NEG_LC_temp))

    with tqdm(total=len(NEG_GC), unit=" spectrums", colour="green", desc="\t  NEG_GC_Exp") as pbar:
        NEG_GC_temp = NEG_GC[NEG_GC['PREDICTED'] == "false"]
        pbar.update(len(NEG_GC_temp))

    return POS_LC_temp,POS_LC_In_Silico_temp,POS_GC_temp,POS_GC_In_Silico_temp,NEG_LC_temp,NEG_LC_In_Silico_temp,NEG_GC_temp,NEG_GC_In_Silico_temp