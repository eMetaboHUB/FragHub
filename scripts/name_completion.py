from tqdm import tqdm

def names_completion(CONCATENATE_DF):
    """
    Perform name completion for a given DataFrame.

    :param CONCATENATE_DF: The DataFrame containing the 'NAME' column to perform name completion on.
    :return: The modified DataFrame with name completion applied.

    """
    # Définir le nom de la barre de progression
    tqdm.pandas(total=len(CONCATENATE_DF), colour="green", unit=" row", desc="{:>80}".format("updating names"))

    # Appliquer la transformation par groupe avec une barre de progression
    CONCATENATE_DF['NAME'] = CONCATENATE_DF.groupby('INCHI')['NAME'].progress_transform(lambda group: group.fillna(group.dropna().iloc[0] if group.dropna().size > 0 else ''))

    return CONCATENATE_DF