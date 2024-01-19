from tqdm import tqdm

def names_completion(CONCATENATE_DF):
    """
    Perform name completion for a given DataFrame.

    :param CONCATENATE_DF: The DataFrame containing the 'NAME' column to perform name completion on.
    :return: The modified DataFrame with name completion applied.

    """
    # Set progress bar
    tqdm.pandas(total=len(CONCATENATE_DF), colour="green", unit=" row", desc="{:>70}".format("updating names"))

    # Apply transformation by group with a progress bar
    CONCATENATE_DF['NAME'] = CONCATENATE_DF.groupby('INCHI')['NAME'].progress_transform(lambda group: group.fillna(group.dropna().iloc[0] if group.dropna().size > 0 else ''))

    return CONCATENATE_DF