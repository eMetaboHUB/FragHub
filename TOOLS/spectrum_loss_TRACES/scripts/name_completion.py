from tqdm import tqdm

def names_completion(CONCATENATE_DF):
    """
    Perform name completion for a given DataFrame.
    :param CONCATENATE_DF: The DataFrame containing the 'NAME' column to perform name completion on.
    :return: The modified DataFrame with name completion applied.
    """

    # A progress bar from tqdm library is initialized. The total progress is the length of the DataFrame provided and color is green.
    # Each task done is represented by "row" and it provides a description of "updating names".
    tqdm.pandas(total=len(CONCATENATE_DF), colour="green", unit=" row", desc="{:>70}".format("updating names"))

    # The DataFrame is grouped by the 'INCHI' column. Then for each group, a transformation is applied on the 'NAME' column of the DataFrame.
    # The transformation function uses fillna method to replace the 'NaN' values in the group.
    # It replaces 'NaN' values with the first non 'NaN' value in the group if there is any non 'NaN' value.
    # If there are only 'NaN' values in the group, it would replace 'NaN' with an empty string.
    # The progress_transform method applies the function to each group and also updates the progress bar.
    CONCATENATE_DF['NAME'] = CONCATENATE_DF.groupby('INCHI')['NAME'].progress_transform(lambda group: group.fillna(group.dropna().iloc[0] if group.dropna().size > 0 else ''))

    # The DataFrame with updated 'NAME' column is returned.
    return CONCATENATE_DF
