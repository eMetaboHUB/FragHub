from tqdm import tqdm

def names_completion(CONCATENATE_DF):
    """
    Perform name completion for a given DataFrame.
    :param CONCATENATE_DF: The DataFrame containing the 'NAME' column to perform name completion on.
    :return: The modified DataFrame with name completion applied.
    """
    # Count the number of unique groups (based on 'INCHI')
    num_groups = CONCATENATE_DF['INCHI'].nunique()

    # Initialize tqdm progress bar with total as number of unique groups
    pbar = tqdm(total=num_groups, colour="green", unit=" group", desc="{:>70}".format("updating names"))

    # Define the function to apply to each group
    def fill_group(group):
        filled_group = group.fillna(group.dropna().iloc[0] if group.dropna().size > 0 else '')
        pbar.update()  # update the progress bar
        return filled_group

    # Apply the function to each group
    CONCATENATE_DF['NAME'] = CONCATENATE_DF.groupby('INCHI')['NAME'].transform(fill_group)

    pbar.close()  # close the progress bar

    # Return the updated DataFrame
    return CONCATENATE_DF
