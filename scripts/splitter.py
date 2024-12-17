from tqdm import tqdm
import re

def split_pos_neg(CONCATENATE_DF, progress_callback=None, total_items_callback=None, prefix_callback=None,
                  item_type_callback=None):
    """
    This function splits the given DataFrame into two DataFrames based on the value of the `IONMODE` column
    :param CONCATENATE_DF: Initial DataFrame containing mixed type data
    :param progress_callback: Callable function to report progress.
    :param total_items_callback: Callable function to set the total number of items to be processed.
    :param prefix_callback: Callable function to describe the current task.
    :param item_type_callback: Callable function to define the type of item being processed.
    :return: Two DataFrames, one for positive and one for negative 'IONMODE' values
    """

    # Initial tasks description callbacks
    if prefix_callback:
        prefix_callback("Splitting POS/NEG")
    if item_type_callback:
        item_type_callback("rows")

    # Set total items at the beginning
    if total_items_callback:
        total_items_callback(len(CONCATENATE_DF), 0)

    # Splitting the DataFrame where 'IONMODE' = 'positive'
    POS = CONCATENATE_DF[CONCATENATE_DF['IONMODE'] == 'positive']
    if progress_callback:
        progress_callback(len(POS))  # Update progress with the size of the split positive DataFrame

    # Splitting the DataFrame where 'IONMODE' = 'negative'
    NEG = CONCATENATE_DF[CONCATENATE_DF['IONMODE'] == 'negative']
    if progress_callback:
        progress_callback(len(CONCATENATE_DF))  # Update progress to the total size, indicating completion

    # Return both split DataFrames
    return POS, NEG


def split_LC_GC(POS, NEG):
    """
    :param POS: DataFrame containing positive spectrums
    :param NEG: DataFrame containing negative spectrums
    :return: Four DataFrames containing positive LC spectrums, positive GC spectrums, negative LC spectrums, and negative GC spectrums
    """
    # Creating a progress bar for positive GC spectrums separation
    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>70}".format("POS_GC")) as pbar:
        # Extracting rows in POS DataFrame where 'INSTRUMENTTYPE' contains either 'GC' or 'EI'. Case-insensitive
        POS_GC = POS[POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
        # Updating the progress bar with the number of rows in POS_GC
        pbar.update(len(POS_GC))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Creating a progress bar for positive LC spectrums separation
    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>70}".format("POS_LC")) as pbar:
        # Extracting rows in POS DataFrame where 'INSTRUMENTTYPE' does not contain either 'GC' or 'EI'. Case-insensitive
        POS_LC = POS[~POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
        # Updating the progress bar with the number of rows in POS_LC
        pbar.update(len(POS_LC))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Creating a progress bar for negative GC spectrums separation
    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>70}".format("NEG_GC")) as pbar:
        # Extracting rows in NEG DataFrame where 'INSTRUMENTTYPE' contains either 'GC' or 'EI'. Case-insensitive
        NEG_GC = NEG[NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
        # Updating the progress bar with the number of rows in NEG_GC
        pbar.update(len(NEG_GC))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Creating a progress bar for negative LC spectrums separation
    with tqdm(total=len(POS), unit=" spectrums", colour="green", desc="{:>70}".format("NEG_LC")) as pbar:
        # Extracting rows in NEG DataFrame where 'INSTRUMENTTYPE' does not contain either 'GC' or 'EI'. Case-insensitive
        NEG_LC = NEG[~NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
        # Updating the progress bar with the number of rows in NEG_LC
        pbar.update(len(NEG_LC))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Returning separated POS and NEG DataFrames for LC and GC spectrums
    return POS_LC, POS_GC, NEG_LC, NEG_GC

def exp_in_silico_splitter(POS_LC, POS_GC, NEG_LC, NEG_GC):
    """
    Function to split the provided datasets into various categories based on 'PREDICTED' column value being "true" or "false"

    :param POS_LC: The data frame containing the positive LC data
    :param POS_GC: The data frame containing the positive GC data
    :param NEG_LC: The data frame containing the negative LC data
    :param NEG_GC: The data frame containing the negative GC data
    :return: A tuple containing separate data frames for positive LC in silico, positive GC in silico,
        negative LC in silico, negative GC in silico, positive LC experimental, positive GC experimental,
        negative LC experimental, and negative GC experimental.
    """

    # Create and update progress bar for POS_LC_In_Silico
    # Get rows where PREDICTED column value is "true"
    with tqdm(total=len(POS_LC), unit=" spectrums", colour="green", desc="{:>70}".format("POS_LC_In_Silico")) as pbar:
        POS_LC_In_Silico_temp = POS_LC[POS_LC['PREDICTED'] == "true"]
        pbar.update(len(POS_LC_In_Silico_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Repeat similar process for POS_GC_In_Silico
    with tqdm(total=len(POS_GC), unit=" spectrums", colour="green", desc="{:>70}".format("POS_GC_In_Silico")) as pbar:
        POS_GC_In_Silico_temp = POS_GC[POS_GC['PREDICTED'] == "true"]
        pbar.update(len(POS_GC_In_Silico_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Continue for NEG_LC_In_Silico
    with tqdm(total=len(NEG_LC), unit=" spectrums", colour="green", desc="{:>70}".format("NEG_LC_In_Silico")) as pbar:
        NEG_LC_In_Silico_temp = NEG_LC[NEG_LC['PREDICTED'] == "true"]
        pbar.update(len(NEG_LC_In_Silico_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # And for NEG_GC_In_Silico
    with tqdm(total=len(NEG_GC), unit=" spectrums", colour="green", desc="{:>70}".format("NEG_GC_In_Silico")) as pbar:
        NEG_GC_In_Silico_temp = NEG_GC[NEG_GC['PREDICTED'] == "true"]
        pbar.update(len(NEG_GC_In_Silico_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # The next set of progress bars are for rows where PREDICTED column value is "false"

    # Proceed with POS_LC_Exp
    with tqdm(total=len(POS_LC), unit=" spectrums", colour="green", desc="{:>70}".format("POS_LC_Exp")) as pbar:
        POS_LC_temp = POS_LC[POS_LC['PREDICTED'] == "false"]
        pbar.update(len(POS_LC_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Then with POS_GC_Exp
    with tqdm(total=len(POS_GC), unit=" spectrums", colour="green", desc="{:>70}".format("POS_GC_Exp")) as pbar:
        POS_GC_temp = POS_GC[POS_GC['PREDICTED'] == "false"]
        pbar.update(len(POS_GC_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Continue with NEG_LC_Exp
    with tqdm(total=len(NEG_LC), unit=" spectrums", colour="green", desc="{:>70}".format("NEG_LC_Exp")) as pbar:
        NEG_LC_temp = NEG_LC[NEG_LC['PREDICTED'] == "false"]
        pbar.update(len(NEG_LC_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Finally, proceed with NEG_GC_Exp
    with tqdm(total=len(NEG_GC), unit=" spectrums", colour="green", desc="{:>70}".format("NEG_GC_Exp")) as pbar:
        NEG_GC_temp = NEG_GC[NEG_GC['PREDICTED'] == "false"]
        pbar.update(len(NEG_GC_temp))
        pbar.n = pbar.total
        pbar.refresh()
        pbar.close()

    # Return a tuple of all newly created dataframes.
    return POS_LC_temp, POS_LC_In_Silico_temp, POS_GC_temp, POS_GC_In_Silico_temp, NEG_LC_temp, NEG_LC_In_Silico_temp, NEG_GC_temp, NEG_GC_In_Silico_temp
