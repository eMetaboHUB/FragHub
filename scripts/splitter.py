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


def split_LC_GC(POS, NEG, progress_callback=None, total_items_callback=None, prefix_callback=None,
                item_type_callback=None):
    """
    This function separates the given positive and negative DataFrames into LC and GC spectrums.

    :param POS: DataFrame containing positive spectrums.
    :param NEG: DataFrame containing negative spectrums.
    :param progress_callback: Callable function to update progress information.
    :param total_items_callback: Callable function to define the total number of items to process.
    :param prefix_callback: Callable function to provide the current task context.
    :param item_type_callback: Callable function to define the processed item type.
    :return: Four DataFrames (POS_LC, POS_GC, NEG_LC, NEG_GC).
    """

    # Initial callbacks to describe the process
    if prefix_callback:
        prefix_callback("Splitting LC/GC")
    if item_type_callback:
        item_type_callback("rows")

    # Total rows to process (Positive + Negative)
    total_rows = len(POS) + len(NEG)
    if total_items_callback:
        total_items_callback(total_rows, 0)  # Report total and initialize completed items to 0

    # Splitting positive GC spectrums
    POS_GC = POS[POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
    if progress_callback:
        progress_callback(len(POS_GC))  # Update progress after processing positive GC

    # Splitting positive LC spectrums
    POS_LC = POS[~POS['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
    if progress_callback:
        progress_callback(len(POS))  # Update progress after processing all positive spectrums

    # Splitting negative GC spectrums
    NEG_GC = NEG[NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
    if progress_callback:
        progress_callback(len(POS) + len(NEG_GC))  # Update progress after processing negative GC

    # Splitting negative LC spectrums
    NEG_LC = NEG[~NEG['INSTRUMENTTYPE'].str.contains('GC|EI', case=False)]
    if progress_callback:
        progress_callback(total_rows)  # Indicate that all rows have been processed

    # Returning separated DataFrames
    return POS_LC, POS_GC, NEG_LC, NEG_GC


def exp_in_silico_splitter(POS_LC, POS_GC, NEG_LC, NEG_GC, progress_callback=None, total_items_callback=None,
                           prefix_callback=None, item_type_callback=None):
    """
    Function to split the provided datasets into various categories based on 'PREDICTED' column value being "true" or "false".

    :param POS_LC: The dataframe containing the positive LC data.
    :param POS_GC: The dataframe containing the positive GC data.
    :param NEG_LC: The dataframe containing the negative LC data.
    :param NEG_GC: The dataframe containing the negative GC data.
    :param progress_callback: Callable function to update progress.
    :param total_items_callback: Callable function to set total number of items to process.
    :param prefix_callback: Callable function to provide current task context.
    :param item_type_callback: Callable function to define the type of items being processed.
    :return: Tuple containing separate dataframes for positive/negative LC/GC in silico and experimental analyses.
    """

    # Initial callbacks for task setup
    if prefix_callback:
        prefix_callback("Splitting Experimental/InSilico")
    if item_type_callback:
        item_type_callback("rows")

    total_rows = len(POS_LC) + len(POS_GC) + len(NEG_LC) + len(NEG_GC)
    if total_items_callback:
        total_items_callback(total_rows, 0)  # Initialize the total items

    # Helper function for splitting
    def split_dataframe(df, field, value):
        return df[df[field] == value]

    # Splitting POS_LC (In Silico and Experimental)
    POS_LC_In_Silico_temp = split_dataframe(POS_LC, 'PREDICTED', "true")
    if progress_callback:
        progress_callback(len(POS_LC_In_Silico_temp))

    POS_LC_temp = split_dataframe(POS_LC, 'PREDICTED', "false")
    if progress_callback:
        progress_callback(len(POS_LC))

    # Splitting POS_GC (In Silico and Experimental)
    POS_GC_In_Silico_temp = split_dataframe(POS_GC, 'PREDICTED', "true")
    if progress_callback:
        progress_callback(len(POS_LC) + len(POS_GC_In_Silico_temp))

    POS_GC_temp = split_dataframe(POS_GC, 'PREDICTED', "false")
    if progress_callback:
        progress_callback(len(POS_LC) + len(POS_GC))

    # Splitting NEG_LC (In Silico and Experimental)
    NEG_LC_In_Silico_temp = split_dataframe(NEG_LC, 'PREDICTED', "true")
    if progress_callback:
        progress_callback(len(POS_LC) + len(POS_GC) + len(NEG_LC_In_Silico_temp))

    NEG_LC_temp = split_dataframe(NEG_LC, 'PREDICTED', "false")
    if progress_callback:
        progress_callback(len(POS_LC) + len(POS_GC) + len(NEG_LC))

    # Splitting NEG_GC (In Silico and Experimental)
    NEG_GC_In_Silico_temp = split_dataframe(NEG_GC, 'PREDICTED', "true")
    if progress_callback:
        progress_callback(len(POS_LC) + len(POS_GC) + len(NEG_LC) + len(NEG_GC_In_Silico_temp))

    NEG_GC_temp = split_dataframe(NEG_GC, 'PREDICTED', "false")
    if progress_callback:
        progress_callback(total_rows)  # Progress complete

    # Returning all datasets
    return (POS_LC_temp, POS_LC_In_Silico_temp, POS_GC_temp, POS_GC_In_Silico_temp,
            NEG_LC_temp, NEG_LC_In_Silico_temp, NEG_GC_temp, NEG_GC_In_Silico_temp)

