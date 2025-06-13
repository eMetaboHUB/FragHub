from scripts.set_projects import parameters_dict
import pandas as pd
import numpy as np
import json
import time
import os
import re

def write_msp(spectrum_list, filename, mode, update, output_directory, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Write an MSP file.
    :param spectrum_list: A list of spectra to be written to the MSP file.
    :param filename: The name of the output file.
    :param mode: The mode of writing (either 'w' for writing or 'a' for appending).
    :param update: Boolean value indicating whether to update the existing file or create a new one.
    :param progress_callback: Callback function to report progress.
    :param total_items_callback: Callback function to set the total number of items.
    :param prefix_callback: Callback function to provide the context/task prefix.
    :param item_type_callback: Callback function to specify the type of items (e.g., "spectra").
    :return: None
    """

    # Handle the output path
    output_file_path = f"{output_directory}/MSP/{mode}/{filename}"

    # Set the context of the operation using the prefix_callback
    if prefix_callback:
        prefix_callback(f"Writing {filename} to MSP:")

    # Specify the type of items with the item_type_callback
    if item_type_callback:
        item_type_callback("spectra")

    # Initialize the total items using total_items_callback
    total_spectra = len(spectrum_list)
    if total_items_callback:
        total_items_callback(total_spectra, 0)  # Initialize with 0 progress

    # Open the file for writing spectra
    file_mode = 'a' if update else 'w'
    with open(output_file_path, file_mode, encoding="UTF-8") as f:
        # Iterate through each spectrum in the list
        for index, spectrum in enumerate(spectrum_list):
            try:
                # Write the spectrum to the file
                f.write(spectrum)
                f.write("\n\n")

                # Update progress using the callback after each write
                if progress_callback:
                    progress_callback(index + 1)  # Signal progress (current line)
            except Exception as e:
                # If an error occurs, skip the spectrum but log the issue if necessary
                print(f"Error writing spectrum {index}: {e}")
                continue


def writting_msp(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico, output_directory, update=False, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    This function writes MSP files for positive and negative LC and GC data. It takes in various inputs for different types of data, a profile_name which is the name of the profile and an optional update flag which defaults to False when not provided.

    Args:
      POS_LC: This parameter represents the positive LC data.
      POS_LC_insilico: This parameter represents the insilico data for positive LC.
      POS_GC: This parameter represents the positive GC data.
      POS_GC_insilico: This parameter represents the insilico data for positive GC.
      NEG_LC: This parameter represents the negative LC data.
      NEG_LC_insilico: This parameter represents the insilico data for negative LC.
      NEG_GC: This parameter represents the negative GC data.
      NEG_GC_insilico: This parameter represents the insilico data for negative GC.
      profile_name: This is the name of the profile.
      update: This is an optional parameter which is set to False by default. If set to True, it results in an update operation.

    Returns: None
    """

    # sleep for 0.1 second to avoid too many computations
    time.sleep(0.1)
    # write the positive LC data to the "POS_LC.msp" file
    write_msp(POS_LC, "POS_LC.msp", "POS", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # delete the positive LC data to free up memory
    del POS_LC
    # repeat the process for other types of data
    time.sleep(0.1)
    write_msp(POS_LC_insilico, "POS_LC_insilico.msp", "POS", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_LC_insilico
    time.sleep(0.1)
    write_msp(POS_GC, "POS_GC.msp", "POS", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC
    time.sleep(0.1)
    write_msp(POS_GC_insilico, "POS_GC_insilico.msp", "POS", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_insilico
    time.sleep(0.1)
    write_msp(NEG_LC, "NEG_LC.msp", "NEG", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC
    time.sleep(0.1)
    write_msp(NEG_LC_insilico, "NEG_LC_insilico.msp", "NEG", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_insilico
    time.sleep(0.1)
    write_msp(NEG_GC, "NEG_GC.msp", "NEG", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC
    time.sleep(0.1)
    write_msp(NEG_GC_insilico, "NEG_GC_insilico.msp", "NEG", update, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_insilico


def write_csv(df, filename, mode, update, first_run, output_directory, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    This function writes a pandas DataFrame to a CSV file in chunks.
    :param df: DataFrame - the data to be written to the CSV file.
    :param filename: str - the name of the CSV file to write the data to.
    :param mode: str - the mode in which the file should be opened ('w' - write, 'a' - append).
    :param update: Bool - a flag to indicate whether the file is updated or created from scratch.
    :param first_run: Bool - indicates whether it's the first run (write the headers if True).
    :param profile_name: str - a name related to the CSV file for output organization.
    :param progress_callback: Callable, function to track progress of row writing.
    :param total_items_callback: Callable, function to set the total number of rows.
    :param prefix_callback: Callable, function to indicate the current task's description context.
    :param item_type_callback: Callable, function to indicate the type of elements being processed.
    :return: None.
    """
    # Replace newline characters with semicolons in the PEAKS_LIST column
    if 'PEAKS_LIST' in df.columns:
        df['PEAKS_LIST'] = df['PEAKS_LIST'].str.replace('\n', ';', regex=False)

    # Construct the file path dynamically
    output_file_path = f"{output_directory}/CSV/{mode}/{filename}"

    # Define chunk size for writing DataFrame in parts
    chunk_size = 5000

    # Calculate the total number of chunks
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))

    # Set up task description
    if prefix_callback:
        prefix_callback(f"Writing {filename} to CSV:")

    # Indicate item type
    if item_type_callback:
        item_type_callback("rows")

    # Define total items expected
    total_rows = len(df)
    if total_items_callback:
        total_items_callback(total_rows, 0)  # Initialize total and completed count

    # Process the DataFrame in chunks and write each chunk to the file
    for chunk_index, start in enumerate(range(0, total_rows, chunk_size)):
        # Select the slice of the DataFrame for the current chunk
        df_slice = df.iloc[start:start + chunk_size]

        # Write headers only on the first chunk and during the first run
        if start == 0 and first_run:
            df_slice.to_csv(output_file_path, mode='w', sep="\t", quotechar='"', encoding="UTF-8", index=False)
        else:
            # Append the chunk without writing headers
            df_slice.to_csv(output_file_path, mode='a', sep="\t", quotechar='"', encoding="UTF-8", index=False,
                            header=False)

        # Notify progress of written chunks (via progress_callback)
        if progress_callback:
            processed_rows = min((chunk_index + 1) * chunk_size, total_rows)
            progress_callback(processed_rows)  # Update with processed row count


def writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, first_run, output_directory, update=False, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Write data to CSV files.

    :param POS_LC_df: Pandas DataFrame containing Positive LC data
    :param POS_GC_df: Pandas DataFrame containing Positive GC data
    :param NEG_LC_df: Pandas DataFrame containing Negative LC data
    :param NEG_GC_df: Pandas DataFrame containing Negative GC data
    :param POS_LC_df_insilico: Pandas DataFrame containing Positive LC data (In Silico)
    :param POS_GC_df_insilico: Pandas DataFrame containing Positive GC data (In Silico)
    :param NEG_LC_df_insilico: Pandas DataFrame containing Negative LC data (In Silico)
    :param NEG_GC_df_insilico: Pandas DataFrame containing Negative GC data (In Silico)
    :param first_run: Boolean indicating if it's the first run
    :param profile_name: Name of the profile
    :param update: Boolean indicating if it's an update
    :return: None
    """

    time.sleep(0.1)
    # Call write_csv for Positive LC data and write it to POS_LC.csv file
    write_csv(POS_LC_df, "POS_LC.csv", "POS", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Clear memory used by POS_LC_df dataframe
    del POS_LC_df

    time.sleep(0.1)
    write_csv(POS_GC_df, "POS_GC.csv", "POS", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df

    time.sleep(0.1)
    write_csv(NEG_LC_df, "NEG_LC.csv", "NEG", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df

    time.sleep(0.1)
    write_csv(NEG_GC_df, "NEG_GC.csv", "NEG", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_df

    time.sleep(0.1)
    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico.csv", "POS", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_LC_df_insilico

    time.sleep(0.1)
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico.csv", "POS", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df_insilico

    time.sleep(0.1)
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico.csv", "NEG", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df_insilico

    time.sleep(0.1)
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico.csv", "NEG", update, first_run, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_df_insilico


def write_json(df, filename, mode, output_directory, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Write a DataFrame to a JSON file.
    :param df: The DataFrame to be written.
    :type df: pandas.DataFrame
    :param filename: The name of the output file.
    :type filename: str
    :param mode: The mode for opening the output file, e.g., 'w', 'a'.
    :type mode: str
    :param progress_callback: Callback to report progress.
    :param total_items_callback: Callback to set the total number of items.
    :param prefix_callback: Callback to define a prefix or task description.
    :param item_type_callback: Callback to specify the type of items being processed.
    :return: None
    """

    # Define the output file path
    output_file_path = f"{output_directory}/JSON/{mode}/{filename}"

    # Convert the DataFrame to a list of dictionaries
    json_records = df.to_dict('records')

    # If a prefix_callback is defined, set the current task context
    if prefix_callback:
        prefix_callback(f"Writing {filename} to JSON:")

    # Set the type of items using the item_type_callback
    if item_type_callback:
        item_type_callback("rows")

    # Declare the total number of records to process
    total_records = len(json_records)
    if total_items_callback:
        total_items_callback(total_records, 0)  # Initialize with 0 processed items

    # Open the file to write
    with open(output_file_path, 'w') as f:
        # Write the opening bracket to start the JSON array
        f.write('[\n')

        # Iterate through the converted records
        for i, record in enumerate(json_records):
            # Write a tab for indentation
            f.write('\t')

            # Write the dictionary as a JSON string
            json.dump(record, f)

            # Add a comma at the end, except for the last element
            if i < total_records - 1:
                f.write(',\n')

            # If a progress_callback is defined, call it to update the progress
            if progress_callback:
                progress_callback(i + 1)  # Update with the current index (starting from 1)

        # Write the closing bracket to end the JSON array
        f.write('\n]')


def writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, output_directory, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Write dataframes to JSON files.
    :param POS_LC_df: Positive LC dataframe
    :param POS_GC_df: Positive GC dataframe
    :param NEG_LC_df: Negative LC dataframe
    :param NEG_GC_df: Negative GC dataframe
    :param POS_LC_df_insilico: Positive LC In Silico dataframe
    :param POS_GC_df_insilico: Positive GC In Silico dataframe
    :param NEG_LC_df_insilico: Negative LC In Silico dataframe
    :param NEG_GC_df_insilico: Negative GC In Silico dataframe
    :param profile_name: Name of the profile
    :return: None
    """
    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive LC DataFrame to a JSON file named "POS_LC.json"
    write_json(POS_LC_df, "POS_LC.json", "POS", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Deleting the DataFrame from memory as it's no longer needed
    del POS_LC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive GC DataFrame to a JSON file named "POS_GC.json"
    write_json(POS_GC_df, "POS_GC.json", "POS", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative LC DataFrame to a JSON file named "NEG_LC.json"
    write_json(NEG_LC_df, "NEG_LC.json", "NEG", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative GC DataFrame to a JSON file named "NEG_GC.json"
    write_json(NEG_GC_df, "NEG_GC.json", "NEG", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive In Silico LC DataFrame to a JSON file named "POS_LC_In_Silico.json"
    write_json(POS_LC_df_insilico, "POS_LC_In_Silico.json", "POS", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_LC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive In Silico GC DataFrame to a JSON file named "POS_GC_In_Silico.json"
    write_json(POS_GC_df_insilico, "POS_GC_In_Silico.json", "POS", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative In Silico LC DataFrame to a JSON file named "NEG_LC_In_Silico.json"
    write_json(NEG_LC_df_insilico, "NEG_LC_In_Silico.json", "NEG", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative In Silico GC DataFrame to a JSON file named "NEG_GC_In_Silico.json"
    write_json(NEG_GC_df_insilico, "NEG_GC_In_Silico.json", "NEG", output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Deleting the DataFrame from memory as it's no longer needed
    del NEG_GC_df_insilico
