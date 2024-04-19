from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import time
import os
import re

def write_msp(spectrum_list, filename, mode, update, profile_name):
    """
    Write MSP file.
    :param spectrum_list: A list of spectra to be written to the MSP file.
    :param filename: The name of the output file.
    :param mode: The mode of writing (either 'w' for write or 'a' for append).
    :param update: Boolean value indicating whether to update the existing file or create a new one.
    :param profile_name: The name of the profile being used.
    :return: None
    """

    # Construct the output path
    output_file_path = os.path.join(f"../OUTPUT/{profile_name}/MSP/{mode}", filename)

    # Creating a progress bar for displaying the status of the operation.
    with tqdm(total=len(spectrum_list), unit=" row", colour="green", desc="{:>70}".format(f"writting {filename}")) as pbar:

        # When update is not required (writing in a new file or overwriting the existing file).
        if not update:
            # Open the file in write mode
            with open(output_file_path, 'w', encoding="UTF-8") as f:

                # Iterate through each spectrum in the list
                for spectrum in spectrum_list:
                    try:
                        # Write the spectrum to the file
                        f.write(spectrum)

                        # Add two newline characters after each spectrum
                        f.write("\n\n")
                    except:
                        # In case of any exception, skip the current spectrum
                        continue

                    # Update the progress bar
                    pbar.update()
                pbar.close()
        else:
            # When update is required (appending to the existing file).
            with open(output_file_path, 'a', encoding="UTF-8") as f:
                for spectrum in spectrum_list:
                    try:
                        f.write(spectrum)
                        f.write("\n\n")
                    except:
                        continue

                    pbar.update()
                pbar.close()

def writting_msp(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico, profile_name, update=False):
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
    write_msp(POS_LC, "POS_LC.msp", "POS", update, profile_name)
    # delete the positive LC data to free up memory
    del POS_LC
    # repeat the process for other types of data
    time.sleep(0.1)
    write_msp(POS_LC_insilico, "POS_LC_insilico.msp", "POS", update, profile_name)
    del POS_LC_insilico
    time.sleep(0.1)
    write_msp(POS_GC, "POS_GC.msp", "POS", update, profile_name)
    del POS_GC
    time.sleep(0.1)
    write_msp(POS_GC_insilico, "POS_GC_insilico.msp", "POS", update, profile_name)
    del POS_GC_insilico
    time.sleep(0.1)
    write_msp(NEG_LC, "NEG_LC.msp", "NEG", update, profile_name)
    del NEG_LC
    time.sleep(0.1)
    write_msp(NEG_LC_insilico, "NEG_LC_insilico.msp", "NEG", update, profile_name)
    del NEG_LC_insilico
    time.sleep(0.1)
    write_msp(NEG_GC, "NEG_GC.msp", "NEG", update, profile_name)
    del NEG_GC
    time.sleep(0.1)
    write_msp(NEG_GC_insilico, "NEG_GC_insilico.msp", "NEG", update, profile_name)
    del NEG_GC_insilico

def write_csv(df, filename, mode, update, first_run, profile_name):
    """
    This function writes a pandas DataFrame to a CSV file in chunks.
    :param df: DataFrame - the data to be written to the CSV file
    :param filename: str - the name of CSV file to which the data will be written
    :param mode: str - the mode in which the file should be opened ('w' - write, 'a' - append)
    :param update: Bool - a flag to indicate whether the file should be updated or created from scratch
    :param first_run: Bool - indicates whether it's the first run (headers needs to be written)
    :param profile_name: str - a name related to the CSV file for organizing the output.
    :return: None.
    """
    # Placeholder string for the directory output path
    output_file_path = os.path.join(f"../OUTPUT/{profile_name}/CSV/{mode}", filename)

    # The chunk size is limit for each write operation
    chunk_size = 5000

    # Calculate the number of chunks that will be required to write the entire dataframe
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))

    # tqdm progress bar initialization and settings
    with tqdm(total=num_chunks, unit=" row", colour="green", desc="{:>70}".format(f"writting {filename}")) as pbar:

        # Split the dataframe into chunks and write each chunk separately
        for start in range(0, df.shape[0], chunk_size):
            # Slice the dataframe for writing
            df_slice = df[start:start + chunk_size]

            # If it is the first chunk in a first run, write with headers
            if start == 0 and first_run:
                df_slice.to_csv(output_file_path, mode='w', sep=";", quotechar='"', encoding="UTF-8", index=False)

            # If it is not the first chunk, write without headers (append)
            else:
                df_slice.to_csv(output_file_path, mode='a', header=False, index=False, sep=";", quotechar='"', encoding="UTF-8")

            # Update the progress bar
            pbar.update()
        pbar.close()

def writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, first_run, profile_name, update=False):
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

    # Pause execution for 0.1 seconds
    time.sleep(0.1)
    # Call write_csv for Positive LC data and write it to POS_LC.csv file
    write_csv(POS_LC_df, "POS_LC.csv", "POS", update, first_run, profile_name)
    # Clear memory used by POS_LC_df dataframe
    del POS_LC_df

    # Similar comments can be written for the rest of the calls and deletions below
    # Just replace "POS_LC_df" and "POS_LC.csv" with appropriate values

    time.sleep(0.1)
    write_csv(POS_GC_df, "POS_GC.csv", "POS", update, first_run, profile_name)
    del POS_GC_df

    time.sleep(0.1)
    write_csv(NEG_LC_df, "NEG_LC.csv", "NEG", update, first_run, profile_name)
    del NEG_LC_df

    time.sleep(0.1)
    write_csv(NEG_GC_df, "NEG_GC.csv", "NEG", update, first_run, profile_name)
    del NEG_GC_df

    time.sleep(0.1)
    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico.csv", "POS", update, first_run, profile_name)
    del POS_LC_df_insilico

    time.sleep(0.1)
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico.csv", "POS", update, first_run, profile_name)
    del POS_GC_df_insilico

    time.sleep(0.1)
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico.csv", "NEG", update, first_run, profile_name)
    del NEG_LC_df_insilico

    time.sleep(0.1)
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico.csv", "NEG", update, first_run, profile_name)
    del NEG_GC_df_insilico

def write_json(df, filename, mode, profile_name):
    """
    Write DataFrame to a JSON file.
    :param df: The DataFrame to be written.
    :type df: pandas.DataFrame
    :param filename: The name of the output file.
    :type filename: str
    :param mode: The mode for opening the output file. e.g., 'w', 'a'.
    :type mode: str
    :param profile_name: The name of the profile.
    :type profile_name: str
    :return: None
    """

    # Define the output file path for the JSON file.
    output_file_path = os.path.join(f"../OUTPUT/{profile_name}/JSON/{mode}", filename)

    # Convert the DataFrame to a list of dict records.
    json_records = df.to_dict('records')

    # Open the file with write mode.
    with open(output_file_path, 'w') as f:

        # Create a progress bar for writing records to file
        pbar = tqdm(total=len(json_records), unit=" row", colour="green", desc="{:>70}".format(f"writing {filename}"))

        # Write an opening square bracket to start the JSON array.
        f.write('[\n')

        # Loop through each record in the DataFrame.
        for i, record in enumerate(json_records):

            # Write a tab character for indentation.
            f.write('\t')

            # Dump the current record to the file as a JSON formatted string.
            json.dump(record, f)

            # Add a comma after each record except the last one.
            if i < len(json_records) - 1:
                f.write(',\n')

            # Update the progress bar
            pbar.update()
        pbar.close()

        # Write a closing square bracket to end the JSON array.
        f.write('\n]')

def writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, profile_name):
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
    write_json(POS_LC_df, "POS_LC.json", "POS", profile_name)
    # Deleting the DataFrame from memory as it's no longer needed
    del POS_LC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive GC DataFrame to a JSON file named "POS_GC.json"
    write_json(POS_GC_df, "POS_GC.json", "POS", profile_name)
    del POS_GC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative LC DataFrame to a JSON file named "NEG_LC.json"
    write_json(NEG_LC_df, "NEG_LC.json", "NEG", profile_name)
    del NEG_LC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative GC DataFrame to a JSON file named "NEG_GC.json"
    write_json(NEG_GC_df, "NEG_GC.json", "NEG", profile_name)
    del NEG_GC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive In Silico LC DataFrame to a JSON file named "POS_LC_In_Silico.json"
    write_json(POS_LC_df_insilico, "POS_LC_In_Silico.json", "POS", profile_name)
    del POS_LC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive In Silico GC DataFrame to a JSON file named "POS_GC_In_Silico.json"
    write_json(POS_GC_df_insilico, "POS_GC_In_Silico.json", "POS", profile_name)
    del POS_GC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative In Silico LC DataFrame to a JSON file named "NEG_LC_In_Silico.json"
    write_json(NEG_LC_df_insilico, "NEG_LC_In_Silico.json", "NEG", profile_name)
    del NEG_LC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative In Silico GC DataFrame to a JSON file named "NEG_GC_In_Silico.json"
    write_json(NEG_GC_df_insilico, "NEG_GC_In_Silico.json", "NEG", profile_name)

    # Deleting the DataFrame from memory as it's no longer needed
    del NEG_GC_df_insilico
