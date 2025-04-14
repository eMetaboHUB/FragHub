from scripts.GUI.utils.global_vars import parameters_dict
import json
import os

def remove_files(directory):
    """
    Removes all files (except .gitkeep) in the given directory and its subdirectories.
    """

    for filename in os.listdir(directory):  # iterate through each file in the directory
        file_path = os.path.join(directory, filename)  # create a complete filepath

        if os.path.isfile(file_path):  # if the path is a file
            os.remove(file_path)  # remove the file
        elif os.path.isdir(file_path):  # if the path is a directory
            remove_files(file_path)  # call this function recursively to remove files in subdirectory

def reset_updates(output_directory):
    """
    Resets the updates by deleting the contents of the updates.json file and removing any existing output files.
    """

    json_update_path = os.path.join(output_directory, "updates.json")  # path to the relevant update.json file
    output_path = output_directory  # path to the relevant output directory

    # Reset the json file - Writing an empty json object to the file effectively clears it
    with open(json_update_path, 'w') as f:
        json.dump({}, f)

    # Remove output files
    if os.path.exists(output_path):  # if the output directory exists
        remove_files(output_path)  # call the remove_files function to remove all files in the directory

def init_project(output_directory):
    """
    Initializes the profile with the given profile name.

    :param output_directory: The output directory to initialize.
    :return: None
    """

    # Path to the file where updates will be stored
    updates_file_path = os.path.join(output_directory, "updates.json")

    # Path to the .fraghub file
    fraghub_file_path = os.path.join(output_directory, ".fraghub")

    # List of main directories to be created under the output directory
    main_directories = ['CSV', 'JSON', 'MSP']

    # List of subdirectories to be created under each main directory
    sub_directories = ['NEG', 'POS']

    # Checking if the output directory for the profile already exists
    if not os.path.isdir(output_directory):
        # Creating the output directory since it doesn't exist
        os.makedirs(output_directory)

    # Checking if the updates file already exists
    if not os.path.isfile(updates_file_path):
        # Creating the updates file since it doesn't exist
        with open(updates_file_path, 'w') as fp:
            # Write an empty JSON object to the file as initial content
            json.dump({}, fp)

    # Create an empty .fraghub file if it doesn't already exist
    if not os.path.isfile(fraghub_file_path):
        with open(fraghub_file_path, 'w') as fp:
            # Leave the file empty (just create it)
            pass

    # For each directory in main_directories list...
    for main_dir in main_directories:
        # ...and for each subdirectory in sub_directories list...
        for sub_dir in sub_directories:
            # ...create a full path to the subdirectory
            dir_path = os.path.join(output_directory, main_dir, sub_dir)
            # If the directory path doesn't exist...
            if not os.path.isdir(dir_path):
                # ...create the subdirectory
                os.makedirs(dir_path)

    # Create the DELETED_SPECTRUMS directory
    deleted_spectrums_dir = os.path.join(output_directory, "DELETED_SPECTRUMS")
    if not os.path.isdir(deleted_spectrums_dir):
        os.makedirs(deleted_spectrums_dir)
