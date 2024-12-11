from GUI.utils.global_vars import parameters_dict
import json
import os

def remove_files(directory):
    """
    Removes all files (except .gitkeep) in the given directory and its subdirectories.
    """

    for filename in os.listdir(directory):  # iterate through each file in the directory
        if filename == '.gitkeep':
            continue  # skip .gitkeep files as we don't want to remove these
        file_path = os.path.join(directory, filename)  # create a complete filepath

        if os.path.isfile(file_path):  # if the path is a file
            os.remove(file_path)  # remove the file
        elif os.path.isdir(file_path):  # if the path is a directory
            remove_files(file_path)  # call this function recursively to remove files in subdirectory

def reset_updates(profile_name):
    """
    Resets the updates by deleting the contents of the updates.json file and removing any existing output files.
    """

    json_update_path = rf"../datas/updates/{profile_name}.json"  # path to the relevant update.json file
    ouput_path = os.path.join(parameters_dict["output_directory"],profile_name)  # path to the relevant output directory

    # Reset the json file - Writing an empty json object to the file effectively clears it
    with open(json_update_path, 'w') as f:
        json.dump({}, f)

    # Remove output files
    if os.path.exists(ouput_path):  # if the output directory exists
        remove_files(ouput_path)  # call the remove_files function to remove all files in the directory

def init_project(profile_name):
    """
    Initializes the profile with the given profile name.

    :param profile_name: The name of the profile to initialize.
    :return: None
    """

    # Path to the file where updates will be stored
    updates_file_path = os.path.join("../datas/updates", profile_name + ".json")

    # Directory where outputs will be stored for the profile
    output_directory = os.path.join(parameters_dict["output_directory"],profile_name)

    # List of main directories to be created under the output directory
    main_directories = ['CSV', 'JSON', 'MSP']

    # List of subdirectories to be created under each main directory
    sub_directories = ['NEG', 'POS']

    # Checking if the updates file already exists
    if not os.path.isfile(updates_file_path):
        # Creating the updates file since it doesn't exist
        with open(updates_file_path, 'w') as fp:
            # Write an empty JSON object to the file as initial content
            json.dump({}, fp)

    # Checking if the output directory for the profile already exists
    if not os.path.isdir(output_directory):
        # Creating the output directory since it doesn't exist
        os.makedirs(output_directory)

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
