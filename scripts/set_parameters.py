from tkinter import Tk, Checkbutton, Button, Label, Entry, StringVar, IntVar, LEFT, Frame, RIGHT
from tkinter import ttk  # import ttk module
import json
import os

global parameters_dict
parameters_dict = {}

func_names = [
    'check_minimum_peak_requiered',
    'remove_peak_above_precursormz',
    'reduce_peak_list',
    'normalize_intensity',
    'keep_mz_in_range',
    'check_minimum_of_high_peaks_requiered'
]

def on_done_button_clicked():
    """
    Closes the root window and converts the values in the parameters_dict dictionary to floats.

    :return: None
    """
    for key in [
        'reset_updates',
        'csv',
        'msp',
        'json'
    ] + func_names + [
        'check_minimum_peak_requiered_n_peaks',
        'reduce_peak_list_max_peaks',
        'keep_mz_in_range_from_mz',
        'keep_mz_in_range_to_mz',
        'check_minimum_of_high_peaks_requiered_intensity_percent',
        'check_minimum_of_high_peaks_requiered_no_peaks'
    ]:
        if key in parameters_dict:
            parameters_dict[key] = float(parameters_dict[key].get())

    root.destroy()

def build_window():
    """
    Builds a window for setting filters parameters.

    :return: None
    """
    global root
    root = Tk()
    root.title("FragHub")

    # Create Tab Control
    tabControl = ttk.Notebook(root)

    # New tab
    tab4 = Frame(tabControl)

    # Existing tabs
    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tab3 = Frame(tabControl)

    # Add the new tab to the tabControl
    tabControl.add(tab2, text='Filters Settings')
    tabControl.add(tab3, text='Output Settings')
    tabControl.add(tab4, text='Profils')
    tabControl.add(tab1, text='Update Settings')

    tabControl.pack(expand=1, fill='both')

    # Make a list of profiles
    profile_list = ['Profile 1', 'Profile 2', 'Profile 3']

    # Create a StringVar for the dropdown/listbox
    selected_profile = StringVar()

    # Create the Combobox, link it with the StringVar, and set the options
    profiles_dropdown = ttk.Combobox(tab4, textvariable=selected_profile)
    profiles_dropdown['values'] = [file.replace('.json','') for file in os.listdir("../datas/updates") if file.endswith('.json')]
    selected_profile.set('basic')
    parameters_dict["selected_profile"] = "basic"
    profiles_dropdown.pack()

    # Update parameters_dict when a new selection is made
    def update_parameters_dict(*args):
        parameters_dict['selected_profile'] = selected_profile.get()

    selected_profile.trace_add("write", update_parameters_dict)

    # Create Checkbuttons in the third tab
    parameters_dict['csv'] = StringVar()
    parameters_dict['csv'].set("1.0")
    Checkbutton(tab3, text=".csv", variable=parameters_dict['csv'], onvalue="1.0", offvalue="0.0").pack()

    parameters_dict['msp'] = StringVar()
    parameters_dict['msp'].set("1.0")
    Checkbutton(tab3, text=".msp", variable=parameters_dict['msp'], onvalue="1.0", offvalue="0.0").pack()

    parameters_dict['json'] = StringVar()
    parameters_dict['json'].set("1.0")
    Checkbutton(tab3, text=".json", variable=parameters_dict['json'], onvalue="1.0", offvalue="0.0").pack()

    # Add Checkbox to tab1
    parameters_dict['reset_updates'] = IntVar()
    Checkbutton(tab1, text='Reset updates', variable=parameters_dict['reset_updates'], onvalue=True, offvalue=False, fg='red', font=("Helvetica", 9, 'bold')).pack()

    for func in func_names:
        parameters_dict[func] = StringVar()
        parameters_dict[func].set(True)

        frame = Frame(tab2)
        frame.pack(fill='x')

        frame_func = Frame(frame)
        frame_func.pack(side=LEFT)

        Checkbutton(
            frame_func,
            variable=parameters_dict[func],
            onvalue=True,
            offvalue=False
        ).pack(side=LEFT)

        Label(frame_func, text=func).pack(side=LEFT)

        frame_params = Frame(frame)
        frame_params.pack(side=RIGHT)

        if func == 'check_minimum_peak_requiered':
            Label(frame_params, text="n_peaks:").pack(side=LEFT)
            parameters_dict['check_minimum_peak_requiered_n_peaks'] = IntVar()
            parameters_dict['check_minimum_peak_requiered_n_peaks'].set(3)
            Entry(frame_params, textvariable=parameters_dict['check_minimum_peak_requiered_n_peaks']).pack(side=LEFT)
        elif func == 'reduce_peak_list':
            Label(frame_params, text="max_peaks:").pack(side=LEFT)
            parameters_dict['reduce_peak_list_max_peaks'] = IntVar()
            parameters_dict['reduce_peak_list_max_peaks'].set(500)
            Entry(frame_params, textvariable=parameters_dict['reduce_peak_list_max_peaks']).pack(side=LEFT)
        elif func == 'keep_mz_in_range':
            Label(frame_params, text="from_mz:").pack(side=LEFT)
            parameters_dict['keep_mz_in_range_from_mz'] = IntVar()
            parameters_dict['keep_mz_in_range_from_mz'].set(50)
            Entry(frame_params, textvariable=parameters_dict['keep_mz_in_range_from_mz']).pack(side=LEFT)
            Label(frame_params, text="to_mz:").pack(side=LEFT)
            parameters_dict['keep_mz_in_range_to_mz'] = IntVar()
            parameters_dict['keep_mz_in_range_to_mz'].set(2000)
            Entry(frame_params, textvariable=parameters_dict['keep_mz_in_range_to_mz']).pack(side=LEFT)
        elif func == 'check_minimum_of_high_peaks_requiered':
            Label(frame_params, text="intensity_percent:").pack(side=LEFT)
            parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent'] = IntVar()
            parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent'].set(5)
            Entry(frame_params,
                  textvariable=parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent']).pack(
                side=LEFT)
            Label(frame_params, text="no_peaks:").pack(side=LEFT)
            parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks'] = IntVar()
            parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks'].set(2)
            Entry(frame_params, textvariable=parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks']).pack(
                side=LEFT)

    Button(root, text="Done", command=on_done_button_clicked).pack()

    root.mainloop()

def remove_files(directory):
    """
    Removes all files (except .gitkeep) in the given directory and its subdirectories.
    """
    for filename in os.listdir(directory):
        if filename == '.gitkeep':
            continue  # skip .gitkeep files

        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)  # remove the file
        elif os.path.isdir(file_path):
            remove_files(file_path)  # call this function again with the subdirectory

def reset_updates(profile_name):
    """
    Resets the updates by deleting the contents of the updates.json file and removing any existing output files.
    """
    json_update_path = rf"../datas/updates/{profile_name}.json"
    ouput_path = rf"../OUTPUT/{profile_name}"

    # Reset the json file - Writing empty json object
    with open(json_update_path, 'w') as f:
        json.dump({}, f)

    # Remove output files
    if os.path.exists(ouput_path):
        remove_files(ouput_path)

def init_profile(profile_name):
    """
    Initializes the profile with the given profile name.

    :param profile_name: The name of the profile to initialize.
    :return: None
    """
    # Get the file path
    updates_file_path = os.path.join("../datas/updates", profile_name + ".json")
    output_directory = os.path.join("../OUTPUT", profile_name)

    main_directories = ['CSV', 'JSON', 'MSP']
    sub_directories = ['NEG', 'POS']

    # Check if file doesn't exist
    if not os.path.isfile(updates_file_path):
        # Create the file
        with open(updates_file_path, 'w') as fp:
            # Write an empty JSON object to the file
            json.dump({}, fp)

    if not os.path.isdir(output_directory):
        # Create the directory
        os.makedirs(output_directory)

    # Loop through each main directory
    for main_dir in main_directories:
        # Loop through each sub directory
        for sub_dir in sub_directories:
            # Get the sub directory path
            dir_path = os.path.join(output_directory, main_dir, sub_dir)
            # If the directory doesn't exist, create it
            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
