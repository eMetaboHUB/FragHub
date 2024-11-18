from tkinter import Tk, ttk, Checkbutton, Button, Label, Entry, StringVar, IntVar, LEFT, Frame, RIGHT, PhotoImage, Canvas, filedialog
import json
import os

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
    This function is triggered when the 'Done' button of the GUI is clicked.

    :return: None
    """

    convert_dict_values_to_float([
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
    ])

    root.destroy()

def convert_dict_values_to_float(keys):
    """
    Convert specified dictionary values to float.

    For each key in the provided list, this function checks if the key
    exists in the global parameters_dict dictionary. If the key is present,
    it converts the corresponding value to a float and updates the dictionary
    in place.

    Parameters:
        keys (list): A list of keys to convert values to float in the
        parameters_dict dictionary.
    """
    for key in keys:
        if key in parameters_dict:
            parameters_dict[key] = float(parameters_dict[key].get())

def build_window():
    """
    Builds the main application window for FragHub version 1.0.0. This function initializes the GUI including
    window settings, tab control, frames, buttons, and other widgets necessary for the application functionality.

     - Initializes the Tkinter root and sets its icon and title.
     - Creates a main frame and bottom frame for organizing widgets.
     - Adds several tabs to the window for different settings and options.
     - Configures widgets like buttons, labels, checkbuttons, and combobox for user interactions.
     - Loads a logo into the window and sets up directory selection widgets.
     - Sets up functionalities where user inputs directly update the application's parameters dictionary.

    Returns:
        None
    """
    global root
    root = Tk()

    # Initialize StringVars after creating root
    parameters_dict['input_directory'] = StringVar()
    parameters_dict['output_directory'] = StringVar()

    root.iconbitmap('../FragHub.ico')
    root.title("FragHub 1.2.0")

    # Create a main frame
    main_frame = Frame(root)
    main_frame.pack(fill='both', expand=True)

    # Create Tab Control
    tabControl = ttk.Notebook(main_frame)
    # Enter code here to add tabs...

    # Pack the tab control
    tabControl.pack(side="top", fill='both', expand=True)

    # Create the bottom frame
    bottom_frame = Frame(main_frame)
    bottom_frame.pack(side='top', fill='both', expand=True)

    # Create and place the "Done" button
    done_button = Button(bottom_frame, text="Done", command=on_done_button_clicked, height=3, width=10)
    done_button.place(relx=0.5, rely=0.5, anchor='c')

    # Create and place the logo label at the right of the window
    logo_frame = Frame(bottom_frame)
    logo_frame.pack(side='right', fill='y')

    # Create and position the logo
    logo = PhotoImage(file='../LOGO.png')  # Replace with the path to your image file
    logo = logo.subsample(6, 6)
    label_logo = Label(logo_frame, image=logo)
    label_logo.pack(side='right')

    # New tab
    tab4 = Frame(tabControl)

    # Existing tabs
    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tab3 = Frame(tabControl)
    tab3 = Frame(tabControl)
    tab4 = Frame(tabControl)
    tab_input = Frame(tabControl)
    tab_output = Frame(tabControl)
    # Existing tabs

    # Existing tabs
    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tab3 = Frame(tabControl)

    # Add the new tab to the tabControl
    tabControl.add(tab_input, text='INPUT')
    tabControl.add(tab_output, text='OUTPUT')
    tabControl.add(tab2, text='Filters Settings')
    tabControl.add(tab3, text='Output Settings')
    tabControl.add(tab4, text='Projects')
    tabControl.add(tab1, text='Update Settings')
    tabControl.add(tab_input, text='INPUT')
    add_directory_browsing(tab_input, 'Select input directory', 'input_directory')
    add_directory_browsing(tab_output, 'Select output directory', 'output_directory')
    tabControl.add(tab_output, text='OUTPUT')

    # Make a list of profiles
    profile_list = ['Profile 1', 'Profile 2', 'Profile 3']

    # Create a StringVar for the dropdown/listbox
    selected_profile = StringVar()

    # Create the Combobox, link it with the StringVar, and set the options
    profiles_dropdown = ttk.Combobox(tab4, textvariable=selected_profile)
    profiles_dropdown['values'] = [file.replace('.json', '') for file in os.listdir("../datas/updates") if
                                   file.endswith('.json')]
    selected_profile.set('basic')
    parameters_dict["selected_profile"] = "basic"
    profiles_dropdown.pack()

    # Update parameters_dict when a new selection is made
    def update_parameters_dict(*args):
        """
        Builds and displays the main application window.

        Returns:
            None
        """
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
    Checkbutton(tab1, text='Reset updates', variable=parameters_dict['reset_updates'], onvalue=True, offvalue=False,
                fg='red', font=("Helvetica", 9, 'bold')).pack()

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

    root.mainloop()


def add_directory_browsing(tab, label_text, key):
    """
    Create a frame for directory browsing within a specified tab and configure its elements.
    This includes a label, an entry widget for displaying the selected directory, and a browse button.

    Arguments:
        tab (Frame): The parent frame or tab where the directory browsing components will be created.
        label_text (str): The text for the label describing the directory selection input.
        key (str): The key to access the directory variable in the parameters dictionary.

    Parameters:
        tab (Frame): The parent frame or tab where the directory browsing components will be created.
        label_text (str): The text for the label describing the directory selection input.
        key (str): The key to access the directory variable in the parameters dictionary.

    """
    # Create a frame for directory browsing in the tab and pack it
    directory_frame = Frame(tab)
    directory_frame.pack(side='top', fill='x', pady=10)

    # Add a label to describe the directory selection
    directory_label = Label(directory_frame, text=label_text)
    directory_label.pack(side='top', padx=5)

    # Entry to show the selected directory
    directory_var = parameters_dict[key]
    directory_entry = Entry(directory_frame, textvariable=directory_var, width=40)
    directory_entry.pack(side='top', padx=5)

    # Button to browse for directory
    browse_button = Button(directory_frame, text="Browse", command=lambda: browse_directory(directory_var, key))
    browse_button.pack(side='top', padx=5)


def browse_directory(directory_var, key):
    """
    Sets the specified directory to the directory_var and updates parameters_dict with the selected directory's absolute path.

    Parameters:
        directory_var (tkinter.StringVar): The variable that stores the directory path.
        key (str): The key used to update parameters_dict.
    """
    directory = filedialog.askdirectory()
    if directory:
        directory_var.set(directory)
        # Update parameters_dict with the absolute path
        parameters_dict[key] = directory

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
    ouput_path = rf"../OUTPUT/{profile_name}"  # path to the relevant output directory

    # Reset the json file - Writing an empty json object to the file effectively clears it
    with open(json_update_path, 'w') as f:
        json.dump({}, f)

    # Remove output files
    if os.path.exists(ouput_path):  # if the output directory exists
        remove_files(ouput_path)  # call the remove_files function to remove all files in the directory

def init_profile(profile_name):
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
