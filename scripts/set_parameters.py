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
    for func in func_names:
        parameters_dict[func] = float(parameters_dict[func].get())

    for key in [
        'reset_updates',
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
    root.title("Set filters parameters")

    # Create Tab Control
    tabControl = ttk.Notebook(root)
    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tabControl.add(tab1, text='Update Settings')
    tabControl.add(tab2, text='Filters Settings')
    tabControl.pack(expand=1, fill='both')

    # Add Checkbox to tab1
    parameters_dict['reset_updates'] = IntVar()
    Checkbutton(tab1, text='Reset updates', variable=parameters_dict['reset_updates'], onvalue=True, offvalue=False).pack()

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
    Removes all files in the given directory and its subdirectories.
    """
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)  # remove the file
        elif os.path.isdir(file_path):
            remove_files(file_path)  # call this function again

def reset_updates():
    """
    Resets the updates by deleting the contents of the updates.json file and removing any existing output files.
    """
    json_update_path = r"../data/updates.json"
    ouput_path = r"../OUTPUT"

    # Reset the json file - Writing empty json object
    with open(json_update_path, 'w') as f:
        json.dump({}, f)

    # Remove output files
    if os.path.exists(ouput_path):
        remove_files(ouput_path)

