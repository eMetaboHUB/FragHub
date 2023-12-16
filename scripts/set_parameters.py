from tkinter import Tk, Checkbutton, Button, Label, Entry, StringVar, IntVar, LEFT, Frame, RIGHT

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
    for func in func_names:
        parameters_dict[func] = float(parameters_dict[func].get())

    for key in [
        'check_minimum_peak_requiered_n_peaks',
        'reduce_peak_list_max_peaks',
        'keep_mz_in_range_from_mz',
        'keep_mz_in_range_to_mz',
        'check_minimum_of_high_peaks_requiered_intensity_percent',
        'check_minimum_of_high_peaks_requiered_no_peaks'
    ]:
        if key in parameters_dict:
            parameters_dict[key] = float(parameters_dict[key].get())

    parameters_dict['check_minimum_peak_requiered'] = 1.0

    root.destroy()


def build_window():
    global root
    root = Tk()
    root.title("Set filters parameters")

    for func in func_names:
        parameters_dict[func] = StringVar()
        parameters_dict[func].set(True)

        frame = Frame(root)
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