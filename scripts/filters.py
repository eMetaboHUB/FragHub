import pandas as pd
import re


def remove_peaks_above_precursor_mz(spectrum):
    if spectrum != None:
        Precursor_MZ = None
        if re.search(r"PRECURSORMZ: (.*)",spectrum):
            try:
                Precursor_MZ = float(re.search(r"PRECURSORMZ: (.*)",spectrum).group(1)) + 5.0
            except:
                return None

        peak_list = None
        if re.search(r"(NUM PEAKS: \d*\n)(((\d+\.?\d*)\s+(\d+\.?\d*)\n)*)",spectrum):
            peak_list = re.search(r"(NUM PEAKS: \d*\n)(((\d+\.?\d*)\s+(\d+\.?\d*)\n)*)",spectrum).group(2)

        if peak_list != None and Precursor_MZ != None:

            peaks_dict = {"MZ":[], "intensity":[]}

            peaks = re.findall(r"((\d+\.?\d*)\s+(\d+\.?\d*))",peak_list)

            for matchs in peaks:
                peaks_dict["MZ"].append(matchs[1])
                peaks_dict["intensity"].append(matchs[2])

            peaks_df = pd.DataFrame(peaks_dict)
            peaks_df['MZ'] = peaks_df['MZ'].astype(float)
            peaks_df['intensity'] = peaks_df['intensity'].astype(float)

            peaks_df = peaks_df[peaks_df['MZ'] < Precursor_MZ]

            if not peaks_df.empty:
                # Réécriture de la liste de pics
                peak_list_string = peaks_df.to_string(header=False, index=False)
                peak_list_string = "\n".join(["\t".join(line.split()) for line in peak_list_string.split("\n")])

                spectrum = re.sub(r"(NUM PEAKS: \d*\n)(((\d+\.?\d*)\s+(\d+\.?\d*)\n)*)", f"NUM PEAKS: {len(peaks_df)}\n{peak_list_string}\n", spectrum)

                return spectrum
            else:
                return None
    else:
        return None




