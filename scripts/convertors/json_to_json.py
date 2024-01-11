import concurrent.futures
from tqdm import tqdm
import re

global peak_list_json_to_json_pattern
peak_list_json_to_json_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

def flatten_dict(d, parent_key='', flat_dict=None):
    if flat_dict is None:
        flat_dict = {}
    for k, v in d.items():
        new_key = f"{parent_key}.{k}" if parent_key else k
        if isinstance(v, dict):
            flatten_dict(v, new_key, flat_dict)
        else:
            flat_dict[new_key] = v

    return flat_dict

def json_to_json(json_dict):
    json_dict = flatten_dict(json_dict)
    if "spectrum" in json_dict:
        spectrum_string = json_dict["spectrum"]

        # Find all matches and convert them to floats
        spectrum_pairs = [[float(g) for g in match] for match in re.findall(peak_list_json_to_json_pattern, spectrum_string)]

        json_dict["spectrum"] = spectrum_pairs

    return json_dict

def json_to_json_processing(FINAL_JSON):
    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(FINAL_JSON), unit=" spectrums", colour="green", desc="{:>80}".format("converting MSP spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(FINAL_JSON), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = FINAL_JSON[i:i + chunk_size]
            results = list(executor.map(json_to_json, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final