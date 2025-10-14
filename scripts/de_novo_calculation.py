from scripts.GUI.utils.global_vars import parameters_dict
from scripts.calculate_maximized_chunk_size import *
from scripts.globals_vars import atoms_of_life
import concurrent.futures
from math import floor
from numba import jit
import pandas as pd
import numpy as np
import re

ppm_tol = parameters_dict.get('de_novo_ppm_tolerance', 5)


@jit(nopython=True, nogil=True)
def _generate_combinations_numba(element_index, current_mass, composition_array,
                                 element_masses, max_counts,
                                 min_mass, max_mass, max_remaining_masses):
    """
    Numba-compatible recursive function for combination generation.
    Returns a list of NumPy arrays, each array being a valid composition.
    """
    if current_mass + max_remaining_masses[element_index] < min_mass:
        # The return type must be consistent, so we return an empty list
        return [np.empty((0,), dtype=np.int64)]

    # Base case: we have considered all elements
    if element_index == len(element_masses):
        if current_mass >= min_mass:
            # Return the found composition in a list
            return [composition_array.copy()]
        else:
            return [np.empty((0,), dtype=np.int64)]

    # Initialize a list to store the results of this branch
    results = []

    element_mass = element_masses[element_index]

    for count in range(max_counts[element_index] + 1):
        new_mass = current_mass + count * element_mass
        if new_mass > max_mass:
            break

        composition_array[element_index] = count

        # Recursive call
        sub_results = _generate_combinations_numba(element_index + 1, new_mass, composition_array,
                                                   element_masses, max_counts,
                                                   min_mass, max_mass, max_remaining_masses)

        # Add the valid compositions found in the sub-branch
        for res in sub_results:
            if res.size > 0: # Check if the list is not the empty marker
                results.append(res)

    # Reset the composition for backtracking
    composition_array[element_index] = 0
    return results

# --- FUNCTION #2: THE SIMPLE CASE (CHNOPS) ---
def annotate_simple_chnops_case(spectrum_row, atoms, ppm_tolerance):
    """
    Ultra-optimized version with Numba for combinatorial search.
    """
    max_composition_global, peaks, proton_mass, elements_to_test = _prepare_annotation_data(spectrum_row, atoms)

    if not peaks:
        return {}

    # Sort elements (unchanged, always crucial)
    elements_to_test.sort(key=lambda el: atoms[el], reverse=True)

    # Create a mapping table for the index of 'C', 'H', 'N', 'P'
    # for the plausibility function.
    c_idx, h_idx, n_idx, p_idx = -1, -1, -1, -1
    for i, el in enumerate(elements_to_test):
        if el == 'C': c_idx = i
        elif el == 'H': h_idx = i
        elif el == 'N': n_idx = i
        elif el == 'P': p_idx = i

    # Convert data to NumPy arrays for Numba
    element_masses_np = np.array([atoms[el] for el in elements_to_test], dtype=np.float64)

    # The plausibility function is applied in Python to the Numba results
    def _is_plausible_formula(comp_array):
        c = comp_array[c_idx] if c_idx != -1 else 0
        if c == 0: return False, None # At least one carbon is required

        h = comp_array[h_idx] if h_idx != -1 else 0
        n = comp_array[n_idx] if n_idx != -1 else 0
        p = comp_array[p_idx] if p_idx != -1 else 0

        dbe = c - (h / 2) + (n / 2) + (p / 2) + 1
        if dbe < 0 or dbe != floor(dbe):
            return False, None
        if not (0.2 <= h / c <= 3.0):
            return False, None
        return True, int(dbe)

    annotation_results = {}
    for mz, intensity in peaks:
        target_neutral_mass = mz - proton_mass
        min_mass = target_neutral_mass * (1 - ppm_tolerance / 1_000_000)
        max_mass = target_neutral_mass * (1 + ppm_tolerance / 1_000_000)

        # Create composition limits for this peak (in NumPy)
        max_counts_list = []
        for el in elements_to_test:
            max_count = int(max_mass // atoms[el])
            max_counts_list.append(min(max_composition_global.get(el, 0), max_count))
        max_counts_np = np.array(max_counts_list, dtype=np.int64)

        # Pre-calculate remaining masses (in NumPy)
        max_remaining_masses_np = np.zeros(len(elements_to_test) + 1, dtype=np.float64)
        for i in range(len(elements_to_test) - 1, -1, -1):
            el_max_mass = max_counts_np[i] * element_masses_np[i]
            max_remaining_masses_np[i] = max_remaining_masses_np[i + 1] + el_max_mass

        # Start the search with Numba
        initial_composition = np.zeros(len(elements_to_test), dtype=np.int64)
        valid_compositions_raw = _generate_combinations_numba(
            0, 0.0, initial_composition,
            element_masses_np, max_counts_np,
            min_mass, max_mass, max_remaining_masses_np
        )

        # Process and format results in Python
        results_for_peak = []
        for comp_array in valid_compositions_raw:
            if comp_array.size == 0: continue

            is_plausible, dbe = _is_plausible_formula(comp_array)
            if is_plausible:
                # Reconstruct the dictionary and calculate the final mass
                current_mass = np.sum(comp_array * element_masses_np)
                ion_composition = {elements_to_test[i]: comp_array[i] for i in range(len(elements_to_test))}
                ion_composition['H'] += 1

                formula = generate_hill_formula_string(ion_composition) + '+'
                error = ((current_mass - target_neutral_mass) / target_neutral_mass) * 1_000_000
                results_for_peak.append({
                    'formula': formula,
                    'calculated_mass': round(current_mass, 6),
                    'error_ppm': round(error, 2),
                    'dbe': dbe
                })

        if results_for_peak:
            annotation_results[mz] = results_for_peak

    return annotation_results


# --- FUNCTION #3: THE COMPLEX CASE (ALL OTHER ATOMS) ---
def annotate_complex_case(spectrum_row, atoms, ppm_tolerance):
    """
    Version for complex cases, optimized with Numba and algorithmic pruning.
    """
    valences = {'C': 4, 'N': 3, 'P': 3, 'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'Na': 1, 'K': 1}
    max_composition_global, peaks, proton_mass, elements_to_test = _prepare_annotation_data(spectrum_row, atoms)

    if not peaks:
        return {}

    # --- Sort elements by mass (always crucial) ---
    elements_to_test.sort(key=lambda el: atoms[el], reverse=True)

    # --- Prepare data for Numba and plausibility ---
    element_masses_np = np.array([atoms[el] for el in elements_to_test], dtype=np.float64)

    # Create indexes for the plausibility function
    el_indices = {el: i for i, el in enumerate(elements_to_test)}
    c_idx = el_indices.get('C', -1)
    h_idx = el_indices.get('H', -1)
    o_idx = el_indices.get('O', -1)
    n_idx = el_indices.get('N', -1)
    p_idx = el_indices.get('P', -1)
    s_idx = el_indices.get('S', -1)

    def _is_plausible_complex_formula_np(comp_array, calculated_mass):
        # Nitrogen rule
        safe_elements_for_nitrogen_rule = {'C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}
        contains_special_atoms = any(
            elements_to_test[i] not in safe_elements_for_nitrogen_rule
            for i, count in enumerate(comp_array) if count > 0
        )
        if not contains_special_atoms:
            nitrogen_type_atom_count = (comp_array[n_idx] if n_idx != -1 else 0) + \
                                       (comp_array[p_idx] if p_idx != -1 else 0)
            if (round(calculated_mass) % 2) != (nitrogen_type_atom_count % 2):
                return False, None

        # DBE rule
        dbe_contribution = sum(
            count * (valences.get(elements_to_test[i], 2) - 2)
            for i, count in enumerate(comp_array) if count > 0
        )
        dbe = 1 + dbe_contribution / 2
        if dbe < 0 or dbe != floor(dbe):
            return False, None

        # Senior ratios
        c = comp_array[c_idx] if c_idx != -1 else 0
        if c > 0:
            h = comp_array[h_idx] if h_idx != -1 else 0
            if not (0.2 <= h / c <= 3.0): return False, None
            if o_idx != -1 and (comp_array[o_idx] / c > 1.2): return False, None
            if n_idx != -1 and (comp_array[n_idx] / c > 1.3): return False, None
            if p_idx != -1 and (comp_array[p_idx] / c > 0.3): return False, None
            if s_idx != -1 and (comp_array[s_idx] / c > 0.8): return False, None

        return True, int(dbe)

    annotation_results = {}
    for mz, intensity in peaks:
        target_neutral_mass = mz - proton_mass
        min_mass = target_neutral_mass * (1 - ppm_tolerance / 1_000_000)
        max_mass = target_neutral_mass * (1 + ppm_tolerance / 1_000_000)

        # Composition limits per peak (in NumPy)
        max_counts_list = [min(max_composition_global.get(el, 0), int(max_mass // atoms[el])) for el in elements_to_test]
        max_counts_np = np.array(max_counts_list, dtype=np.int64)

        # Pre-calculate remaining masses (in NumPy)
        max_remaining_masses_np = np.zeros(len(elements_to_test) + 1, dtype=np.float64)
        for i in range(len(elements_to_test) - 1, -1, -1):
            max_remaining_masses_np[i] = max_remaining_masses_np[i + 1] + (max_counts_np[i] * element_masses_np[i])

        # Start the search with Numba
        initial_composition = np.zeros(len(elements_to_test), dtype=np.int64)
        valid_compositions_raw = _generate_combinations_numba(
            0, 0.0, initial_composition, element_masses_np, max_counts_np,
            min_mass, max_mass, max_remaining_masses_np
        )

        # Process and format results in Python
        results_for_peak = []
        for comp_array in valid_compositions_raw:
            if comp_array.size == 0 or (c_idx != -1 and comp_array[c_idx] == 0):
                continue

            current_mass = np.sum(comp_array * element_masses_np)
            is_plausible, dbe = _is_plausible_complex_formula_np(comp_array, current_mass)

            if is_plausible:
                ion_composition_dict = {elements_to_test[i]: comp_array[i] for i in range(len(elements_to_test))}
                ion_composition_dict['H'] = (ion_composition_dict.get('H', 0)) + 1

                formula = generate_hill_formula_string(ion_composition_dict) + '+'
                error = ((current_mass - target_neutral_mass) / target_neutral_mass) * 1_000_000
                results_for_peak.append({
                    'formula': formula,
                    'calculated_mass': round(current_mass, 6),
                    'error_ppm': round(error, 2),
                    'dbe': dbe
                })

        if results_for_peak:
            annotation_results[mz] = results_for_peak

    return annotation_results


def generate_hill_formula_string(composition):
    """
    Generates a Hill system formatted chemical formula string from a composition dict.
    Example: {'C': 6, 'H': 12, 'O': 6} -> "C6H12O6"
    """
    if not composition:
        return ""

    formula_str = ""
    # Create a copy without zero-count elements
    comp = {k: v for k, v in composition.items() if v > 0}

    # Handle C and H first as per Hill system convention
    if 'C' in comp:
        count = comp.pop('C')
        formula_str += f"C{count if count > 1 else ''}"
    if 'H' in comp:
        count = comp.pop('H')
        formula_str += f"H{count if count > 1 else ''}"

    # Append the rest of the elements sorted alphabetically
    for element in sorted(comp.keys()):
        count = comp[element]
        formula_str += f"{element}{count if count > 1 else ''}"

    return formula_str


def _prepare_annotation_data(spectrum_row, atoms):
    """
    Helper function: parses a spectrum row (pandas.Series) and returns essential data.
    """
    precursor_formula = spectrum_row.get("FORMULA", "")
    max_composition = {}

    # Use regex to find all element-count pairs in the formula string
    for element, count_str in re.findall(r'([A-Z][a-z]*)(\d*)', precursor_formula):
        if element in atoms:
            max_composition[element] = int(count_str) if count_str else 1

    peaks_list_str = spectrum_row.get("PEAKS_LIST", "")
    if not peaks_list_str:
        return None, None, None, None

    # Handle both semicolon-separated and newline-separated peak formats
    if ';' in peaks_list_str:
        peaks = [(float(mz), float(i)) for p in peaks_list_str.strip().split(';') for mz, i in [p.split()]]
    else:
        peaks = [(float(mz), float(i)) for l in peaks_list_str.strip().split('\n') for mz, i in [l.split()]]

    proton_mass = atoms.get('H')  # Safer access
    elements_to_test = list(max_composition.keys())

    return max_composition, peaks, proton_mass, elements_to_test


# --- MAIN DISPATCHER FUNCTION ---
def parse_and_annotate_spectrum(spectrum_row, atoms, ppm_tolerance):
    """
    Parses the formula from a spectrum row, determines if it's a simple (CHNOPS)
    or complex case, and calls the appropriate annotation function.
    'spectrum_row' is expected to be a pandas Series (a row from a DataFrame).
    """
    precursor_formula = spectrum_row.get("FORMULA", "")
    if not precursor_formula:
        raise ValueError("The 'FORMULA' column is missing or empty in the spectrum row.")

    present_elements = set(el for el, count in re.findall(r'([A-Z][a-z]*)(\d*)', precursor_formula))
    simple_atoms = {'C', 'H', 'N', 'O', 'P', 'S'}

    if present_elements.issubset(simple_atoms):
        return annotate_simple_chnops_case(spectrum_row, atoms, ppm_tolerance)
    else:
        return annotate_complex_case(spectrum_row, atoms, ppm_tolerance)


def process_single_spectrum(spectrum_dict):
    """
    Processes a single spectrum (represented by a dictionary).
    This is the function that will be executed by each thread.
    """
    # Step 1: Calculate annotations
    # Note: ensure that 'parse_and_annotate_spectrum' and 'atoms_of_life'
    # are accessible in this scope.
    annotations = parse_and_annotate_spectrum(spectrum_dict, atoms_of_life, ppm_tolerance=ppm_tol)

    # Step 2: Generate the new PEAKS_LIST string
    original_peak_list = spectrum_dict.get('PEAKS_LIST', "")

    def _generate_updated_peak_list_string(peak_list_str, annots):
        if not annots or not peak_list_str or pd.isna(peak_list_str):
            return peak_list_str
        lines = str(peak_list_str).strip().split('\n')
        keys = list(annots.keys())
        new_lines = []
        for line in lines:
            parts = line.split()
            if len(parts) < 2: continue
            mz_str, intensity_str = parts[0], parts[1]
            mz_val = float(mz_str)
            closest_key = min(keys, key=lambda k: abs(k - mz_val))
            formula_str = ""
            if abs(closest_key - mz_val) < 1e-5:
                data = annots.get(closest_key, [{}])[0]
                if 'formula' in data:
                    formula_str = f"{data['formula']}/{data['error_ppm']}"
            new_lines.append(f"{mz_str} {intensity_str} {formula_str}")
        return "\n".join(new_lines)

    updated_list = _generate_updated_peak_list_string(original_peak_list, annotations)

    # Update the dictionary with the new results
    spectrum_dict['annotation_results'] = annotations
    spectrum_dict['PEAKS_LIST'] = updated_list

    return spectrum_dict

def de_novo_calculation(spectrum_df, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Calculates de novo annotations using a ThreadPoolExecutor and batch processing.
    """
    # Configure callbacks
    if prefix_callback:
        prefix_callback("Calculating de novo formulas")
    if item_type_callback:
        item_type_callback("spectra")

    # 1. Convert the DataFrame into a list of dictionaries for processing
    records = spectrum_df.to_dict('records')
    total_items = len(records)

    if total_items_callback:
        total_items_callback(total_items, 0)

    # 2. Calculate the batch size (chunks)
    chunk_size = calculate_maximized_chunk_size(data_list=records)

    # List to collect the final results
    final_results = []
    processed_items = 0

    # 3. Loop over the list, processing one batch at a time
    for i in range(0, total_items, chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Extract the current batch
            chunk = records[i:i + chunk_size]

            # Execute the processing function on each dictionary in the batch
            results = list(executor.map(process_single_spectrum, chunk))

        # Add the results of the processed batch to the final list
        final_results.extend(results)

        # Update the number of processed items
        processed_items += len(chunk)

        # 4. Update the progress bar after each batch
        if progress_callback:
            progress_callback(processed_items)

    # 5. Convert the list of dictionaries back into a DataFrame
    final_df = pd.DataFrame(final_results)

    # Remove the 'annotation_results' column before returning the result.
    # 'errors='ignore'' prevents an error if the column is not found.
    return final_df.drop(columns=['annotation_results'], errors='ignore')

