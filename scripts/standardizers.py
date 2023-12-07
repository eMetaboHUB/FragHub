
import re

def remove_no_inchikey(spectrum):
    """
    Remove No Inchikey

    This method removes the spectra that do not have an Inchikey value.

    :param spectrum: The input spectrum
    :return: The spectrum without the spectra with no Inchikey value

    """
    if spectrum != None:
        if re.search("INCHIKEY: \n|INCHIKEY: None\n",spectrum):
            return None
        else:
            return spectrum

def harmonize_adduct(spectrum):
    """
    :param spectrum: the spectrum containing the precursor information
    :return: the spectrum with the harmonized precursor adduct
    """
    if spectrum != None:
        if re.search("((^|\n)(PRECURSORTYPE:)) (.*)\n",spectrum):
            adduct = re.search("((^|\n)(PRECURSORTYPE:)) (.*)\n",spectrum).group(4)
            if adduct != "None":
                if "[" not in adduct or "]" not in adduct:
                    if re.search("((^|\n)(IONMODE:)) (.*)\n",spectrum):
                        ionmode = re.search("((^|\n)(IONMODE:)) (.*)\n",spectrum).group(4)
                        if ionmode != "None":
                            if ionmode.lower().startswith("p"):
                                adduct = "["+adduct+"]"+"+\n"
                                spectrum = re.sub("(^|\n)(PRECURSORTYPE:) (.*)\n",f"\nPRECURSORTYPE: {adduct}",spectrum)
                                return spectrum
                            elif ionmode.lower().startswith("n"):
                                adduct = "[" + adduct + "]" + "-\n"
                                spectrum = re.sub("((^|\n)(PRECURSORTYPE:)) (.*)\n", f"\nPRECURSORTYPE: {adduct}", spectrum)
                                return spectrum
                        else:
                            return spectrum
                    else:
                        return spectrum
                elif re.search("((PRECURSORTYPE:)) ((.*\+)\*)\n",spectrum):
                    adduct = re.search("((PRECURSORTYPE:)) ((.*\+)\*)\n", spectrum).group(4)
                    spectrum = re.sub("((PRECURSORTYPE:)) ((.*\+)\*)\n",f"PRECURSORTYPE: {adduct}\n",spectrum)
                    return spectrum
                elif re.search("((PRECURSORTYPE:)) ((.*\-)\*)\n",spectrum):
                    adduct = re.search("((PRECURSORTYPE:)) ((.*\-)\*)\n", spectrum).group(4)
                    spectrum = re.sub("((PRECURSORTYPE:)) ((.*\-)\*)\n",f"PRECURSORTYPE: {adduct}\n",spectrum)
                    return spectrum
                else:
                    return spectrum
            else:
                return spectrum
        else:
            return spectrum
    else:
        return spectrum

def remove_no_mass(spectrum):
    """
    Remove spectra with no mass information from the given spectrum.

    :param spectrum: The spectrum to remove spectra with no mass information from.
    :return: The spectrum with no mass spectra removed.
    """
    if spectrum != None:
        if re.search("PRECURSORMZ: None\n",spectrum) and re.search("PARENTMASS: None\n",spectrum):
            return None
        else:
            return spectrum

def correct_ionmode(spectrum):
    """
    Corrects the ion mode of the given spectrum.

    :param spectrum: The spectrum to correct the ion mode.
    :type spectrum: str
    :return: The corrected spectrum.
    :rtype: str
    """
    if spectrum != None:
        if re.search("\nIONMODE: n/a",spectrum):
            if re.search("PRECURSORTYPE: (.*)\-\n",spectrum):
                spectrum = re.sub("\nIONMODE: n/a","\nIONMODE: negative",spectrum)
                spectrum = re.sub("CHARGE: (.*)\n", "CHARGE: -1\n", spectrum)
                return spectrum
            elif re.search("PRECURSORTYPE: (.*)\+\n",spectrum):
                spectrum = re.sub("\nIONMODE: n/a","\nIONMODE: positive",spectrum)
                spectrum = re.sub("CHARGE: (.*)\n", "CHARGE: 1\n", spectrum)
                return spectrum
        else:
            return spectrum

def harmonize_retention_time(spectrum):
    """
    :param spectrum: A string representation of a spectrum.
    :return: A modified version of the spectrum with the retention time harmonized.

    The function `harmonize_retention_time` takes a spectrum as input and modifies it by harmonizing the retention time. If the input spectrum is not None, the function checks if the retention
    * time is set to 0. If so, it replaces the 0 with `None`. If the retention time is already present in the spectrum (in the format "RETENTIONTIME: value"), the function checks if the
    * value is a valid float. If it is, the spectrum is returned unchanged. If the value is not a valid float, it is replaced with `None` and the modified spectrum is returned.
    """
    if spectrum != None:
        if re.search("RETENTIONTIME: 0\n",spectrum):
            spectrum = re.sub("RETENTIONTIME: 0\n", "RETENTIONTIME: None\n", spectrum)
            return spectrum
        elif re.search("(^|\n)(RETENTIONTIME:) (.*)\n",spectrum):
            RT = re.search("((^|\n)(RETENTIONTIME:)) (.*)\n", spectrum).group(4)
            try:
                RT_test = float(RT)
                return spectrum
            except:
                spectrum = re.sub("((^|\n)(RETENTIONTIME:)) (.*)\n", "\nRETENTIONTIME: None\n", spectrum)
                return spectrum

def harmonize_ms_level(spectrum):
    """
    :param spectrum: The original string representation of the spectrum with MSLEVEL information.
    :return: The harmonized string representation of the spectrum with updated MSLEVEL information.

    This method takes a spectrum string and performs a series of substitutions on the MSLEVEL information in the string.
    The substitutions change the MSLEVEL values to a standard format, ensuring consistency in the representation.

    Here is an overview of the substitutions made:
    - "MSLEVEL: MS" is replaced with "MSLEVEL: 1"
    - "MSLEVEL: MS1" is replaced with "MSLEVEL: 1"
    - "MSLEVEL: MS2" is replaced with "MSLEVEL: 2"
    - "MSLEVEL: MS3" is replaced with "MSLEVEL: 3"
    - "MSLEVEL: MS4" is replaced with "MSLEVEL: 4"
    - "MSLEVEL: 2-MS4 Composite" is replaced with "MSLEVEL: 4"
    - "MSLEVEL: 2-MS5 Composite" is replaced with "MSLEVEL: 5"

    The original spectrum string is not modified, and the updated string with harmonized MSLEVEL information is returned.
    """
    if spectrum != None:
        spectrum = re.sub("MSLEVEL: MS\n", "MSLEVEL: 1\n", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: MS1","MSLEVEL: 1", spectrum,flags=re.I)
        spectrum = re.sub("MSLEVEL: MS2", "MSLEVEL: 2", spectrum,flags=re.I)
        spectrum = re.sub("MSLEVEL: MS3", "MSLEVEL: 3", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: MS4", "MSLEVEL: 4", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: 2-MS4 Composite", "MSLEVEL: 4", spectrum, flags=re.I)
        spectrum = re.sub("MSLEVEL: 2-MS5 Composite", "MSLEVEL: 5", spectrum, flags=re.I)

    return spectrum

# def harmonize_collisionenergy(spectrum):
#     """
#     :param spectrum: A string representing the input spectrum.
#     :return: The harmonized spectrum with updated collision energy.
#     """
#     if spectrum != None:
#         if re.search("(^|\n)(COLLISIONENERGY:) ([0-9]*)((?: )?)((?:\()?)([a-zA-Z]*)((?:\))?)\n", spectrum):
#             collisionenergy =  re.search("(^|\n)(COLLISIONENERGY:) ([0-9]*)((?: )?)((?:\()?)([a-zA-Z]*)((?:\))?)\n", spectrum)
#             if collisionenergy.group(3) == '':
#                 return spectrum
#             elif collisionenergy.group(3).isnumeric():
#                 if collisionenergy.group(6) == '':
#                     return spectrum
#                 elif collisionenergy.group(6).isalpha():
#                     fragmentation = collisionenergy.group(6)
#                     # deleting alpha caracters
#                     spectrum = re.sub("(^|\n)(COLLISIONENERGY:) ([0-9]*)((?: )?)((?:\()?)([a-zA-Z]*)((?:\))?)\n",f"\nCOLLISIONENERGY: {collisionenergy.group(3)}\n",spectrum)
#                     # check if fragmentation mode already exist
#                     if re.search("(^|\n)(FRAGMENTATIONMODE:) (.*)\n", spectrum):
#                         fragmentationmode = re.search("(^|\n)(FRAGMENTATIONMODE:) (.*)\n", spectrum).group(3)
#                         if fragmentationmode == "None":
#                             spectrum = re.sub("FRAGMENTATIONMODE: None",f"FRAGMENTATIONMODE: {fragmentation}",spectrum)
#
#     return spectrum

def harmonize_syns(spectrum):
    """
    Replaces "$:00in-source" with "None" in the given spectrum.

    :param spectrum: The spectrum to be harmonized.
    :type spectrum: str
    :return: The harmonized spectrum.
    :rtype: str
    """
    if spectrum != None:
        if re.search("\$:00in-source",spectrum):
            spectrum = re.sub("\$:00in-source","None",spectrum)

    return spectrum

def harmonize_formula(spectrum):
    """
    Harmonizes the formula in the given spectrum.

    :param spectrum: The spectrum to harmonize.
    :return: The harmonized spectrum.
    """
    if spectrum != None:
        if re.search("FORMULA: \[(.*)\](\+|-)\n",spectrum):
            FORMULA = re.search("FORMULA: \[(.*)\](\+|-)\n",spectrum).group(1)
            spectrum = re.sub("FORMULA: (.*)\n",f"FORMULA: {FORMULA}\n",spectrum)
        elif re.search("FORMULA: (.*)(\+|-)\n",spectrum):
            FORMULA = re.search("FORMULA: (.*)(\+|-)\n",spectrum).group(1)
            spectrum = re.sub("FORMULA: (.*)\n", f"FORMULA: {FORMULA}\n",spectrum)
        elif re.search(r"FORMULA: N\\A",spectrum,flags=re.I):
            spectrum = re.sub(r"FORMULA: N\\A", "FORMULA: None", spectrum, flags=re.I)

    return spectrum

def harmonize_empties(spectrum):
    """
    :param spectrum: The spectrum to be harmonized by replacing empty values with None.
    :return: The harmonized spectrum.

    """
    if spectrum != None:
        if re.search(": \n",spectrum):
            spectrum = re.sub(": \n",": None\n",spectrum)
    return spectrum

def predicted_correction(spectrum):
    """
    :param spectrum: The spectrum to be checked for predicted correction.
    :return: The spectrum with the predicted correction flag updated.

    This method takes in a spectrum and checks if it has a predicted correction flag. If the spectrum is not None and does not match the pattern "FILENAME: MSMS_Public.*", the method removes
    * any existing predicted correction value and checks if the spectrum contains keywords such as "in-silico", "insilico", "predicted", "theoretical", or "Annotation level-3" (ignoring
    * case). If any of these keywords are found, the predicted correction flag is set to true; otherwise, it is set to false.

    If the spectrum matches the pattern "FILENAME: MSMS_Public.*", the method removes both the file name and predicted correction value from the spectrum. It then follows the same logic
    * as described above to update the predicted correction flag.

    Finally, the method returns the updated spectrum with the predicted correction flag.
    """
    if spectrum != None:
        if not re.search("FILENAME: MSMS_Public.*",spectrum):
            temp_spectrum = re.sub("PREDICTED: (.*)\n", "", spectrum)
            if re.search("in-silico|insilico|predicted|theoretical|Annotation level-3",temp_spectrum,flags=re.I):
                spectrum = re.sub("PREDICTED: .*\n","PREDICTED: true\n",spectrum)
            else:
                spectrum = re.sub("PREDICTED: .*\n", "PREDICTED: false\n", spectrum)
        else:
            temp_spectrum = re.sub("FILENAME: (.*)\n", "", spectrum)
            temp_spectrum = re.sub("PREDICTED: (.*)\n", "", temp_spectrum)
            if re.search("in-silico|insilico|predicted|theoretical|Annotation level-3", temp_spectrum, flags=re.I):
                spectrum = re.sub("PREDICTED: .*\n", "PREDICTED: true\n", spectrum)
            else:
                spectrum = re.sub("PREDICTED: .*\n", "PREDICTED: false\n", spectrum)

    return spectrum

def harmonize_db_informations(spectrum):
    """
    :param spectrum: A string representing the spectrum information.
    :return: The harmonized spectrum information.

    This method takes a spectrum string and harmonizes the database information in it. It checks for two patterns: GNPS CAS and MSMS CAS. If the pattern is found, it extracts the spectrum
    * ID and adds it to the spectrum string.

    Example usage:
    spectrum = "SPECTRUMID: None\nDB#=12345, origin=GNPS"
    harmonize_db_informations(spectrum)  # returns "SPECTRUMID: 12345\nDB#=12345, origin=GNPS"
    """
    if spectrum != None:
        # GNPS CAS
        if re.search("(DB#=(.*)(;|,)((?: )?)(origin=(.*)))",spectrum):
            matche = re.search("(DB#=(.*)(;|,)((?: )?)(origin=(.*)))",spectrum)
            spectrum_id = matche.group(2)
            # adding spectrum id
            spectrum = re.sub("SPECTRUMID: None\n",f"SPECTRUMID: {spectrum_id}\n",spectrum)
        # MSMS CAS
        elif re.search("(spec_id=(.*)(;|,)((?: )?)(origin=(.*)))", spectrum):
            matche = re.search("(spec_id=(.*)(;|,)((?: )?)(origin=(.*)))", spectrum)
            spectrum_id = matche.group(2)
            # adding spectrum id
            spectrum = re.sub("SPECTRUMID: None\n", f"SPECTRUMID: {spectrum_id}\n", spectrum)

    return spectrum

def harmonize_fields_values(spectrum):
    """
    :param spectrum: The original spectrum data that needs to be harmonized.
    :return: The harmonized spectrum data.
    """
    # spectrum = remove_no_inchikey(spectrum)
    # spectrum = remove_no_mass(spectrum)
    spectrum = harmonize_adduct(spectrum)
    spectrum = correct_ionmode(spectrum)
    # print("harmonize_fields_values",spectrum)
    spectrum = harmonize_retention_time(spectrum)
    spectrum = harmonize_ms_level(spectrum)
    # spectrum = harmonize_collisionenergy(spectrum)
    spectrum = harmonize_syns(spectrum)
    spectrum = harmonize_formula(spectrum)
    spectrum = harmonize_empties(spectrum)
    spectrum = predicted_correction(spectrum)
    spectrum = harmonize_db_informations(spectrum)

    return spectrum

def harmonize_fields_names(spectrum):
    """
    :param spectrum: The input spectrum string.
    :return: The harmonized spectrum string with standardized field names and sorted fields.
    """
    if spectrum is not None:
        expected_fields = ["FRAGHUBID","SYNON","INCHIKEY","INSTRUMENT","FORMULA","SMILES","INCHI","COMMENT","IONIZATION","RESOLUTION","FRAGMENTATIONMODE","NAME","SPECTRUMID","PRECURSORTYPE","MSLEVEL",
                           "INSTRUMENTTYPE","IONMODE","COLLISIONENERGY","PARENTMASS","PRECURSORMZ","CHARGE","NUM PEAKS","PREDICTED","RETENTIONTIME","FILENAME"]

        spectrum = re.sub("COMPOUND_NAME:","NAME:",spectrum,flags=re.I)
        spectrum = re.sub("PRECURSOR_MZ:", "PRECURSORMZ:", spectrum, flags=re.I)
        spectrum = re.sub("ADDUCT:", "PRECURSORTYPE:", spectrum, flags=re.I)
        spectrum = re.sub("INCHIKEY: \n", "INCHIKEY: None\n", spectrum)
        spectrum = re.sub("((^|\n)(.*?):) \n", "\n",spectrum)
        spectrum = re.sub("\n{2,}", "\n", spectrum)
        fields = re.finditer("(^|\n)(.*?):",spectrum)
        fields_names = [matche.group(2) for matche in fields]


        # Remove undesired punctuation
        for field in fields_names:
            spectrum = spectrum.replace(field,re.sub("\!|\(|\)|\[|\]|\{|\}|\;|\:|\'|\\|\,|\<|\>|\.|\/|\?|\@|\#|\$|\%|\^|\&|\*|\~|\+","",field))

        fields_names = [re.sub("\!|\(|\)|\[|\]|\{|\}|\;|\'|\\|\,|\<|\>|\.|\/|\?|\@|\#|\$|\%|\^|\&|\*|\~|\+","",field) for field in fields_names]

        for field in fields_names:
            if field not in expected_fields: # Si champ dans le spectre pas voulu, on le supprime
                spectrum = re.sub(rf"(^|\n){field}:.*\n","\n",spectrum)
        for field in expected_fields:
            if field not in fields_names: # Si un champ voulu est manquant, on le rajoute
                spectrum = field+": None\n"+spectrum

        # Sort fields
        SPECTRUM = ""

        SPECTRUM = SPECTRUM + re.search("((^|\n)(FRAGHUBID: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(FILENAME: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PREDICTED: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(SPECTRUMID: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(RESOLUTION: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(SYNON: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(CHARGE: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PARENTMASS: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(IONIZATION: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(MSLEVEL: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(FRAGMENTATIONMODE: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(NAME: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PRECURSORMZ: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(PRECURSORTYPE: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INSTRUMENTTYPE: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INSTRUMENT: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(SMILES: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INCHI: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(INCHIKEY: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(COLLISIONENERGY: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(FORMULA: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(RETENTIONTIME: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(IONMODE: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(COMMENT: (.*)\n))", spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("((^|\n)(NUM PEAKS: (.*)\n))",spectrum).group(3)
        SPECTRUM = SPECTRUM + re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum).group(2)

        return SPECTRUM