
import re

def remove_no_inchikey(spectrum):
    """
    .. function:: remove_no_inchikey(spectrum)

        This method removes spectra with no InChIKey.

        :param spectrum: The input spectrum.
        :type spectrum: any

        :return: The spectrum with InChIKey if it exists, or None if no InChIKey is found.
        :rtype: any or None

    """
    if spectrum != None:
        if re.search("INCHIKEY: \n|INCHIKEY: None\n",spectrum):
            return None
        else:
            return spectrum

def harmonize_adduct(spectrum):
    """
    :param spectrum: The spectrum string
    :return: The modified spectrum string with harmonized adduct
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
    Removes spectra with no mass information from the input spectrum.

    :param spectrum: The input spectrum.
    :return: The modified spectrum with no mass information.
    """
    if spectrum != None:
        if re.search("PRECURSORMZ: None\n",spectrum) and re.search("PARENTMASS: None\n",spectrum):
            return None
        else:
            return spectrum

def correct_ionmode(spectrum):
    """

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
    Harmonize the retention time of a spectrum.

    :param spectrum: The input spectrum.
    :type spectrum: str
    :return: The spectrum with harmonized retention time.
    :rtype: str
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
    :param spectrum: A string representing the input spectrum
    :return: The spectrum with harmonized MS level

    This method takes a spectrum string as input and replaces the MS level values in the spectrum with harmonized values. The harmonization is done by replacing specific MS level patterns
    * with their corresponding values.

    The MS level patterns that are replaced are:
    - "MSLEVEL: MS" is replaced with "MSLEVEL: 1"
    - "MSLEVEL: MS1" is replaced with "MSLEVEL: 1"
    - "MSLEVEL: MS2" is replaced with "MSLEVEL: 2"
    - "MSLEVEL: MS3" is replaced with "MSLEVEL: 3"
    - "MSLEVEL: MS4" is replaced with "MSLEVEL: 4"
    - "MSLEVEL: 2-MS4 Composite" is replaced with "MSLEVEL: 4"
    - "MSLEVEL: 2-MS5 Composite" is replaced with "MSLEVEL: 5"

    The harmonized spectrum is then returned as output.
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
    Replace "$:00in-source" with "None" in the given spectrum.

    :param spectrum: The input spectrum string.
    :return: The spectrum string with "$:00in-source" replaced by "None".
    """
    if spectrum != None:
        if re.search("\$:00in-source",spectrum):
            spectrum = re.sub("\$:00in-source","None",spectrum)

    return spectrum

def harmonize_formula(spectrum):
    """
    Method to harmonize the formula in the given spectrum.

    :param spectrum: The spectrum string.
    :return: The harmonized spectrum string.
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
    :param spectrum: The input spectrum string that may contain empty values represented by ": \n".
    :return: The spectrum string with all empty values replaced by ": None\n".
    """
    if spectrum != None:
        if re.search(": \n",spectrum):
            spectrum = re.sub(": \n",": None\n",spectrum)
    return spectrum

def predicted_correction(spectrum):
    """
    :param spectrum: The input spectrum string.
    :return: The corrected spectrum string.

    This method checks if the spectrum is predicted or not and corrects the spectrum string accordingly. If the spectrum is not None, it first checks if the spectrum contains the string
    * "FILENAME: MSMS_Public". If it doesn't, it removes the line starting with "PREDICTED: " from the spectrum string. Then it checks if the spectrum contains any of the following keywords
    *: "in-silico", "insilico", "predicted", "theoretical", "Annotation level-3" (case-insensitive). If it does, it replaces the line starting with "PREDICTED: " with "PREDICTED: true".
    * Otherwise, it replaces the line with "PREDICTED: false".

    If the spectrum contains "FILENAME: MSMS_Public", it removes the line starting with "FILENAME: " and "PREDICTED: " from the spectrum string, and performs the same check for keywords
    *.

    Finally, it returns the corrected spectrum string.
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
    Harmonizes the database information in the given spectrum.

    :param spectrum: The spectrum with database information.
    :return: The modified spectrum.
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
    :param spectrum: the spectrum data to be harmonized
    :return: the harmonized spectrum data

    This method takes a spectrum data as input and performs a series of operations to harmonize the values of various fields in the spectrum. The input spectrum is modified in-place and
    * the modified spectrum is returned.
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
    :param spectrum: the spectrum to be harmonized
    :return: the harmonized spectrum with fields in a specific order and modified field names
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