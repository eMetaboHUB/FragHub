import ijson
import json
import os
import re

def load_spectrum_list_from_msp(msp_file_path, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Load a spectrum list from a given MSP (Mass Spectral Peak) file.

    Arguments:
    msp_file_path (str): The path to the MSP file.
    progress_callback (callable, optional): A function to update the progress percentage.
    total_items_callback (callable, optional): A function to set the total number of spectra (e.g., initialize a progress bar).
    prefix_callback (callable, optional): A function to update the prefix dynamically.

    Returns:
    spectrum_list (List[str]): The list of spectra read from the file, where each spectrum is represented as a string.
    """
    # Étape 1 : Compter le nombre de spectres
    num_spectra = sum(1 for line in open(msp_file_path, 'r', encoding="UTF-8") if line.strip() == '')

    # Étape 2 : Exécuter le callback pour initialiser les totaux si défini
    if total_items_callback:
        total_items_callback(num_spectra, 0)  # Total items = num_spectra, Completed = 0

    # Mise à jour dynamique avec le préfixe si défini
    if prefix_callback:
        prefix_callback(f"loading [{os.path.basename(msp_file_path)}]:")

    if item_type_callback:
        item_type_callback("spectra")

    # Étape 3 : Obtenir le nom du fichier
    filename = os.path.basename(msp_file_path)

    # Étape 4 : Initialiser une liste pour stocker les spectres
    spectrum_list = []

    # Étape 5 : Initialiser un buffer pour contenir temporairement les données
    buffer = [f"FILENAME: {filename}\n"]  # initialisation du buffer avec le nom du fichier

    # Étape 6 : Ouvrir le fichier et lire ligne par ligne
    with open(msp_file_path, 'r', encoding="UTF-8") as file:

        # Initialisez un compteur pour suivre les progrès
        processed_spectra = 0

        for line in file:

            # Si la ligne est vide, cela signifie que le spectre est terminé
            if line.strip() == '':
                # Vérifier si le buffer contient des données
                if buffer:
                    # Joignez toutes les données du buffer en une seule chaîne
                    spectrum = '\n'.join(buffer)
                    spectrum = re.sub(
                        r"FILENAME: .*\n",
                        f"FILENAME: {filename}\n",
                        spectrum,
                        flags=re.IGNORECASE
                    )
                    # Ajouter le spectre à la liste
                    spectrum_list.append(spectrum)

                    # Réinitialiser le buffer
                    buffer = [f"FILENAME: {filename}\n"]

                    # Mise à jour de la progression
                    processed_spectra += 1

                    if progress_callback:
                        progress_callback(processed_spectra)  # Met à jour directement avec le spectre traité

            # Si la ligne n'est pas vide, ajoutez-la au buffer
            else:
                buffer.append(line.strip())

    # Retourne la liste des spectres extraits
    return spectrum_list


def load_spectrum_list_from_mgf(mgf_file_path, progress_callback=None, total_items_callback=None, prefix_callback=None,
                                item_type_callback=None):
    """
    Load a spectrum list from a given MGF (Mascot Generic Format) file, with support for progress callbacks.
    :param mgf_file_path: The path to the MGF file.
    :param progress_callback: A function to update the progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items processed (optional).
    :return: The list of spectra read from the file. Each spectrum is represented as a string.
    """
    # Count the total number of spectra (lines with 'END IONS')
    num_spectra = sum(1 for line in open(mgf_file_path, 'r', encoding="UTF-8") if line.strip() == 'END IONS')

    # Extract file name
    filename = os.path.basename(mgf_file_path)

    # Set total items via callback if provided
    if total_items_callback:
        total_items_callback(num_spectra, 0)  # total = num_spectra, completed = 0

    # Update the prefix dynamically via callback if provided
    if prefix_callback:
        prefix_callback(f"Loading [{filename}]:")

    # Specify the type of items being processed via callback if provided
    if item_type_callback:
        item_type_callback("spectra")

    # Initialize variables for parsing
    spectrum_list = []
    buffer = [f"FILENAME={filename}"]
    processed_items = 0

    with open(mgf_file_path, 'r', encoding="UTF-8") as file:
        for line in file:
            # Detect end of a spectrum
            if line.strip() == 'END IONS':
                if buffer:
                    # Process the buffer into a single spectrum string
                    spectrum = '\n'.join(buffer)
                    spectrum = re.sub(r"FILENAME=.*\n", f"FILENAME={filename}\n", spectrum, flags=re.IGNORECASE)
                    spectrum_list.append(spectrum)
                    buffer = [f"FILENAME={filename}"]  # Reset the buffer for the next spectrum

                # Update progress via callback if provided
                processed_items += 1
                if progress_callback:
                    progress_callback(processed_items)
            else:
                # Accumulate lines in the buffer
                buffer.append(line.strip())

    # Return the list of spectra
    return spectrum_list


def load_spectrum_list_json(json_file_path, progress_callback=None, total_items_callback=None, prefix_callback=None,
                            item_type_callback=None):
    """
    This function loads a list of spectra from a given JSON file, with support for progress callbacks based on the number of items.

    :param json_file_path: A string representing the path to the JSON file containing the spectra.
    :param progress_callback: A function to update the progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items processed (optional).
    :return: A generator yielding each spectrum (a dictionary) from the JSON file, one at a time.
    """
    # Extract the filename from the given path
    filename = os.path.basename(json_file_path)

    # Dynamically update the prefix if provided
    if prefix_callback:
        prefix_callback(f"loading [{filename}]:")

    # Update the item type being processed if provided
    if item_type_callback:
        item_type_callback("spectra")

    # Open the JSON file for reading
    with open(json_file_path, 'r', encoding="UTF-8") as file:
        # Use ijson to count the total number of items upfront
        total_items = sum(1 for _ in ijson.items(file, 'item'))

    # Re-set the total number of items via the total_items_callback
    if total_items_callback:
        total_items_callback(total_items, 0)  # total = total_items, completed = 0

    # Re-open the JSON file to retrieve items using a generator
    with open(json_file_path, 'r', encoding="UTF-8") as file:
        # Create a generator to yield each 'item' in the JSON file using ijson
        spectra = ijson.items(file, 'item')

        # Progress variable to track the number of processed items
        processed_items = 0

        # Loop over each spectrum in the file
        for spectrum in spectra:
            # Add the originating filename to the spectrum dictionary
            spectrum["filename"] = filename

            # Update progress if progress_callback is provided
            processed_items += 1
            if progress_callback:
                progress_callback(processed_items)

            # Yield the current spectrum, allowing the function to be used as a generator
            yield spectrum




def load_spectrum_list_json_2(json_file_path, progress_callback=None, total_items_callback=None, prefix_callback=None,
                              item_type_callback=None):
    """
    Cette fonction charge une liste de spectres à partir d'un fichier JSON donné, avec gestion de la progression
    basée sur le nombre d'éléments via des callbacks.

    :param json_file_path: Une chaîne représentant le chemin vers le fichier JSON contenant les spectres.
    :param progress_callback: Une fonction pour mettre à jour la progression (optionnel).
    :param total_items_callback: Une fonction pour définir le total des éléments à traiter (optionnel).
    :param prefix_callback: Une fonction pour définir dynamiquement le préfixe de l'opération (optionnel).
    :param item_type_callback: Une fonction pour spécifier le type d'éléments traités (optionnel).
    :return: Un générateur retournant chaque spectre (dictionnaire) du fichier JSON, un à la fois.
    """

    # Extraire le nom de fichier à partir du chemin donné
    filename = os.path.basename(json_file_path)

    # Définir le préfixe dynamique via callback si fourni
    if prefix_callback:
        prefix_callback(f"loading [{filename}]:")

    # Spécifier le type d'éléments traités via callback
    if item_type_callback:
        item_type_callback("spectra")

    # Calculer le nombre total de spectres dans le fichier pour définir `total_items`
    with open(json_file_path, 'r', encoding="UTF-8") as file:
        total_items = sum(1 for line in file if line.strip())  # Compter les lignes non vides

    # Définir le total via `total_items_callback` si fourni
    if total_items_callback:
        total_items_callback(total_items, 0)  # total = total_items, completed = 0

    # Réinitialiser le fichier pour commencer la lecture des éléments
    with open(json_file_path, 'r', encoding="UTF-8") as file:
        processed_items = 0  # Variable pour suivre le nombre d'éléments traités

        # Lire ligne par ligne
        for line in file:
            if line.strip():  # Vérifiez que la ligne n'est pas vide
                # Analyser l'objet JSON à partir de la ligne actuelle
                spectrum = json.loads(line.strip())

                # Ajouter le nom du fichier d'origine au spectre
                spectrum["filename"] = filename

                # Mettre à jour la progression via `progress_callback` si défini
                processed_items += 1
                if progress_callback:
                    progress_callback(processed_items)

                # Produire le spectre actuel
                yield spectrum
