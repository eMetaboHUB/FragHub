from tqdm import tqdm
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
        prefix_callback(f"loading [{os.path.basename(msp_file_path)}]: ")

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

def load_spectrum_list_from_mgf(mgf_file_path):
    """
    Load a spectrum list from a given MGF (Mascot Generic Format) file.
    :param mgf_file_path: The path to the MGF file.
    :return: The list of spectra read from the file. Each spectrum is represented as a string.
    """
    # Count the total number of spectra
    num_spectra = sum(1 for line in open(mgf_file_path, 'r', encoding="UTF-8") if line.strip() == 'END IONS')

    filename = os.path.basename(mgf_file_path)
    spectrum_list = []
    buffer = [f"FILENAME={filename}"]

    with tqdm(total=num_spectra, unit="spectra", colour="green", desc="{:>70}".format(f"loading [{filename}]")) as pbar:
        with open(mgf_file_path, 'r', encoding="UTF-8") as file:
            for line in file:
                if line.strip() == 'END IONS':
                    if buffer:
                        spectrum = '\n'.join(buffer)
                        spectrum = re.sub(r"FILENAME=.*\n", f"FILENAME={filename}\n", spectrum, flags=re.IGNORECASE)
                        spectrum_list.append(spectrum)
                        buffer = [f"FILENAME={filename}"]
                    pbar.update(1)  # Update progress bar each time a spectrum is being processed
                else:
                    buffer.append(line.strip())

    # return the list of spectra
    return spectrum_list

def load_spectrum_list_json_2(json_file_path, progress_callback=None, total_items_callback=None, prefix_callback=None,
                              item_type_callback=None):
    """
    Cette fonction charge une liste de spectres à partir d'un fichier JSON donné, avec gestion de la progression via des callbacks.

    :param json_file_path: Une chaîne représentant le chemin vers le fichier JSON contenant les spectres.
    :param progress_callback: Une fonction pour mettre à jour la progression (optionnel).
    :param total_items_callback: Une fonction pour définir le total des éléments ou octets à traiter (optionnel).
    :param prefix_callback: Une fonction pour définir dynamiquement le préfixe de l'opération (optionnel).
    :param item_type_callback: Une fonction pour spécifier le type d'éléments traités (optionnel).
    :return: Un générateur retournant chaque spectre (dictionnaire) du fichier JSON, un à la fois.
    """

    # Extraire le nom de fichier à partir du chemin donné
    filename = os.path.basename(json_file_path)

    # Calculer la taille totale du fichier en octets pour suivre la progression
    total_bytes = os.path.getsize(json_file_path)

    # Définir le total via callback si fourni
    if total_items_callback:
        total_items_callback(total_bytes, 0)  # total = total_bytes, completed = 0

    # Définir le préfixe dynamique via callback si fourni
    if prefix_callback:
        prefix_callback(f"loading [{filename}]: ")

    # Spécifier le type d'éléments traités (optionnel)
    if item_type_callback:
        item_type_callback("spectra")

    # Variable pour suivre les octets traités
    processed_bytes = 0

    # Ouvrir le fichier JSON pour lecture
    with open(json_file_path, 'r', encoding="UTF-8") as file:
        # Lire ligne par ligne
        for line in file:
            if line.strip():  # Vérifiez que la ligne n'est pas vide
                # Analyser l'objet JSON à partir de la ligne actuelle
                spectrum = json.loads(line.strip())

                # Ajouter le nom du fichier d'origine au spectre
                spectrum["filename"] = filename

                # Calculer la taille en octets traités
                processed_bytes += len(line)

                # Mettre à jour le callback de progression (si défini)
                if progress_callback:
                    progress_callback(processed_bytes)

                # Produire le spectre actuel
                yield spectrum




def load_spectrum_list_json_2(json_file_path):
    """
    Cette fonction charge une liste de spectres à partir d'un fichier JSON donné.

    :param json_file_path: Une chaîne représentant le chemin vers le fichier JSON contenant les spectres.

    :return: Un générateur retournant chaque spectre (dictionnaire) du fichier JSON, un à la fois.
    """

    # Extraire le nom de fichier à partir du chemin donné
    filename = os.path.basename(json_file_path)

    # Calculer la taille totale du fichier en octets pour la progression
    total_bytes = os.path.getsize(json_file_path)

    with open(json_file_path, 'r', encoding="UTF-8") as file:
        progress = tqdm(total=total_bytes, unit="B", unit_scale=True, colour="green",
                        desc="{:>70}".format(f"loading [{filename}]"))

        # Lire ligne par ligne
        for line in file:
            if line.strip():
                try:
                    # Analyser l'objet JSON de chaque ligne
                    spectrum = json.loads(line.strip())

                    # Ajouter le nom de fichier d'origine dans le dictionnaire du spectre
                    spectrum["filename"] = filename

                    # Mettre à jour la barre de progression
                    progress.update(len(line))

                    # Produire le spectre actuel
                    yield spectrum
                except json.JSONDecodeError as e:
                    print(f"Erreur de décodage JSON: {e}")

        progress.close()
