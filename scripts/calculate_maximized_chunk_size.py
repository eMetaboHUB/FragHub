import psutil
import math
import os

def calculate_maximized_chunk_size(data_size: int, estimate_item_size: int = 1024) -> int:
    """
    Calcule la taille optimale des chunks en fonction des ressources système pour utiliser pleinement le CPU et la mémoire.

    :param data_size: Nombre total d'éléments dans les données (par exemple, taille de `spectrum_list`).
    :param estimate_item_size: Taille approximative en mémoire d'un élément (spectre) en octets. Par défaut : 1 KB.
    :return: Taille optimale d'un chunk pour maximiser l'utilisation des ressources.
    """

    # Obtenez le nombre maximal de threads logiques (CPU disponibles)
    cpu_count = os.cpu_count()  # Nombre de cœurs logiques
    # Mémoire disponible (en octets)
    available_memory = psutil.virtual_memory().available

    # Calcul maximal : combien d'éléments peuvent tenir dans la mémoire par cœur
    max_chunk_memory_per_cpu = available_memory // cpu_count  # Mémoire disponible pour chaque cœur
    chunk_items_by_memory = max_chunk_memory_per_cpu // estimate_item_size  # Nombre d'items par chunk (taille mémoire)

    # Nombre d'items pour une répartition "équilibrée"
    chunk_items_by_cpu = math.ceil(data_size / cpu_count)  # Diviser équitablement les données par cœur

    # Taille optimale est le minimum entre la mémoire disponible et la charge équilibrée par cœur
    chunk_size = min(chunk_items_by_memory, chunk_items_by_cpu)

    return max(1, chunk_size)  # Toujours retourner au moins 1