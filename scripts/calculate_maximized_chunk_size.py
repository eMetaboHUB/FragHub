import globals_vars
import math
import sys


def calculate_maximized_chunk_size(data_list: list) -> int:
    """
    Calcule la taille optimale des chunks en fonction des ressources système pour utiliser pleinement le CPU et la mémoire.

    :param data_list: Liste d'éléments à traiter (ex: spectrum_list).
                      La taille du premier élément sera utilisée pour estimer la taille moyenne des items.
    :return: Taille optimale d'un chunk pour maximiser l'utilisation des ressources.
    """
    # Obtenez le nombre maximal de threads logiques (CPU disponibles)
    cpu_count = globals_vars.cpu_count
    # Mémoire disponible (en octets)
    available_memory = globals_vars.available_memory

    # Estimer la taille moyenne d'un élément en mémoire (en octets) d'après le premier élément
    if not data_list:
        raise ValueError("La liste des données est vide. Impossible de calculer le chunk optimal.")

    estimate_item_size = sys.getsizeof(data_list[0])  # Taille mémoire du premier élément (en octets)

    # Calcul maximal : combien d'éléments peuvent tenir dans la mémoire par cœur
    max_chunk_memory_per_cpu = available_memory // cpu_count  # Mémoire disponible pour chaque cœur
    chunk_items_by_memory = max_chunk_memory_per_cpu // estimate_item_size  # Nombre possible d'items par chunk (mémoire max)

    # Nombre d'items par répartition équilibrée (nombre total d'éléments divisé par le nombre de CPU)
    chunk_items_by_cpu = math.ceil(len(data_list) / cpu_count)  # Diviser équitablement entre les cœurs

    # Taille optimale est le minimum entre la mémoire disponible (pour chaque CPU) et la répartition équilibrée
    chunk_size = min(chunk_items_by_memory, chunk_items_by_cpu)

    return max(1, chunk_size)  # Toujours retourner au moins 1 élément par chunk
