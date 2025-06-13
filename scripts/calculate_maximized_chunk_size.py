import scripts.globals_vars
import math
import sys


def calculate_maximized_chunk_size(data_list: list) -> int:
    """
    Calculates the optimal chunk size based on system resources to fully utilize the CPU and memory.

    :param data_list: List of elements to process (e.g., spectrum_list).
                      The size of the first element will be used to estimate the average size of the items.
    :return: Optimal chunk size to maximize resource utilization.
    """
    # Get the maximum number of logical threads (available CPUs)
    cpu_count = scripts.globals_vars.cpu_count
    # Available memory (in bytes)
    available_memory = scripts.globals_vars.available_memory

    # Estimate the average size of an element in memory (in bytes) based on the first element

    estimate_item_size = sys.getsizeof(data_list[0])  # Memory size of the first element (in bytes)

    # Maximum calculation: how many elements can fit in memory per core
    max_chunk_memory_per_cpu = available_memory // cpu_count  # Memory available per core
    chunk_items_by_memory = max_chunk_memory_per_cpu // estimate_item_size  # Possible number of items per chunk (maximum memory)

    # Number of items based on balanced distribution (total number of elements divided by the number of CPUs)
    chunk_items_by_cpu = math.ceil(len(data_list) / cpu_count)  # Distribute items evenly across cores

    # Optimal size is the minimum between available memory (per CPU) and balanced distribution
    chunk_size = min(chunk_items_by_memory, chunk_items_by_cpu)

    return max(1, chunk_size)  # Always return at least 1 item per chunk
