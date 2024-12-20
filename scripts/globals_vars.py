import psutil
import os


global available_memory
available_memory = psutil.virtual_memory().available

global cpu_count
cpu_count = os.cpu_count()  # Nombre de cœurs logiques