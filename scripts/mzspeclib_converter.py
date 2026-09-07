import os
import sys
import traceback
from mzspeclib.backends.msp import MSPSpectralLibrary
from mzspeclib.backends.json import JSONSpectralLibraryWriter
import time

def convert_msp_to_mzspeclib_json(input_msp, output_json, progress_cb=None):
    try:
        lib = MSPSpectralLibrary(input_msp)
        writer = JSONSpectralLibraryWriter(output_json)

        count = 0
        skipped = 0
        # On itere manuellement sur l'index plutot que "for spec in lib" :
        # un spectre individuel invalide (ex: IONMODE "NOT FOUND" que mzspeclib
        # ne sait pas convertir en polarite) ne doit pas interrompre tout le fichier.
        for record in lib.index:
            try:
                spec = lib.get_spectrum(record.number)
            except Exception as spec_e:
                skipped += 1
                print(f"[mzSpecLib] Spectre ignore (#{record.number}) dans {os.path.basename(input_msp)}: {spec_e}")
                if progress_cb:
                    progress_cb(1)
                continue

            writer.write_spectrum(spec)
            count += 1
            if progress_cb:
                progress_cb(1)

        writer.close()
        if skipped:
            print(f"[mzSpecLib] {skipped} spectre(s) ignore(s) dans {os.path.basename(input_msp)} (non convertibles).")
        return count
    except Exception as e:
        traceback.print_exc()
        raise e

def convert_all_msp(output_directory, step_cb=None, prefix_cb=None, progress_cb=None, total_items_cb=None):
    msp_dir = os.path.join(output_directory, "MSP")
    json_dir = os.path.join(output_directory, "mzSpecLib")
    
    if not os.path.exists(msp_dir):
        return
        
    msp_files = []
    for root, dirs, files in os.walk(msp_dir):
        for f in files:
            if f.endswith('.msp'):
                msp_files.append(os.path.join(root, f))
    
    if not msp_files:
        return
        
    if step_cb:
        step_cb("Writing mzSpecLib JSON")

    total_count = 0
    # On calcule le total des spectres
    # (ce bloc prend un peu de temps mais permet la barre de progression)
    for msp_path in msp_files:
        with open(msp_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.upper().startswith("NAME:"):
                    total_count += 1
    
    if total_items_cb:
        total_items_cb(total_count)

    global_count = 0
    def global_progress_cb(count):
        nonlocal global_count
        global_count += 1
        if progress_cb and global_count % 100 == 0:
            progress_cb(global_count)

    for msp_path in msp_files:
        rel_path = os.path.relpath(msp_path, msp_dir)
        json_path = os.path.join(json_dir, rel_path)
        json_path = json_path.replace('.msp', '.mzSpecLib.json')
        
        os.makedirs(os.path.dirname(json_path), exist_ok=True)
        
        if prefix_cb:
            prefix_cb(f"Writing {os.path.basename(json_path)} to disk...")
            
        convert_msp_to_mzspeclib_json(msp_path, json_path, progress_cb=global_progress_cb)

    if progress_cb:
        progress_cb(global_count)
