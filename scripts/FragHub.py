import asyncio
import socketio
import uvicorn
import multiprocessing
import traceback
import mzspeclib_converter
import time
from contextlib import asynccontextmanager
from fastapi import FastAPI, BackgroundTasks
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, ConfigDict

import sys
import os
import fraghub_rust
import rdkit_worker  # Ensure PyInstaller bundles it for Rust to use

parameters_dict = {}

if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


is_cli_mode = "--cli" in sys.argv
is_quiet_mode = "--quiet" in sys.argv

# Redirection forcée de la sortie standard et d'erreur vers un fichier (UNIQUEMENT si packagé avec PyInstaller)
if getattr(sys, 'frozen', False) and not is_cli_mode:
    log_path = os.path.join(os.path.expanduser("~"), "fraghub_debug.txt")
    sys.stdout = open(log_path, 'w')
    sys.stderr = sys.stdout

if not is_quiet_mode:
    print(f"--- Démarrage de FragHub ---")
    print(f"CWD: {os.getcwd()}")
    print(f"sys.argv: {sys.argv}")

# Variables globales
loop = None
last_emit_time = 0
EMIT_THROTTLE = 0.05  # 20ms pour une fluidité maximale

# Configuration Socket.IO et FastAPI
sio = socketio.AsyncServer(async_mode='asgi', cors_allowed_origins='*')

@asynccontextmanager
async def lifespan(app: FastAPI):
    global loop
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = asyncio.get_event_loop()
    yield

app = FastAPI(lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

socket_app = socketio.ASGIApp(sio, other_asgi_app=app)

class FragHubParams(BaseModel):
    input_directory: list
    input_db_names: dict = Field(default_factory=dict)
    output_directory: str
    normalize_intensity: float
    remove_peak_above_precursormz: float
    check_minimum_peak_requiered: float
    check_minimum_peak_requiered_n_peaks: float
    reduce_peak_list: float
    reduce_peak_list_max_peaks: float
    remove_spectrum_under_entropy_score: float
    remove_spectrum_under_entropy_score_value: float
    keep_mz_in_range: float
    keep_mz_in_range_from_mz: float
    keep_mz_in_range_to_mz: float
    check_minimum_of_high_peaks_requiered: float
    check_minimum_of_high_peaks_requiered_intensity_percent: float
    check_minimum_of_high_peaks_requiered_no_peaks: float
    calculate_de_novo: float
    de_novo_ppm_tolerance: float
    csv: float
    msp: float
    json_enabled: float = Field(alias='json')
    mzspeclib_json: float = Field(alias='mzspeclib_json', default=1.0)
    reset_updates: float

    model_config = ConfigDict(populate_by_name=True)

def emit_to_frontend(event, data):
    global loop, last_emit_time

    # Throttling UNIQUEMENT sur progress — tous les autres événements
    # (total_items, prefix, step, completion, deletion…) passent toujours.
    # Throttler total_items causait des resets perdus → barre qui saute à 100%.
    if event == 'progress':
        current_time = time.time()
        if (current_time - last_emit_time) < EMIT_THROTTLE:
            return
        last_emit_time = current_time

    if loop:
        try:
            asyncio.run_coroutine_threadsafe(sio.emit(event, data), loop)
        except Exception:
            pass

# --- CALLBACKS ---
deferred_deletion = None
deferred_report = False
current_total_items = 0  # <-- Nouvelle variable globale

def progress_callback(*args):
    global last_emit_time, current_total_items
    if global_stop_flag:
        raise Exception("Process stopped by user.")
    if args:
        # Si on atteint 100%, on remet le chrono à zéro pour FORCER l'envoi
        if args[0] >= current_total_items:
            last_emit_time = 0

        emit_to_frontend('progress', args[0])

def total_items_callback(*args):
    global last_emit_time, current_total_items
    last_emit_time = 0
    if args:
        current_total_items = args[0]  # <-- On enregistre le maximum
        emit_to_frontend('total_items', args[0])

def prefix_callback(*args):
    global last_emit_time
    last_emit_time = 0
    if args: emit_to_frontend('prefix', args[0])

def item_type_callback(*args):
    if args: emit_to_frontend('item_type', args[0])

def step_callback(*args):
    global deferred_report
    if args:
        if args[0] == "Generating Report" and parameters_dict.get("mzspeclib_json", 0.0) == 1.0:
            deferred_report = True
            return
        emit_to_frontend('step', args[0])

def deletion_callback(*args):
    global deferred_deletion
    if args:
        if parameters_dict.get("mzspeclib_json", 0.0) == 1.0:
            deferred_deletion = args[0]
            return
        emit_to_frontend('deletion', args[0])

def completion_callback(*args):
    if args:
        if len(args) > 1:
            emit_to_frontend('completion', {'message': args[0], 'report_path': args[1]})
        else:
            emit_to_frontend('completion', {'message': args[0], 'report_path': None})

# --- GESTION DU STOP ---
global_stop_flag = False
def get_stop_flag(): return global_stop_flag

@app.get("/stop-analysis")
async def stop_analysis():
    global global_stop_flag
    global_stop_flag = True
    return {"status": "stopped"}

@app.get("/health")
async def health_check():
    return {"status": "ok"}

@app.get("/report")
async def get_report(path: str):
    if os.path.exists(path):
        return FileResponse(path)
    return {"error": "Report not found"}

@app.get("/init-data")
async def init_data():
    fraghub_rust.load_internal_databases(BASE_DIR)
    return {"status": "loaded"}

# --- EXÉCUTION ---
def execute_main_safely():
    try:
        # We need to capture the completion data from Rust so it doesn't close the GUI early
        rust_completion_data = None
        def proxy_completion_callback(msg, report_path=""):
            nonlocal rust_completion_data
            rust_completion_data = (msg, report_path)

        # Appel direct de la fonction Rust
        fraghub_rust.main_orchestrator(
            parameters_dict,
            progress_callback,
            total_items_callback,
            prefix_callback,
            item_type_callback,
            step_callback,
            proxy_completion_callback,
            deletion_callback,
            get_stop_flag
        )
        
        # Si l'utilisateur a cliqué sur stop ou s'il y a eu une erreur précoce dans Rust
        if rust_completion_data and rust_completion_data[0] == "Process stopped by user":
            if completion_callback:
                completion_callback(rust_completion_data[0], rust_completion_data[1])
            return

        # ── Intégration de mzspeclib pour la conversion propre en mzSpecLib.json ──
        if parameters_dict.get("mzspeclib_json", 0.0) == 1.0:
            try:
                mzspeclib_converter.convert_all_msp(
                    parameters_dict["output_directory"],
                    step_cb=step_callback,
                    prefix_cb=prefix_callback,
                    progress_cb=progress_callback,
                    total_items_cb=total_items_callback
                )
            except Exception as cv_e:
                if completion_callback:
                    completion_callback(f"Warning: mzSpecLib JSON conversion failed: {cv_e}")
                return # Si on retourne ici, on n'affiche pas le TOTAL TIME

        # On envoie maintenant les étapes qu'on avait différées
        if deferred_deletion:
            emit_to_frontend('deletion', deferred_deletion)
        if deferred_report:
            emit_to_frontend('step', "Generating Report")

        # Une fois tout terminé, on envoie la vraie fin de process au GUI
        if completion_callback and rust_completion_data:
            completion_callback(rust_completion_data[0], rust_completion_data[1])

    except Exception as e:
        # Gestion de l'interruption utilisateur spécifiée dans Rust
        if str(e) == "Process stopped by user.":
            if deletion_callback:
                deletion_callback("\n-- PROCESS INTERRUPTED BY USER --")
            if completion_callback:
                completion_callback("Process stopped by user.")
        else:
            traceback.print_exc()
            if completion_callback:
                completion_callback(f"Error: {str(e)}")

@app.post("/run-analysis")
async def run_analysis(params: FragHubParams, background_tasks: BackgroundTasks):
    global global_stop_flag
    global_stop_flag = False
    params_data = params.model_dump(by_alias=True)
    for key, value in params_data.items():
        parameters_dict[key] = value
    background_tasks.add_task(execute_main_safely)
    return {"status": "started"}

if __name__ == "__main__":
    multiprocessing.freeze_support()
    
    if is_cli_mode:
        import argparse
        parser = argparse.ArgumentParser(description="FragHub CLI Mode")
        parser.add_argument("--cli", action="store_true", help="Enable CLI mode")
        parser.add_argument("--quiet", action="store_true", help="Silence all progress bars and only print final success and total time")
        parser.add_argument("--input_directory", nargs='+', required=True, help="List of input files/directories")
        parser.add_argument("--output_directory", type=str, required=True, help="Output directory path")
        
        # Filtres on/off (yes/no)
        parser.add_argument("--normalize_intensity", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--remove_peak_above_precursormz", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--check_minimum_peak_requiered", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--reduce_peak_list", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--remove_spectrum_under_entropy_score", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--keep_mz_in_range", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--check_minimum_of_high_peaks_requiered", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--calculate_de_novo", type=str, choices=['yes', 'no'], default='no')
        parser.add_argument("--csv", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--msp", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--json", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--mzspeclib_json", type=str, choices=['yes', 'no'], default='yes')
        parser.add_argument("--reset_updates", type=str, choices=['yes', 'no'], default='no')
        
        # Valeurs numériques (seuils, tolérances)
        parser.add_argument("--check_minimum_peak_requiered_n_peaks", type=float, default=3.0)
        parser.add_argument("--reduce_peak_list_max_peaks", type=float, default=500.0)
        parser.add_argument("--remove_spectrum_under_entropy_score_value", type=float, default=0.5)
        parser.add_argument("--keep_mz_in_range_from_mz", type=float, default=50.0)
        parser.add_argument("--keep_mz_in_range_to_mz", type=float, default=2000.0)
        parser.add_argument("--check_minimum_of_high_peaks_requiered_intensity_percent", type=float, default=5.0)
        parser.add_argument("--check_minimum_of_high_peaks_requiered_no_peaks", type=float, default=2.0)
        parser.add_argument("--de_novo_ppm_tolerance", type=float, default=10.0)
        
        args = parser.parse_args()
        
        # Résolution intelligente des chemins (fichiers et dossiers)
        resolved_files = []
        for path in args.input_directory:
            if os.path.isfile(path):
                resolved_files.append(os.path.abspath(path))
            elif os.path.isdir(path):
                # Parcourt tous les sous-dossiers à la recherche de fichiers valides
                for root, _, files in os.walk(path):
                    for file in files:
                        if file.lower().endswith(('.msp', '.mgf', '.csv', '.json')):
                            resolved_files.append(os.path.abspath(os.path.join(root, file)))
            else:
                print(f"[WARNING] Le chemin spécifié est introuvable : {path}")

        if not resolved_files:
            print("\n[ERROR] Aucun fichier MS valide (.msp, .mgf, .csv, .json) trouvé dans les chemins spécifiés.")
            sys.exit(1)
            
        # Hydratation du parameters_dict avec les fichiers résolus et conversion yes/no
        args_dict = vars(args)
        
        # Conversion transparente yes/no en 1.0/0.0 pour Rust
        yes_no_keys = [
            'normalize_intensity', 'remove_peak_above_precursormz', 'check_minimum_peak_requiered',
            'reduce_peak_list', 'remove_spectrum_under_entropy_score', 'keep_mz_in_range',
            'check_minimum_of_high_peaks_requiered', 'calculate_de_novo', 'csv', 'msp', 'json', 'mzspeclib_json', 'reset_updates'
        ]
        for key in yes_no_keys:
            args_dict[key] = 1.0 if args_dict[key] == 'yes' else 0.0
            
        args_dict['input_directory'] = resolved_files
        args_dict['input_db_names'] = {}
        parameters_dict.update(args_dict)
        
        if not args.quiet:
            print("\n========================================")
            print("          FragHub CLI Mode Actif          ")
            print("========================================\n")
        
        # Chargement initial des bases
        fraghub_rust.load_internal_databases(BASE_DIR)
        
        global_start_time = time.time()
        
        # Variables pour la barre de progression en console
        cli_total_items = [0]
        cli_current_prefix = [""]
        cli_start_time = [0.0]
        
        def cli_progress_callback(*cb_args):
            if args.quiet or not cb_args: return
            val = cb_args[0]
            total = cli_total_items[0]
            start_t = cli_start_time[0]
            if total > 0:
                percent = val / total
                bar_length = 30
                filled_length = int(bar_length * percent)
                bar = '█' * filled_length + '░' * (bar_length - filled_length)
                
                # Calculs de temps et vitesse
                elapsed = time.time() - start_t
                speed = val / elapsed if elapsed > 0 else 0
                remaining = (total - val) / speed if speed > 0 else 0
                
                def format_time(seconds):
                    if seconds == float('inf') or seconds < 0: return "--:--"
                    m, s = divmod(int(seconds), 60)
                    h, m = divmod(m, 60)
                    if h > 0: return f"{h:02d}:{m:02d}:{s:02d}"
                    return f"{m:02d}:{s:02d}"
                
                # Couleurs ANSI (Bleu clair = en cours, Vert = terminé)
                color = '\033[96m' if val < total else '\033[92m'
                reset = '\033[0m'
                
                stats = f"[{format_time(elapsed)} < {format_time(remaining)}, {speed:.1f} it/s]"
                
                # \033[K permet d'effacer le reste de la ligne pour éviter les artefacts visuels
                print(f"\r  {color}[{bar}]{reset} {percent*100:.1f}% | {val}/{total} {stats}\033[K", end="", flush=True)
                if val >= total:
                    print() # Nouvelle ligne à 100%
                    
        def cli_total_items_callback(*cb_args):
            if cb_args:
                cli_total_items[0] = cb_args[0]
                cli_start_time[0] = time.time() # Reset chrono au début de chaque étape
            
        def cli_prefix_callback(*cb_args):
            if args.quiet or not cb_args: return
            prefix = cb_args[0]
            cli_current_prefix[0] = prefix
            # Flèche de tâche et texte en blanc/gris
            print(f"\n\033[1;37m▶ {prefix}\033[0m")
            
        def cli_item_type_callback(*cb_args):
            pass # Non nécessaire en CLI
            
        def cli_step_callback(*cb_args):
            if args.quiet or not cb_args: return
            step = cb_args[0]
            # Étape principale en jaune et en gras
            print(f"\n\033[1;33m=== {step} ===\033[0m")
            
        def cli_completion_callback(*cb_args):
            if args.quiet or not cb_args: return
            msg = cb_args[0]
            # Message de fin en vert
            print(f"\n\033[1;32m✔ {msg}\033[0m\n")
            
        def cli_deletion_callback(*cb_args):
            if args.quiet or not cb_args: return
            msg = cb_args[0]
            # Rapports de suppression en rouge/orange
            print(f"\033[38;5;208m[REPORT] {msg}\033[0m")
            
        def cli_get_stop_flag():
            return False # Pas d'interruption via UI en mode CLI (l'utilisateur fera Ctrl+C)
            
        try:
            # Lancement de l'orchestrateur de façon synchrone dans le thread principal
            fraghub_rust.main_orchestrator(
                parameters_dict,
                cli_progress_callback,
                cli_total_items_callback,
                cli_prefix_callback,
                cli_item_type_callback,
                cli_step_callback,
                cli_completion_callback,
                cli_deletion_callback,
                cli_get_stop_flag
            )
        except Exception as e:
            print(f"\n[ERROR] Une erreur s'est produite : {e}")
            sys.exit(1)
            
        # Calcul et affichage du temps total
        total_time = time.time() - global_start_time
        def format_total(seconds):
            m, s = divmod(int(seconds), 60)
            h, m = divmod(m, 60)
            if h > 0: return f"{h}h {m}m {s}s"
            if m > 0: return f"{m}m {s}s"
            return f"{s}s"
            
        if args.quiet:
            print(f"FragHub processed successfully in {format_total(total_time)}.")
        else:
            print(f"\n\033[1;32m★ All tasks completed successfully in {format_total(total_time)} ★\033[0m\n")
            
        sys.exit(0)
    else:
        # Lancement normal du GUI (FastAPI / Uvicorn / WebSockets)
        uvicorn.run(socket_app, host="127.0.0.1", port=8000)