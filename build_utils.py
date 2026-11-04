"""Outils communs aux scripts de build FragHub (build_backend.py, build_frontend.py, build_all.py).

Ce fichier n'est pas à lancer directement.
"""
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUST_DIR = ROOT / "scripts" / "fraghub_rust"
SCRIPTS_DIR = ROOT / "scripts"
GUI_DIR = ROOT / "GUI"
VENV_DIR = ROOT / ".venv"
IS_WINDOWS = sys.platform == "win32"

DRY_RUN = "--dry-run" in sys.argv

_BLUE, _GREEN, _RED, _BOLD, _RESET = "\033[94m", "\033[92m", "\033[91m", "\033[1m", "\033[0m"
if IS_WINDOWS:
    os.system("")  # active les couleurs ANSI dans la console Windows


class BuildError(Exception):
    pass


def step(title):
    print(f"\n{_BOLD}{_BLUE}=== {title} ==={_RESET}", flush=True)


def ok(msg):
    print(f"{_GREEN}✔ {msg}{_RESET}", flush=True)


def fail(msg):
    print(f"\n{_RED}{_BOLD}✘ {msg}{_RESET}", flush=True)


def venv_bin(name):
    """Chemin d'un exécutable du venv du projet (.venv/bin/x ou .venv\\Scripts\\x.exe)."""
    folder = VENV_DIR / ("Scripts" if IS_WINDOWS else "bin")
    exe = folder / (name + (".exe" if IS_WINDOWS else ""))
    return exe


def ensure_venv():
    """Crée/synchronise le venv avec uv s'il est absent."""
    if venv_bin("python").exists():
        return
    uv = shutil.which("uv")
    if not uv:
        raise BuildError(
            "Environnement Python introuvable (.venv) et 'uv' n'est pas installé.\n"
            "Installez uv (https://docs.astral.sh/uv/) puis relancez, ou lancez 'uv sync' à la racine du projet."
        )
    run([uv, "sync"], cwd=ROOT, title="Création de l'environnement Python (uv sync)")


# Dossiers où cargo / npm sont souvent installés, au cas où le PATH ne les contient pas
# (script lancé depuis PyCharm, VS Code, un double-clic…).
_EXTRA_PATHS = [Path.home() / ".cargo" / "bin", Path("/opt/homebrew/bin"), Path("/usr/local/bin")]


def require_tool(name, hint):
    path = shutil.which(name)
    if not path:
        for folder in _EXTRA_PATHS:
            candidate = shutil.which(name, path=str(folder))
            if candidate:
                # rend aussi le dossier visible pour les sous-processus (maturin appelle cargo, npm appelle node)
                os.environ["PATH"] = str(folder) + os.pathsep + os.environ.get("PATH", "")
                path = candidate
                break
    if not path:
        if DRY_RUN:
            print(f"  (attention : '{name}' introuvable sur cette machine)")
            return name
        raise BuildError(f"'{name}' introuvable. {hint}")
    return path


def run(cmd, cwd, title=None, env_extra=None):
    """Lance une commande, affiche sa sortie en direct, s'arrête si elle échoue."""
    cmd = [str(c) for c in cmd]
    if title:
        print(f"→ {title}")
    print(f"  $ (cd {cwd}) {' '.join(cmd)}", flush=True)
    if DRY_RUN:
        return
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    start = time.time()
    result = subprocess.run(cmd, cwd=cwd, env=env)
    if result.returncode != 0:
        raise BuildError(f"La commande a échoué (code {result.returncode}) : {' '.join(cmd)}")
    ok(f"terminé en {format_duration(time.time() - start)}")


def format_duration(seconds):
    m, s = divmod(int(seconds), 60)
    return f"{m} min {s:02d} s" if m else f"{s} s"


def sync_project_version():
    """Recopie la version de pyproject.toml dans tout le projet avant de builder (voir sync_version.py)."""
    import sync_version
    version = sync_version.read_source_version()
    step(f"Version du projet : {version}")
    _, problems = sync_version.sync(version, check_only=DRY_RUN, quiet=True)
    if problems:
        raise BuildError("Synchronisation de la version impossible (voir les ⚠ ci-dessus).")
    ok("version vérifiée (rien écrit en --dry-run)" if DRY_RUN else "version synchronisée")


def main_wrapper(func):
    """Exécute une fonction de build avec gestion d'erreur et chronomètre global."""
    if "-h" in sys.argv or "--help" in sys.argv:
        print((func.__doc__ or "").strip())
        print("\nOption : --dry-run  affiche les commandes sans les exécuter.")
        return 0
    start = time.time()
    if DRY_RUN:
        print("(mode --dry-run : les commandes sont affichées, pas exécutées)")
    try:
        sync_project_version()
        func()
    except BuildError as e:
        fail(str(e))
        return 1
    except KeyboardInterrupt:
        fail("Interrompu par l'utilisateur.")
        return 130
    print(f"\n{_GREEN}{_BOLD}★ Build terminé en {format_duration(time.time() - start)} ★{_RESET}")
    return 0
