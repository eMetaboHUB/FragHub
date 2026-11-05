"""Recompile le back-end FragHub : moteur Rust (maturin) puis exécutable Python (PyInstaller).

Usage :  python build_backend.py            (Rust + PyInstaller)
         python build_backend.py --rust-only (juste le module Rust, pour tester en dev)
         python build_backend.py --dry-run   (affiche les commandes sans rien lancer)

Résultat : scripts/dist/FragHub_Backend/
"""
import subprocess
import sys

from build_utils import (BuildError, DRY_RUN, IS_WINDOWS, RDKIT_PREFIX, RUST_DIR, SCRIPTS_DIR, VENV_DIR, ensure_venv,
                         main_wrapper, require_tool, run, step, venv_bin)

# Commande d'installation de RDKit C++, rappelée si on ne le trouve pas
_RDKIT_INSTALL = ("mamba create -p \"{0}\" -c conda-forge librdkit-dev=2026.09.1 libboost-headers=1.90.0"
                  "\n  (micromamba s'utilise avec la même syntaxe ; voir la section 7 de la documentation)"
                  .format(RDKIT_PREFIX))


def rdkit_lib_dirs():
    """Dossiers de RDKit C++ : (en-têtes, DLL ou .so). Erreur claire s'il n'est pas installé."""
    base = RDKIT_PREFIX / "Library" if IS_WINDOWS else RDKIT_PREFIX
    include = base / "include" / "rdkit"
    libs = base / ("bin" if IS_WINDOWS else "lib")
    if not include.is_dir():
        raise BuildError(f"RDKit C++ introuvable dans {RDKIT_PREFIX}.\n"
                         f"Installez-le avec :\n  {_RDKIT_INSTALL}\n"
                         "ou indiquez son emplacement avec la variable d'environnement RDKIT_PREFIX.")
    return include, libs


def build_rust():
    step("Back-end 1/2 : moteur Rust (maturin develop --release)")
    require_tool("cargo", "Installez Rust : https://rustup.rs")
    rdkit_lib_dirs()
    ensure_venv()
    maturin = venv_bin("maturin")
    cmd = [maturin, "develop", "--release"] if maturin.exists() else [venv_bin("python"), "-m", "maturin", "develop", "--release"]
    # VIRTUAL_ENV : indique à maturin d'installer le module dans le .venv du projet
    env = {"VIRTUAL_ENV": str(VENV_DIR), "RDKIT_PREFIX": str(RDKIT_PREFIX)}
    if sys.platform == "darwin":
        # macOS : pas de "strip" après l'édition de liens. Avec les outils Apple récents (SDK 27),
        # le strip par défaut de Cargo laisse la table des symboles mal alignée et macOS refuse de charger
        # le module ("mis-aligned LINKEDIT string pool"). Sans effet sur le code compilé ni sur les résultats.
        env["CARGO_PROFILE_RELEASE_STRIP"] = "false"
    run(cmd, cwd=RUST_DIR, title="Compilation et installation de fraghub_rust dans .venv", env_extra=env)
    if sys.platform == "darwin" and not DRY_RUN:
        # Vérifie tout de suite que macOS accepte de charger le module (sinon l'erreur n'apparaîtrait qu'au lancement)
        check = subprocess.run([str(venv_bin("python")), "-c", "import fraghub_rust"], capture_output=True, text=True)
        if check.returncode != 0:
            raise BuildError("Le module Rust compilé ne se charge pas :\n" + check.stderr.strip()[-1500:])


def build_pyinstaller():
    step("Back-end 2/2 : exécutable Python (PyInstaller)")
    ensure_venv()
    spec = SCRIPTS_DIR / "FragHub_Backend.spec"
    if not spec.exists():
        raise BuildError(f"Fichier introuvable : {spec}")
    _, rdkit_libs = rdkit_lib_dirs()
    # Le .spec ajoute les DLL de RDKit C++ au PATH pour que PyInstaller les copie dans l'exécutable
    run([venv_bin("python"), "-m", "PyInstaller", "-y", spec.name], cwd=SCRIPTS_DIR,
        title="Création de scripts/dist/FragHub_Backend",
        env_extra={"RDKIT_PREFIX": str(RDKIT_PREFIX)})
    # Vérification : sans les DLL de RDKit, l'exécutable ne démarrerait pas sur une machine sans mamba
    internal = SCRIPTS_DIR / "dist" / "FragHub_Backend" / "_internal"
    if IS_WINDOWS and not DRY_RUN and not any(internal.glob("RDKitGraphMol*.dll")):
        raise BuildError(f"Les DLL de RDKit C++ n'ont pas été copiées dans {internal} (cherchées dans {rdkit_libs}).")
    if sys.platform == "darwin" and not DRY_RUN:
        check_macos_bundle(internal, rdkit_libs)


def check_macos_bundle(internal, rdkit_libs):
    """macOS : les bibliothèques de RDKit C++ doivent être copiées dans l'exécutable, et aucun binaire
    ne doit encore pointer vers le dossier de RDKit C++ (sinon l'exécutable ne marcherait que sur cette machine)."""
    if not any(internal.glob("libRDKitGraphMol*.dylib")):
        raise BuildError(f"Les bibliothèques de RDKit C++ n'ont pas été copiées dans {internal} (cherchées dans {rdkit_libs}).")
    prefix = str(RDKIT_PREFIX)
    suspects = []
    for f in internal.rglob("*"):
        if f.is_symlink() or not f.is_file() or f.suffix not in (".so", ".dylib"):
            continue
        # otool -l liste les chemins de chargement (LC_RPATH, LC_LOAD_DYLIB) d'un binaire Mach-O
        out = subprocess.run(["otool", "-l", str(f)], capture_output=True, text=True).stdout
        if prefix in out:
            suspects.append(f.relative_to(internal).as_posix())
    if suspects:
        raise BuildError("Ces binaires de l'exécutable pointent encore vers " + prefix + " :\n  "
                         + "\n  ".join(suspects[:10]))
    print(f"  RDKit C++ bien embarqué, aucun chemin vers {prefix}")


def build_backend():
    """Recompile le back-end FragHub : Rust (maturin) puis PyInstaller.

Usage :  python build_backend.py [--rust-only] [--dry-run]"""
    build_rust()
    if "--rust-only" not in sys.argv:
        build_pyinstaller()


if __name__ == "__main__":
    sys.exit(main_wrapper(build_backend))
