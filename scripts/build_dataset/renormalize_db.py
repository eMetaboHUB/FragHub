# -*- coding: utf-8 -*-
"""
Re-normalise la base interne (datas/internal_databases.parquet) avec EXACTEMENT le traitement RDKit
appliqué aux spectres : neutralisation, InChI et InChIKey standard de la molécule neutralisée, SMILES
normalisé (tautomère canonique quand il garde le même InChI).

Colonnes réécrites : INCHI, INCHIKEY, SMILES, FORMULA, EXACTMASS, AVERAGEMASS. NAME et ontologies
sont recopiés. Le résultat est écrit dans un NOUVEAU dossier : la base d'origine n'est pas modifiée.

Le traitement est long (plusieurs heures). Il affiche sa progression ; Ctrl+C l'arrête, et le relancer
reprend au fichier suivant (les fichiers déjà terminés sont gardés). L'ordinateur (Windows ou macOS) ne
se met pas en veille pendant le calcul.

Usage (depuis la racine du projet) :
  Windows : .venv\\Scripts\\python.exe scripts\\build_dataset\\renormalize_db.py [--processes N] [dossier_sortie]
  macOS   : .venv/bin/python scripts/build_dataset/renormalize_db.py [--processes N] [dossier_sortie]

--processes N : répartit les 16 fichiers entre N processus qui tournent en même temps (chacun avec
                nombre_de_coeurs / N threads). Par défaut 1. Sous Windows, plusieurs processus peuvent
                aller plus vite qu'un seul processus avec beaucoup de threads.
Par défaut, la sortie est scratch/renormalize_db/internal_databases.parquet ; une fois vérifiée, elle
remplace datas/internal_databases.parquet.
"""
import argparse
import os
import subprocess
import sys
import time

# Même réglage que scripts/FragHub.py : sous Windows, indiquer le dossier des DLL de RDKit C++.
if sys.platform == 'win32' and not getattr(sys, 'frozen', False):
    _rdkit_bin = os.path.join(os.environ.get('RDKIT_PREFIX') or os.path.join(os.path.expanduser('~'), 'rdkit_cpp'), 'Library', 'bin')
    if os.path.isdir(_rdkit_bin):
        os.add_dll_directory(_rdkit_bin)

import fraghub_rust


_caffeinate = None


def keep_awake(on):
    """Empêche la mise en veille tant que le calcul tourne : la veille dépend du clavier et de la souris,
    pas de l'activité du processeur. Windows : SetThreadExecutionState ; macOS : commande caffeinate."""
    global _caffeinate
    if sys.platform == 'win32':
        import ctypes
        ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
        ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | (ES_SYSTEM_REQUIRED if on else 0))
    elif sys.platform == 'darwin':
        if on and _caffeinate is None:
            # -i : pas de veille ; -w : s'arrête tout seul quand ce processus se termine
            _caffeinate = subprocess.Popen(["caffeinate", "-i", "-w", str(os.getpid())])
        elif not on and _caffeinate is not None:
            _caffeinate.terminate()
            _caffeinate = None


root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
src = os.path.join(root, "datas", "internal_databases.parquet")

parser = argparse.ArgumentParser(description="Re-normalise la base interne de FragHub.")
parser.add_argument("dossier_sortie", nargs="?", default=os.path.join(root, "scratch", "renormalize_db", "internal_databases.parquet"))
parser.add_argument("--processes", type=int, default=1, help="nombre de processus (défaut : 1)")
parser.add_argument("--worker", default=None, help=argparse.SUPPRESS)  # usage interne : fichiers d'un processus
args = parser.parse_args()
dst = args.dossier_sortie

# --- Processus de travail (lancé par le processus principal) : traite seulement ses fichiers ---
if args.worker is not None:
    files = [f for f in args.worker.split(",") if f]
    keep_awake(True)
    try:
        fraghub_rust.renormalize_internal_database(src, dst, files, os.environ.get("FRAGHUB_RENORM_PREFIX", "[P] "))
    except KeyboardInterrupt:
        sys.exit(1)
    finally:
        keep_awake(False)
    sys.exit(0)

# --- Processus principal ---
if os.path.abspath(src) == os.path.abspath(dst):
    sys.exit("ERREUR : le dossier de sortie doit être différent de la base d'origine.")
if os.path.exists(dst) and os.listdir(dst):
    print(f"{dst} existe déjà : les fichiers déjà terminés sont gardés, le traitement reprend au fichier suivant.")
    print("(Pour tout recommencer, supprimez ce dossier.)")

print(f"Base d'origine : {src}")
print(f"Sortie         : {dst}")
t = time.time()
keep_awake(True)
try:
    if args.processes <= 1:
        rows, with_structure, key_changed, smiles_changed, failed, duplicates = fraghub_rust.renormalize_internal_database(src, dst)
        print(f"\nTerminé en {(time.time() - t) / 60:.1f} min.")
        print(f"  {rows} lignes, dont {with_structure} avec une structure")
        print(f"  {key_changed} InChIKey modifiées, {smiles_changed} SMILES modifiés")
        print(f"  {failed} molécules illisibles par RDKit (valeurs d'origine gardées)")
        print(f"  {duplicates} InChIKey présentes sur plusieurs lignes (FragHub garde la dernière lue)")
    else:
        n = args.processes
        # Fichiers restant à faire, répartis par taille (le plus gros d'abord vers le processus le moins chargé)
        todo = [f for f in sorted(os.listdir(src)) if f.endswith(".parquet") and not os.path.exists(os.path.join(dst, f))]
        todo.sort(key=lambda f: -os.path.getsize(os.path.join(src, f)))
        buckets, loads = [[] for _ in range(n)], [0] * n
        for f in todo:
            k = loads.index(min(loads))
            buckets[k].append(f)
            loads[k] += os.path.getsize(os.path.join(src, f))
        threads = max(1, (os.cpu_count() or 1) // n)
        print(f"{n} processus de {threads} threads ; fichiers à faire : {len(todo)}")
        procs = []
        for k, files in enumerate(buckets):
            if not files:
                continue
            env = dict(os.environ, RAYON_NUM_THREADS=str(threads), FRAGHUB_RENORM_PREFIX=f"[P{k + 1}] ")
            print(f"  P{k + 1} : {', '.join(files)}")
            procs.append(subprocess.Popen([sys.executable, os.path.abspath(__file__), "--worker", ",".join(files), dst], env=env))
        try:
            codes = [p.wait() for p in procs]
        except KeyboardInterrupt:
            for p in procs:
                p.wait()
            sys.exit("\nArrêté (Ctrl+C). Relancer le script reprend au fichier suivant.")
        if any(codes):
            sys.exit("\nUn processus s'est arrêté avec une erreur (voir ci-dessus). Relancer le script reprend au fichier suivant.")
        # Tous les fichiers sont faits : ce dernier passage ne recalcule rien et compte les InChIKey en double
        duplicates = fraghub_rust.renormalize_internal_database(src, dst)[5]
        print(f"\nTerminé en {(time.time() - t) / 60:.1f} min (détail par fichier ci-dessus).")
        print(f"  {duplicates} InChIKey présentes sur plusieurs lignes (FragHub garde la dernière lue)")
except KeyboardInterrupt:
    sys.exit("\nArrêté (Ctrl+C). Relancer le script reprend au fichier suivant.")
finally:
    keep_awake(False)
