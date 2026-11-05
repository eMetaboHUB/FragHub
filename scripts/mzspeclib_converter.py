import io
import json
import multiprocessing
import os
import re
import sys
import traceback
import warnings
from mzspeclib.backends.msp import MSPSpectralLibrary
from mzspeclib.backends.json import JSONSpectralLibraryWriter
import time

# mzspeclib avertit pour chaque spectre dont l'énergie de collision n'est pas un nombre simple
# ("NOT FOUND", "-60,-20", "[15.0, 60.0]", "-70 eV"…). Le spectre est quand même converti (sans ce champ) :
# message sans intérêt, répété des milliers de fois. Seul ce message est masqué, les autres restent visibles.
# Placé ici (au chargement du module), le filtre s'applique aussi aux processus de conversion parallèle.
warnings.filterwarnings("ignore", message=r"Failed to parse .* for COLLISIONENERGY", category=UserWarning)

# ---------------------------------------------------------------------------
# Conversion MSP -> mzSpecLib JSON
#
# Version parallèle : chaque fichier MSP est découpé en morceaux (aux frontières de spectres),
# convertis en même temps sur tous les cœurs avec mzspeclib, puis recollés dans l'ordre.
# Le fichier obtenu est IDENTIQUE, octet pour octet, à celui de l'ancienne conversion séquentielle
# (convert_msp_to_mzspeclib_json, gardée plus bas comme référence et comme solution de repli) :
# - mêmes spectres, découpés avec les mêmes règles que l'index de mzspeclib (create_index) ;
# - chaque spectre est lu par le même code mzspeclib (_parse) avec le même numéro (donc la même
#   « library spectrum key ») et mis en forme par le même writer (write_spectrum) ;
# - le texte JSON est recollé exactement comme json.dump(liste, indent=2, sort_keys=True).
# En plus d'aller plus vite, on n'a plus tout le fichier en mémoire avant l'écriture.
# ---------------------------------------------------------------------------

# Règle de mzspeclib (backends/msp.py) : un spectre commence à une ligne "Name:", "NAME:",
# "Compound:" ou "COMPOUND:" (début de ligne, après suppression des espaces de fin).
_LEADER = re.compile(r"(Name|NAME|Compound|COMPOUND)\s*:")
_CHUNK_BYTES = 8 * 1024 * 1024  # taille visée d'un morceau (~ quelques milliers de spectres)


def _read_text_lines(path, start, end):
    """Lignes du fichier entre deux positions (en octets), lues comme mzspeclib les lit :
    UTF-8 et fins de ligne universelles (open_stream(..., encoding='utf8', newline=None))."""
    with open(path, "rb") as f:
        f.seek(start)
        data = f.read(end - start)
    return io.TextIOWrapper(io.BytesIO(data), encoding="utf8", newline=None)


def _spectrum_buffers(path, start, end, header_state):
    """Découpe un morceau en spectres, avec les mêmes règles que create_index + _buffer_from_stream
    de mzspeclib : lignes sans espaces de fin, lignes vides ignorées, lignes avant le premier
    spectre du fichier ignorées. Renvoie la liste des spectres (chacun = liste de lignes)."""
    buffers = []
    current = None
    state = "header" if header_state else "body"
    for line in _read_text_lines(path, start, end):
        line = line.rstrip()
        if state == "header":
            if _LEADER.match(line):
                state = "body"
            else:
                continue
        if len(line) == 0:
            continue
        if _LEADER.match(line):
            if current:
                buffers.append(current)
            current = []
        if current is None:  # ne se produit qu'au début d'un morceau sans en-tête
            current = []
        current.append(line)
    if current:
        buffers.append(current)
    return buffers


def _chunk_boundaries(path):
    """Positions (en octets) où couper le fichier : chacune est le début d'une ligne de début
    de spectre. Le premier morceau commence au début du fichier (lignes d'en-tête comprises)."""
    size = os.path.getsize(path)
    bounds = [0]
    with open(path, "rb") as f:
        target = _CHUNK_BYTES
        while target < size:
            f.seek(target)
            f.readline()  # aller au début de la ligne suivante
            pos = f.tell()
            found = None
            while pos < size:
                raw = f.readline()
                # Même lecture que mzspeclib (UTF-8, fins de ligne universelles) pour tester la ligne
                text = raw.decode("utf8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
                first = text.split("\n", 1)[0]
                if _LEADER.match(first.rstrip()):
                    found = pos
                    break
                pos = f.tell()
            if found is None or found <= bounds[-1]:
                break
            bounds.append(found)
            target = found + _CHUNK_BYTES
    bounds.append(size)
    return bounds


# --- Travail fait dans les processus parallèles --------------------------------

_worker_cache = {}


def _worker_tools(path):
    tools = _worker_cache.get(path)
    if tools is None:
        _worker_cache.clear()
        lib = MSPSpectralLibrary(path, create_index=False)  # pas d'index : on fournit nous-mêmes les spectres
        writer = JSONSpectralLibraryWriter(io.StringIO())
        tools = _worker_cache[path] = (lib, writer)
    return tools


def _count_task(args):
    path, k, start, end = args
    return path, k, len(_spectrum_buffers(path, start, end, header_state=(k == 0)))


def _convert_task(args):
    """Convertit un morceau. Écrit le texte JSON de ses spectres dans un fichier temporaire, déjà
    mis en forme comme dans json.dump(liste, indent=2, sort_keys=True) (niveau 1, séparateurs ",\\n  ")."""
    path, k, start, end, first_index, part_path = args
    lib, writer = _worker_tools(path)
    spectra = writer.buffer["spectra"]
    count = skipped = 0
    with open(part_path, "w", encoding="ascii", newline="\n") as out:
        for i, buf in enumerate(_spectrum_buffers(path, start, end, header_state=(k == 0))):
            number = first_index + i
            try:
                # Même lecture (_parse avec le même numéro) et même mise en forme que la version séquentielle.
                spec = lib._parse(buf, number)
                writer.write_spectrum(spec)
            except Exception as spec_e:
                skipped += 1
                print(f"[mzSpecLib] Spectre ignore (#{number}) dans {os.path.basename(path)}: {spec_e}")
                continue
            payload = spectra.pop()
            if count:
                out.write(",\n  ")
            out.write(json.dumps(payload, indent=2, sort_keys=True).replace("\n", "\n  "))
            count += 1
    return path, k, count, skipped


def _assemble(json_path, parts):
    """Écrit le fichier final à partir des morceaux (dans l'ordre), comme le faisait le writer :
    ouverture en mode texte 'wt' (même traduction des fins de ligne), liste JSON indentée,
    ou objet seul s'il n'y a qu'un spectre (option simplify de mzspeclib)."""
    total = sum(c for _, c in parts)
    with open(json_path, "wt") as handle:
        if total == 0:
            handle.write("[]")
        elif total == 1:
            text = next(open(p, "r", encoding="ascii", newline="\n").read() for p, c in parts if c)
            handle.write(text.replace("\n  ", "\n"))  # retour au niveau 0 (objet seul)
        else:
            handle.write("[\n  ")
            first = True
            for p, c in parts:
                if not c:
                    continue
                if not first:
                    handle.write(",\n  ")
                with open(p, "r", encoding="ascii", newline="\n") as f:
                    while True:
                        block = f.read(4 * 1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
                first = False
            handle.write("\n]")
    for p, _ in parts:
        try:
            os.remove(p)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Ancienne conversion séquentielle (référence exacte, et repli si le parallélisme échoue)
# ---------------------------------------------------------------------------

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
                # La lecture (parsing MSP) ET l'ecriture (JSON) peuvent echouer
                # independamment sur un spectre precis (ex: bug connu de mzspeclib
                # sur RT qui pose une cle "ERROR" invalide) : les deux etapes sont
                # couvertes ici pour qu'un seul spectre foireux ne fasse pas
                # avorter tout le fichier.
                spec = lib.get_spectrum(record.number)
                writer.write_spectrum(spec)
            except Exception as spec_e:
                skipped += 1
                print(f"[mzSpecLib] Spectre ignore (#{record.number}) dans {os.path.basename(input_msp)}: {spec_e}")
                if progress_cb:
                    progress_cb(1)
                continue

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


def _convert_all_sequential(jobs, prefix_cb, global_progress_cb):
    for msp_path, json_path in jobs:
        if prefix_cb:
            prefix_cb(f"Writing {os.path.basename(json_path)} to disk...")
        convert_msp_to_mzspeclib_json(msp_path, json_path, progress_cb=global_progress_cb)


def convert_all_msp(output_directory, step_cb=None, prefix_cb=None, progress_cb=None, total_items_cb=None, workers=None):
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

    jobs = []
    for msp_path in msp_files:
        rel_path = os.path.relpath(msp_path, msp_dir)
        json_path = os.path.join(json_dir, rel_path)
        json_path = json_path.replace('.msp', '.mzSpecLib.json')
        os.makedirs(os.path.dirname(json_path), exist_ok=True)
        jobs.append((msp_path, json_path))

    # Même réglage que l'ancien rdkit_worker.run_parallel (éprouvé dans l'exécutable PyInstaller)
    workers = workers or max(1, multiprocessing.cpu_count() - 1)
    try:
        pool = multiprocessing.Pool(processes=workers)
    except Exception as e:
        print(f"[mzSpecLib] Conversion parallèle impossible ({e}) : conversion séquentielle.")
        pool = None

    if pool is None:
        total_count = 0
        for msp_path, _ in jobs:
            with open(msp_path, 'r', encoding='utf-8', errors='ignore') as f:
                total_count += sum(1 for line in f if line.upper().startswith("NAME:"))
        if total_items_cb:
            total_items_cb(total_count)
        global_count = 0

        def global_progress_cb(count):
            nonlocal global_count
            global_count += 1
            if progress_cb and global_count % 100 == 0:
                progress_cb(global_count)
        _convert_all_sequential(jobs, prefix_cb, global_progress_cb)
        if progress_cb:
            progress_cb(global_count)
        return

    with pool:
        # 1. Découpage de chaque fichier, puis comptage des spectres de chaque morceau (en parallèle).
        #    Le nombre de spectres avant un morceau donne le numéro de son premier spectre.
        bounds = {p: _chunk_boundaries(p) for p, _ in jobs}
        count_tasks = [(p, k, b[k], b[k + 1]) for p, _ in jobs for b in [bounds[p]] for k in range(len(b) - 1)]
        counts = {}
        for p, k, n in pool.imap_unordered(_count_task, count_tasks):
            counts[(p, k)] = n
        n_chunks = {p: len(bounds[p]) - 1 for p, _ in jobs}
        # Cas particulier de mzspeclib : un fichier sans aucune ligne "Name:" compte quand même
        # un « spectre » (tout le fichier), lu à partir du début.
        totals = {p: sum(counts[(p, k)] for k in range(n_chunks[p])) for p, _ in jobs}

        if total_items_cb:
            total_items_cb(sum(totals.values()))

        done = 0
        for msp_path, json_path in jobs:
            if prefix_cb:
                prefix_cb(f"Writing {os.path.basename(json_path)} to disk...")
            if totals[msp_path] == 0:
                # Fichier sans début de spectre reconnu : comportement exact de l'ancienne conversion
                convert_msp_to_mzspeclib_json(msp_path, json_path)
                continue
            b = bounds[msp_path]
            tasks, first = [], 0
            for k in range(n_chunks[msp_path]):
                tasks.append((msp_path, k, b[k], b[k + 1], first, f"{json_path}.part{k:05d}"))
                first += counts[(msp_path, k)]
            results = {}
            skipped = 0
            try:
                for p, k, c, s in pool.imap_unordered(_convert_task, tasks):
                    results[k] = c
                    skipped += s
                    done += counts[(msp_path, k)]
                    if progress_cb:
                        progress_cb(done)  # peut lever "Process stopped by user." (bouton stop)
            except BaseException:
                # Arrêt ou erreur : on arrête les processus et on ne laisse pas de fichiers temporaires
                pool.terminate()
                pool.join()
                for t in tasks:
                    try:
                        os.remove(t[5])
                    except OSError:
                        pass
                raise
            _assemble(json_path, [(tasks[k][5], results[k]) for k in range(len(tasks))])
            if skipped:
                print(f"[mzSpecLib] {skipped} spectre(s) ignore(s) dans {os.path.basename(msp_path)} (non convertibles).")

    if progress_cb:
        progress_cb(done)
