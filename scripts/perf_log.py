"""Journal de performances de FragHub (temporaire, pour comparer les optimisations).

À chaque étape du traitement (les mêmes que les barres de progression), on note :
le nom de l'étape et de la tâche, le nombre d'éléments traités, la durée, la mémoire utilisée
(à la fin et au maximum pendant la tâche) et le nombre moyen de cœurs CPU occupés.

À la fin de chaque run, trois fichiers sont écrits dans le sous-dossier PERF_LOGS du dossier de sortie
(ou dans le dossier indiqué par la variable d'environnement FRAGHUB_PERF_DIR) :
- perf_<date>.csv : une ligne par tâche (séparateur « ; », s'ouvre directement dans Excel) ;
- perf_<date>.txt : le même tableau lisible, avec les fichiers d'entrée, les paramètres et les
  messages de suppression ;
- historique.csv  : une ligne par run, pour comparer les runs entre eux.
Attention : l'option reset_updates vide le dossier de sortie au début du run, donc aussi les journaux
des runs précédents (celui du run en cours est écrit à la fin, il est conservé).

Pour nommer un run (ex. « avant-optim-msp »), définir la variable d'environnement FRAGHUB_PERF_LABEL.

Pour tout retirer plus tard : supprimer ce fichier et les lignes marquées « # PERF_LOG » dans FragHub.py.
"""
import os
import threading
import time
from datetime import datetime

try:
    import psutil
    _PROC = psutil.Process()
except Exception:  # psutil absent : on garde les durées, sans la mémoire ni le CPU
    psutil = None
    _PROC = None

_SAMPLE_PERIOD = 0.5  # secondes entre deux mesures de mémoire
_MIN_GAP = 0.2  # un temps mort entre deux tâches plus long que ça apparaît dans le journal
_lock = threading.Lock()
_state = None  # run en cours (None = aucun)
_init = None  # dernier chargement des bases internes (fait au démarrage, avant le run)


def _rss_mb():
    if _PROC is None:
        return None
    try:
        return _PROC.memory_info().rss / (1024 * 1024)
    except Exception:
        return None


def _peak_mb():
    """Pic de mémoire du processus depuis son lancement (Windows : peak_wset ; ailleurs : None)."""
    if _PROC is None:
        return None
    try:
        peak = getattr(_PROC.memory_info(), "peak_wset", None)
        return peak / (1024 * 1024) if peak else None
    except Exception:
        return None


def _cpu_s():
    if _PROC is None:
        return None
    try:
        t = _PROC.cpu_times()
        return t.user + t.system
    except Exception:
        return None


class _Task:
    def __init__(self, step, total, prefix, item_type):
        self.step = step
        self.prefix = prefix
        self.item_type = item_type
        self.total = total
        self.last_value = None
        self.progress_events = 0
        self.t0 = time.perf_counter()
        self.cpu0 = _cpu_s()
        self.rss_max = _rss_mb()
        self.t1 = None
        self.rss_end = None
        self.cpu1 = None

    def sample(self):
        rss = _rss_mb()
        if rss is not None and (self.rss_max is None or rss > self.rss_max):
            self.rss_max = rss

    def close(self):
        self.sample()
        self.t1 = time.perf_counter()
        self.rss_end = _rss_mb()
        self.cpu1 = _cpu_s()

    @property
    def duration(self):
        return (self.t1 or time.perf_counter()) - self.t0

    @property
    def cores(self):
        if self.cpu0 is None or self.cpu1 is None or self.duration <= 0:
            return None
        return (self.cpu1 - self.cpu0) / self.duration


class _Run:
    def __init__(self, params):
        self.params = dict(params or {})
        self.started_at = datetime.now()
        self.t0 = time.perf_counter()
        self.rss_start = _rss_mb()
        self.step = ""
        self.tasks = []
        self.current = None
        self.last_end = self.t0  # fin de la dernière tâche (pour repérer le temps passé entre deux tâches)
        self.deletions = []
        self.status = None  # dernier message de fin reçu (terminé, arrêté, erreur…)
        self.stop = threading.Event()
        self.sampler = threading.Thread(target=self._sample_loop, daemon=True)
        self.sampler.start()

    def _sample_loop(self):
        # Mesure la mémoire régulièrement, pour attraper les pics entre deux messages de progression
        while not self.stop.wait(_SAMPLE_PERIOD):
            with _lock:
                if self.current is not None:
                    self.current.sample()

    def close_current(self):
        if self.current is not None:
            self.current.close()
            self.last_end = self.current.t1
            self.current = None

    def new_task(self, total=None, prefix=None, item_type=None):
        self.close_current()
        # Temps passé hors de toute tâche (copies, libérations de mémoire…) : on le note aussi
        gap = time.perf_counter() - self.last_end
        if gap >= _MIN_GAP:
            g = _Task(self.step, None, "(entre deux tâches)", None)
            g.t0 = self.last_end
            g.cpu0 = None
            g.close()
            self.tasks.append(g)
        self.current = _Task(self.step, total, prefix, item_type)
        self.tasks.append(self.current)


# ---------------------------------------------------------------------------
# Fonctions appelées depuis FragHub.py (chacune ne fait rien si aucun run n'est en cours)
# ---------------------------------------------------------------------------

def record_init(duration_s):
    """Chargement des bases internes (au démarrage de l'application, avant le run)."""
    global _init
    _init = (duration_s, _rss_mb())


def start(params):
    global _state
    with _lock:
        if _state is not None:
            _state.stop.set()
        _state = _Run(params)


def on_step(name):
    with _lock:
        if _state is None:
            return
        _state.close_current()
        _state.step = str(name)


# Rust annonce chaque tâche par trois messages (nombre d'éléments, libellé, unité), dans un ordre
# qui varie selon les modules, et parfois par un libellé seul (ex. "Writing deletion logs...").
# Règle : un message complète la tâche en cours tant qu'elle n'a reçu aucune progression et que
# l'information manque encore ; sinon il ouvre une nouvelle tâche.

def _fresh(task):
    return task is not None and task.progress_events == 0


def on_total(total):
    with _lock:
        if _state is None:
            return
        c = _state.current
        if _fresh(c) and c.total is None:
            c.total = total
        else:
            _state.new_task(total=total)


def on_prefix(prefix):
    with _lock:
        if _state is None:
            return
        c = _state.current
        if _fresh(c) and c.prefix is None:
            c.prefix = str(prefix)
        else:
            _state.new_task(prefix=str(prefix))


def on_item_type(item_type):
    with _lock:
        if _state is None:
            return
        c = _state.current
        if _fresh(c) and c.item_type is None:
            c.item_type = str(item_type)
        else:
            _state.new_task(item_type=str(item_type))


def on_progress(value):
    with _lock:
        if _state is None or _state.current is None:
            return
        _state.current.progress_events += 1
        _state.current.last_value = value


def on_deletion(message):
    with _lock:
        if _state is None:
            return
        _state.deletions.append(str(message).strip())


def on_completion(message):
    with _lock:
        if _state is None:
            return
        _state.status = str(message).strip()


def _fmt(v, digits=1):
    if v is None:
        return ""
    return f"{v:.{digits}f}"


def _input_files(params):
    files = []
    for p in params.get("input_directory") or []:
        try:
            size = os.path.getsize(p) / (1024 * 1024)
        except OSError:
            size = None
        files.append((p, size))
    return files


def finish(status=None):
    """Écrit les fichiers du journal. Renvoie le chemin du fichier détaillé (ou None)."""
    global _state
    with _lock:
        s = _state
        _state = None
        if s is None:
            return None
        s.stop.set()
        s.close_current()
    try:
        return _write(s, status or s.status or "terminé (sans message de fin)")
    except Exception as e:  # le journal ne doit jamais faire échouer FragHub
        print(f"[PERF_LOG] Impossible d'écrire le journal : {e}")
        return None


def _write(s, status):
    # Dossier du journal : <dossier de sortie>/PERF_LOGS (ou FRAGHUB_PERF_DIR s'il est défini,
    # ou ~/FragHub_perf_logs si le dossier de sortie est inconnu ou impossible à créer)
    out_dir = os.environ.get("FRAGHUB_PERF_DIR")
    if not out_dir and s.params.get("output_directory"):
        out_dir = os.path.join(s.params["output_directory"], "PERF_LOGS")
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError:
            out_dir = None
    if not out_dir:
        out_dir = os.path.join(os.path.expanduser("~"), "FragHub_perf_logs")
    os.makedirs(out_dir, exist_ok=True)
    label = os.environ.get("FRAGHUB_PERF_LABEL", "").strip()
    stamp = s.started_at.strftime("%Y%m%d-%H%M%S")
    base = f"perf_{stamp}" + (f"_{label}" if label else "")
    total_s = time.perf_counter() - s.t0
    rss_max = max([t.rss_max for t in s.tasks if t.rss_max is not None] or [0]) or None
    peak = _peak_mb()
    files = _input_files(s.params)
    input_mb = sum(sz for _, sz in files if sz)

    # --- CSV détaillé : une ligne par tâche ---
    header = ["n", "etape", "tache", "nb_elements", "unite", "traites", "debut_s", "duree_s",
              "rss_fin_mo", "rss_max_mo", "coeurs_cpu_moyen"]
    rows = []
    if _init is not None:
        rows.append(["0", "Démarrage", "Chargement des bases internes (au lancement de l'app)", "", "", "",
                     "", _fmt(_init[0], 2), _fmt(_init[1], 0), "", ""])
    for i, t in enumerate(s.tasks, 1):
        rows.append([str(i), t.step, (t.prefix or "").rstrip(":").strip(), "" if t.total is None else str(t.total), t.item_type or "",
                     "" if t.last_value is None else str(t.last_value),
                     _fmt(t.t0 - s.t0, 2), _fmt(t.duration, 2), _fmt(t.rss_end, 0), _fmt(t.rss_max, 0),
                     _fmt(t.cores, 1)])

    def csv_line(values):
        return ";".join('"' + v.replace('"', '""') + '"' if (";" in v or '"' in v) else v for v in values)

    def fr(row, numeric_cols):
        # Virgule décimale dans les colonnes numériques : Excel en français lit alors de vrais nombres
        return [v.replace(".", ",") if i in numeric_cols else v for i, v in enumerate(row)]

    csv_path = os.path.join(out_dir, base + ".csv")
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:  # utf-8-sig : accents corrects dans Excel
        f.write("\r\n".join(csv_line(r) for r in [header] + [fr(r, {6, 7, 8, 9, 10}) for r in rows]) + "\r\n")

    # --- Texte lisible ---
    lines = [
        f"FragHub — journal de performances",
        f"Date      : {s.started_at:%Y-%m-%d %H:%M:%S}" + (f"   (run « {label} »)" if label else ""),
        f"Statut    : {status}",
        f"Durée     : {total_s:.1f} s",
        f"Mémoire   : max pendant le run {_fmt(rss_max, 0)} Mo"
        + (f", pic du processus depuis son lancement {_fmt(peak, 0)} Mo" if peak else ""),
        f"Entrée    : {len(files)} fichier(s), {input_mb:.0f} Mo",
    ]
    for p, sz in files:
        lines.append(f"            {os.path.basename(p)} ({_fmt(sz, 0)} Mo)")
    lines.append(f"Sortie    : {s.params.get('output_directory', '')}")
    params = {k: v for k, v in s.params.items() if k not in ("input_directory", "input_db_names", "output_directory")}
    lines.append("Paramètres: " + ", ".join(f"{k}={v}" for k, v in sorted(params.items())))
    lines.append("")
    widths = [max(len(r[c]) for r in [header] + rows) for c in range(len(header))]
    widths[2] = min(widths[2], 60)
    for r in [header] + rows:
        cells = [r[c][:widths[c]].ljust(widths[c]) if c in (1, 2, 4) else r[c].rjust(widths[c]) for c in range(len(r))]
        lines.append("  ".join(cells))
    if s.deletions:
        lines.append("")
        lines.append("Messages de suppression :")
        lines.extend("  " + d for d in s.deletions)
    txt_path = os.path.join(out_dir, base + ".txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    # --- Historique : une ligne par run ---
    hist_path = os.path.join(out_dir, "historique.csv")
    hist_header = ["date", "label", "statut", "duree_totale_s", "rss_max_mo", "nb_fichiers", "entree_mo", "detail"]
    new = not os.path.exists(hist_path)
    with open(hist_path, "a", encoding="utf-8-sig" if new else "utf-8", newline="") as f:
        if new:
            f.write(csv_line(hist_header) + "\r\n")
        f.write(csv_line(fr([f"{s.started_at:%Y-%m-%d %H:%M:%S}", label, status, _fmt(total_s, 1), _fmt(rss_max, 0),
                             str(len(files)), _fmt(input_mb, 0), base + ".csv"], {3, 4, 6})) + "\r\n")

    print(f"[PERF_LOG] Journal écrit : {txt_path}")
    return txt_path
