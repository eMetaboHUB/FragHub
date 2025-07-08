import sys
import os
import traceback
import ctypes
import time
import platform
from threading import Thread

# Configuration de BASE_DIR
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

# Imports PyQt6
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QProgressBar, QTabWidget, QMessageBox, QSplashScreen,
    QSpacerItem, QSizePolicy, QTextEdit
)
from PyQt6.QtGui import QFont, QPixmap, QIcon, QPainter, QColor, QPen
from PyQt6.QtCore import (
    Qt, pyqtSignal, QSize, QTimer, QRectF, QObject, QThread
)

# Imports des composants UI (inchangés)
try:
    from scripts.GUI.tabs.tab_input import InputTab
    from scripts.GUI.tabs.tab_output import OutputTab
    from scripts.GUI.tabs.tab_filters import FiltersTab
    from scripts.GUI.tabs.tab_output_settings import OutputSettingTab
    from scripts.GUI.tabs.tab_projects import ProjectsTab
    from scripts.progress_window import ProgressWindow
except ImportError as e:
    print(f"Erreur critique d'importation des composants UI: {e}")
    # Dans une vraie application, afficher une QMessageBox ici avant de quitter
    sys.exit(f"Erreur critique d'importation: {e}")

def show_error_message(parent, title, message, on_close=None): # Ajouter le paramètre on_close=None
    """Affiche une erreur dans une QMessageBox élargie avec un QTextEdit."""
    msg_box = QMessageBox(parent)
    msg_box.setIcon(QMessageBox.Icon.Critical)
    msg_box.setWindowTitle(title)
    msg_box.setText("An error occurred during execution.")
    msg_box.setStandardButtons(QMessageBox.StandardButton.Ok)

    # Option 2: Ajouter un QTextEdit (plus de contrôle sur la taille/apparence)
    text_area = QTextEdit()
    text_area.setText(message)
    text_area.setReadOnly(True)
    text_area.setStyleSheet("font-family: Consolas, monospace; font-size: 10pt;")
    text_area.setMinimumSize(700, 350) # Ajuster la taille minimale

    try:
        grid_layout = msg_box.layout()
        grid_layout.addWidget(text_area, grid_layout.rowCount(), 0, 1, grid_layout.columnCount())
        msg_box.setMinimumSize(750, 450) # Ajuster la taille de la boîte de dialogue
    except Exception as layout_e:
        print(f"Warning: Could not add QTextEdit to QMessageBox layout ({layout_e}). Using setDetailedText as fallback.")
        msg_box.setDetailedText(message) # Solution de repli

    # Afficher la boîte de dialogue et attendre que l'utilisateur clique sur un bouton
    clicked_button = msg_box.exec()

    # --- NOUVEAU : Exécuter le callback APRÈS la fermeture de la boîte, si OK a été cliqué ---
    if clicked_button == QMessageBox.StandardButton.Ok and on_close:
        try:
            on_close() # Appeler la fonction passée en argument
        except Exception as callback_e:
            # Loguer une erreur si le callback lui-même échoue
            print(f"ERREUR: Exception dans le callback on_close de show_error_message: {callback_e}")
            traceback.print_exc()


class StreamCapturer:
    """
    Capture le flux de sortie (stdout / stderr) et le stocke dans un buffer.
    """

    def __init__(self, original_stream, signal=None):
        self.original_stream = original_stream  # stdout ou stderr d'origine
        self.signal = signal  # Signal optionnel pour transmettre les messages capturés
        self.buffer = []  # Stockage des messages capturés

    def write(self, message):
        # Stocker le message dans le buffer
        self.buffer.append(message)
        # Si un signal est défini, l'émettre
        if self.signal:
            self.signal.emit(message)
        # Envoyer également le message dans le flux d'origine pour debugger via terminal
        self.original_stream.write(message)

    def flush(self):
        # Nécessaire pour éviter les erreurs lors d'appels système via sys.stdout/sys.stderr
        self.original_stream.flush()

    def get_captured_text(self):
        # Retourne tout le contenu capturé
        return ''.join(self.buffer)

    def clear(self):
        # Réinitialise le buffer
        self.buffer = []


# --- Widget Spinner Personnalisé (Inchangé par rapport à votre version corrigée) ---
class InfiniteSpinnerWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._angle = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._update_angle)
        self._timer.setInterval(25)  # ~40 FPS

        self._spinner_color = QColor(Qt.GlobalColor.white)
        self._pen_width = 3.0

        self.setFixedSize(40, 40)

    def _update_angle(self):
        """Met à jour l'angle pour faire tourner l'animation dans le sens horaire."""
        self._angle = (self._angle - 10) % 360  # Sens horaire
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        padding = self._pen_width / 2.0
        rect_x = padding
        rect_y = padding
        rect_w = float(self.width()) - self._pen_width
        rect_h = float(self.height()) - self._pen_width

        if rect_w <= 0 or rect_h <= 0:
            return

        draw_rect = QRectF(rect_x, rect_y, rect_w, rect_h)
        pen = QPen(self._spinner_color, self._pen_width, Qt.PenStyle.SolidLine)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)

        start_angle = self._angle * 16
        span_angle = 90 * 16
        painter.drawArc(draw_rect, start_angle, span_angle)

    def start_animation(self):
        if not self._timer.isActive():
            self._timer.start()

    def stop_animation(self):
        if self._timer.isActive():
            self._timer.stop()

    def set_color(self, color: QColor):
        self._spinner_color = color
        self.update()

    def hideEvent(self, event):
        self.stop_animation()
        super().hideEvent(event)

    def showEvent(self, event):
        self.start_animation()
        super().showEvent(event)

# --- Fin du Widget Spinner Personnalisé ---


# --- Classe LoadingSplashScreen (Inchangée) ---
class LoadingSplashScreen(QWidget):
    #... (Identique à votre version)...
    def __init__(self, icon_pixmap: QPixmap, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.WindowType.SplashScreen | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)

        self.icon_label = QLabel()
        if not icon_pixmap.isNull():
            self.icon_label.setPixmap(icon_pixmap)
        else:
            self.icon_label.setText("Icon Error")
        self.icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.icon_label)

        layout.addSpacerItem(QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed))

        self.spinner = InfiniteSpinnerWidget(self)
        spinner_layout = QHBoxLayout()
        spinner_layout.addStretch()
        spinner_layout.addWidget(self.spinner)
        spinner_layout.addStretch()
        layout.addLayout(spinner_layout)

        layout.addSpacerItem(QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed))

        self.message_label = QLabel("Loading FragHub, please wait...")
        self.message_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.message_label.setStyleSheet("QLabel { color : white; font-weight: bold; }")
        layout.addWidget(self.message_label)

        self.setLayout(layout)
        self.adjustSize()

        if QApplication.primaryScreen():
            center_point = QApplication.primaryScreen().availableGeometry().center()
            self.move(center_point.x() - self.width() // 2, center_point.y() - self.height() // 2)

    def showMessage(self, message, font_size=28, alignment=Qt.AlignmentFlag.AlignCenter, color=Qt.GlobalColor.white):
        """Affiche un message dans le splash screen avec une taille de police configurable."""
        self.message_label.setText(message)
        self.message_label.setAlignment(alignment)
        self.message_label.setStyleSheet(f"""
            QLabel {{
                color: white;
                font-weight: bold;
                font-size: {font_size}px;  /* Appliquer la taille spécifiée */
            }}
        """)

    def closeEvent(self, event):
        self.spinner.stop_animation()
        super().closeEvent(event)
# --- Fin de LoadingSplashScreen ---


# --- Définition de AppUserModelID (Inchangé) ---
if platform.system() == "Windows":
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")


# --- Classe MainWindow (Inchangée) ---
class MainWindow(QMainWindow):
    #... (Identique à votre version)...
    update_progress_signal = pyqtSignal(int)

    # Signal émis lorsqu'une erreur se produit dans le thread worker
    error_occurred_signal = pyqtSignal(str)
    # Signal émis lorsque la tâche (MAIN) est terminée (succès ou erreur)
    task_finished_signal = pyqtSignal()

    def __init__(self, main_function_ref):
        super().__init__()
        self.setWindowFlags(
            Qt.WindowType.Window | Qt.WindowType.WindowMinimizeButtonHint | Qt.WindowType.WindowCloseButtonHint
        )
        self.MAIN_function = main_function_ref
        self.setWindowTitle("FragHub 1.3.2")
        self.setGeometry(100, 100, 1280, 720)

        main_layout = QVBoxLayout()
        banner = QLabel()
        # Utilisation de os.path.join pour la compatibilité
        icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
        pixmap = QPixmap(icon_path)
        if not pixmap.isNull():
             pixmap_scaled = pixmap.scaled(
                 200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
             )
             banner.setPixmap(pixmap_scaled)
        else:
             print(f"Warning: Could not load banner image at {icon_path}")
             banner.setText("FragHub Logo") # Fallback text
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)

        self.tabs = QTabWidget()
        self.input_tab = InputTab() # Garder une référence si nécessaire
        self.output_tab = OutputTab()
        self.filters_tab = FiltersTab() # Garder une référence si nécessaire
        self.output_settings_tab = OutputSettingTab() # Garder une référence si nécessaire
        self.projects_tab = ProjectsTab()

        self.tabs.addTab(self.input_tab, "INPUT")
        self.tabs.addTab(self.output_tab, "OUTPUT")
        self.tabs.addTab(self.filters_tab, "Filters settings")
        self.tabs.addTab(self.output_settings_tab, "Output settings")
        self.tabs.addTab(self.projects_tab, "Projects settings")
        main_layout.addWidget(self.tabs)

        # Connexion du signal pour le répertoire de sortie
        self.output_tab.output_directory_changed.connect(self.projects_tab.output_directory_changed_signal)

        self.start_button = QPushButton("START")
        self.start_button.setFixedSize(140, 60)
        self.start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        self.start_button.clicked.connect(self.open_progress_window)
        main_layout.addWidget(self.start_button, alignment=Qt.AlignmentFlag.AlignCenter)

        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        # --- Variables d'état ---
        self.running = False
        self.progress_window = None
        self.thread = None  # Référence au thread worker (type: threading.Thread)
        self.stop_thread_flag = False  # Flag pour demander l'arrêt du thread

        # --- Connexion des signaux aux slots ---
        self.error_occurred_signal.connect(self.handle_execution_error)


    def open_progress_window(self):
        if not self.running:
            # Cacher la fenêtre principale pendant l'exécution
            self.showMinimized()
            self.setEnabled(False)
            # Créer et montrer la fenêtre de progression
            # Assurez-vous que ProgressWindow est bien un QWidget ou QDialog
            self.progress_window = ProgressWindow(parent=self) # parent=None pour la rendre indépendante
            self.progress_window.setWindowFlags(
                Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.Tool
            )
            self.progress_window.show()
            # Démarrer l'exécution dans un thread séparé
            self.start_execution()

    def start_execution(self):
        """Démarre le traitement lorsque l'utilisateur appuie sur START."""
        if self.running:
            print("Un traitement est déjà en cours.")
            return

        # Réinitialisation d'état pour un nouveau traitement
        self.running = True
        self.stop_thread_flag = False  # Réinitialiser le flag d'arrêt au cas où
        if self.progress_window:
            self.progress_window.progress_bar_widget.update_total_items(total=100, completed=0)
            self.progress_window.progress_bar_widget.update_progress_bar(0)

        # Stopper un thread précédent s'il était encore actif
        if self.thread and self.thread.is_alive():
            self.stop_thread_flag = True
            self.thread.join()  # Attendre que le thread se termine proprement
            self.thread = None

        # Lancer un nouveau thread
        self.thread = Thread(target=self.run_main_function, daemon=True, name="MainWorkerThread_FragHub")
        self.thread.start()

    def _cleanup_after_error_or_finish(self):
        """
        Nettoie l'état après une erreur (ou pourrait être utilisé pour une fin normale si nécessaire).
        Ferme la fenêtre de progression, réinitialise l'état et montre la fenêtre principale.
        """
        print("Exécution du nettoyage (_cleanup_after_error_or_finish)...")
        # Fermer la fenêtre de progression si elle existe
        if self.progress_window:
            try:
                print("Fermeture de progress_window depuis le nettoyage.")
                self.progress_window.close()
                self.progress_window = None
            except Exception as e:
                print(f"Erreur lors de la fermeture de progress_window pendant le nettoyage: {e}")

        # Réinitialiser l'état (comme dans handle_task_finished)
        self.running = False
        self.thread = None
        self.stop_thread_flag = False
        print("Affichage de la fenêtre principale après nettoyage.")
        self.show()
        self.activateWindow()


    def run_main_function(self):
        """
        Wrapper exécuté dans le thread worker. Appelle MAIN et gère les signaux.
        NE PAS interagir directement avec l'UI ici (sauf via signaux émis).
        """
        try:
            if not self.progress_window:
                print("Warning: Progress window non initialisée avant l'appel à MAIN.")
                # On pourrait lever une exception ou continuer prudemment

            # Appel de la fonction principale avec les callbacks et le stop_flag
            # Assurez-vous que votre fonction MAIN accepte 'stop_flag'
            self.MAIN_function(
                progress_callback=self.progress_window.update_progress_signal.emit,
                total_items_callback=self.progress_window.update_total_signal.emit,
                prefix_callback=self.progress_window.update_prefix_signal.emit,
                item_type_callback=self.progress_window.update_item_type_signal.emit,
                step_callback=self.progress_window.update_step_signal.emit,
                completion_callback=self.progress_window.completion_callback.emit,
                deletion_callback=self.progress_window.deletion_callback.emit,
                # Passer une fonction lambda qui vérifie le flag
                stop_flag=lambda: self.stop_thread_flag
            )
        except InterruptedError as ie:
            print(f"Execution interrompue par l'utilisateur: {ie}")
            # Émettre un signal d'erreur spécifique ou juste terminer proprement
            # Dans ce cas, on ne considère pas ça comme une "erreur" à afficher à l'utilisateur
            # Le finally s'exécutera pour nettoyer.
        except Exception:  # Capturer TOUTES les autres exceptions
            full_traceback = traceback.format_exc()
            print(f"Erreur pendant l'exécution de MAIN dans le thread:\n{full_traceback}")
            # Émettre le signal d'erreur vers le thread GUI
            self.error_occurred_signal.emit(full_traceback)
        finally:
            # Ce bloc s'exécute toujours (succès, erreur, interruption)
            # Émettre le signal de fin pour que le thread GUI nettoie
            self.task_finished_signal.emit()

    def handle_execution_error(self, traceback_str):
        """
        Slot exécuté dans le thread GUI pour afficher les erreurs venant du worker.
        Après clic sur OK dans le message d'erreur, l'application sera fermée proprement.
        """
        print("handle_execution_error appelé.")
        # NE PAS fermer la fenêtre de progression ici immédiatement.

        # Afficher la boîte de dialogue d'erreur.
        # Passer self.clean_exit comme callback à exécuter après le clic sur OK.
        show_error_message(
            parent=self,  # Le parent peut être self (MainWindow), même si elle est cachée
            title="Execution Error",
            message=traceback_str,
            on_close=self.clean_exit  # Passer la référence à la méthode clean_exit
            # Par défaut, clean_exit(force_quit_app=True) quittera l'application.
        )

    # --- Méthodes pour l'arrêt propre ---

    def clean_exit(self, force_quit_app=True):
        """
        Tente d'arrêter le thread worker proprement et quitte l'application.
        :param force_quit_app: Si True, quitte QApplication à la fin.
        """
        if self.running and self.thread and self.thread.is_alive():
            print("Demande d'arrêt au thread worker...")
            self.stop_thread_flag = True  # Signaler au thread de s'arrêter

            # Attendre un peu que le thread se termine via le stop_flag
            wait_start_time = time.time()
            timeout_seconds = 3.0  # Attendre max 3 secondes
            while self.thread.is_alive() and (time.time() - wait_start_time) < timeout_seconds:
                QApplication.processEvents()  # Permettre à la boucle d'événements de tourner
                time.sleep(0.05)  # Petite pause

            if self.thread.is_alive():
                print(f"WARNING: Thread worker '{self.thread.name}' n'a pas terminé après {timeout_seconds}s.")
                # Ici, on pourrait décider de forcer ou juste continuer la fermeture de l'UI

        # Fermer la fenêtre de progression si elle existe encore
        if self.progress_window:
            try:
                self.progress_window.close()
                self.progress_window = None
            except Exception as e:
                print(f"Erreur lors de la fermeture de progress_window dans clean_exit: {e}")

        self.running = False  # S'assurer que l'état est correct

        if force_quit_app:
            QApplication.instance().quit()
        else:
            # Si on ne quitte pas l'app (ex: juste annuler la tâche),
            # s'assurer que la fenêtre principale est visible
            self.show()
            self.activateWindow()

    def closeEvent(self, event):
        """Gère la fermeture de la fenêtre principale par l'utilisateur (clic sur X)."""
        if self.running:
            reply = QMessageBox.question(self, 'Quit FragHub ?',
                                         "Are you sure you want to quit?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                         QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes:
                self.clean_exit(force_quit_app=True)  # Tenter l'arrêt propre et quitter
                event.accept()  # Accepter la fermeture de la fenêtre

        else:
            # Pas besoin d'appeler clean_exit ici si on quitte juste l'appli
            # QApplication.quit() sera appelé par la fermeture normale
            event.accept()  # Accepter la fermeture


# --- Fin de MainWindow ---


# --- Worker pour les tâches de démarrage (NOUVEAU) ---
class StartupWorker(QObject):
    """Effectue les tâches de démarrage longues dans un thread séparé."""
    # Signal émis quand les tâches sont finies, passe la fonction MAIN importée
    finished = pyqtSignal(object)
    # Signal émis en cas d'erreur
    error = pyqtSignal(str)
    # Signal pour mettre à jour le message du splash screen
    update_splash_message = pyqtSignal(str, int)

    def __init__(self, base_dir):
        super().__init__()
        self._base_dir = base_dir

    def showMessage(self, message, font_size=28, alignment=Qt.AlignmentFlag.AlignCenter, color=Qt.GlobalColor.white):
        """Affiche un message dans le splash screen avec taille de police configurable."""
        self.message_label.setText(message)
        self.message_label.setAlignment(alignment)
        self.message_label.setStyleSheet(f"""
            QLabel {{
                color: white;
                font-weight: bold;
                font-size: {font_size}px;  /* Ajustement dynamique de la taille de la police */
            }}
        """)

    def run_startup_tasks(self):
        """Exécute l'importation et autres tâches longues."""
        try:
            # --- Tâche potentiellement longue 1 : Importation ---
            self.update_splash_message.emit("Loading FragHub, please wait...", 20)  # Message 1
            time.sleep(1)  # Simule un léger retard pour montrer le spinner

            # Importer le composant principal
            from scripts.MAIN import MAIN as imported_main

            # --- Tâche potentiellement longue 2 : Préparation ---
            self.update_splash_message.emit("Initializing main window", 20)  # Message 2
            time.sleep(1)

            # Émettre le signal de fin en cas de succès
            self.finished.emit(imported_main)

        except ImportError as e:
            error_msg = f"Failed to import 'scripts.MAIN': {e}\n{traceback.format_exc()}"
            print(f"ERROR in worker thread: {error_msg}")
            self.error.emit(error_msg)  # Émettre le signal d'erreur
        except Exception as e:
            error_msg = f"Unexpected error during startup tasks: {e}\n{traceback.format_exc()}"
            print(f"ERROR in worker thread: {error_msg}")
            self.error.emit(error_msg)  # Émettre le signal d'erreur


# --- Fin de StartupWorker ---
# --- Exception Hook Global ---
def exception_hook(exctype, value, tb):
    error_message = ''.join(traceback.format_exception(exctype, value, tb))
    sys.stderr.write(f"Uncaught exception:\n{error_message}")
    # Tenter d'afficher une QMessageBox même pour les erreurs non gérées
    try:
        # Vérifier si QApplication existe avant d'essayer de créer un widget
        if QApplication.instance():
             message_box = QMessageBox()
             message_box.setIcon(QMessageBox.Icon.Critical)
             message_box.setWindowTitle("Unhandled Exception")
             message_box.setText("An critical unexpected error occurred!")
             message_box.setDetailedText(error_message)
             # Définir une taille minimale pour la boîte de dialogue
             message_box.setMinimumSize(600, 300)
             # Tenter d'ajouter le texte détaillé dans une zone scrollable (si setDetailedText ne suffit pas)
             try:
                  text_edit = message_box.findChild(QTextEdit)
                  if text_edit:
                      text_edit.setMinimumSize(550, 200)
             except Exception:
                  pass # Ignorer si on ne peut pas trouver/modifier le QTextEdit
             message_box.exec()
        else:
             print("QApplication non disponible, impossible d'afficher QMessageBox pour l'erreur non gérée.")
    except Exception as e:
        print(f"Impossible d'afficher QMessageBox pour l'erreur non gérée: {e}")
    sys.exit(1)



# --- Fonction run_GUI (MODIFIÉE pour utiliser le threading) ---
def run_GUI():
    # Définir l'exception hook global *avant* de créer QApplication
    sys.excepthook = exception_hook

    app_args = sys.argv if hasattr(sys, 'argv') else ['']
    app = QApplication(app_args)

    # Configuration de l'icône de l'application (inchangé)
    if platform.system() == "Darwin":
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.icns")
    else:
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_Python_icon.ico")
    if not os.path.exists(app_icon_path):
        app_icon_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    if os.path.exists(app_icon_path):
        app.setWindowIcon(QIcon(app_icon_path))
    else:
        print(f"Warning: Application icon not found at {app_icon_path} or fallback path.")

    # --- Splash Screen Setup ---
    splash_pix_path = os.path.join(BASE_DIR, "GUI", "assets", "FragHub_icon.png")
    splash_pixmap = QPixmap(splash_pix_path)
    if splash_pixmap.isNull(): splash_pixmap = QPixmap(1, 1); splash_pixmap.fill(Qt.GlobalColor.transparent)

    splash = LoadingSplashScreen(splash_pixmap) # Utiliser la classe placeholder si besoin
    splash.showMessage("Loading FragHub...") # Message initial
    splash.show()

    # --- Threading Setup for Startup Tasks ---
    startup_thread = QThread()
    startup_worker = StartupWorker(BASE_DIR) # Utiliser la classe placeholder si besoin
    startup_worker.moveToThread(startup_thread)
    shared_state = {'main_window': None, 'main_function': None, 'splash_screen': splash}

    # --- Slots Locaux ---
    def on_startup_complete(imported_main_function):
        if not QApplication.instance(): return # Sécurité si l'app a été quittée entre temps
        shared_state['main_function'] = imported_main_function
        try:
            # Création de la MainWindow après succès du worker
            shared_state['main_window'] = MainWindow(main_function_ref=shared_state['main_function'])
            shared_state['main_window'].show()
        except Exception as e:
             print(f"FATAL ERROR during MainWindow creation: {e}")
             on_startup_error(f"Failed to create main window: {e}\n{traceback.format_exc()}")
             return # Important: Ne pas fermer le splash si l'erreur vient ici

        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None

    def on_startup_error(error_message):
        if not QApplication.instance(): return
        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None
        QMessageBox.critical(None, "Fatal Startup Error", f"Could not start FragHub:\n{error_message}")
        app.quit()

    def update_splash(message, font_size=12): # Taille par défaut ajustée
        if shared_state['splash_screen']:
            # Appeler la méthode showMessage de l'objet splash réel
            try:
                 # On suppose que LoadingSplashScreen a une méthode showMessage(msg)
                 shared_state['splash_screen'].showMessage(message) # Ajuster si la méthode a des params différents
            except Exception as e:
                 print(f"Could not update splash message: {e}")


    # --- Connect Signals and Slots ---
    startup_thread.started.connect(startup_worker.run_startup_tasks)
    startup_worker.finished.connect(on_startup_complete)
    startup_worker.error.connect(on_startup_error)
    startup_worker.update_splash_message.connect(update_splash) # Assurez-vous que StartupWorker émet bien ce signal

    # Nettoyage (inchangé)
    startup_worker.finished.connect(startup_thread.quit)
    startup_worker.error.connect(startup_thread.quit) # Quitter aussi en cas d'erreur
    startup_worker.finished.connect(startup_worker.deleteLater)
    startup_worker.error.connect(startup_worker.deleteLater)
    startup_thread.finished.connect(startup_thread.deleteLater)

    # --- Démarrer le thread ---
    startup_thread.start()

    # --- Boucle principale ---
    exit_code = app.exec()
    sys.exit(exit_code)

# --- Bloc d'exécution principal (Inchangé) ---
if __name__ == "__main__":
    # Vous pouvez ajouter des logs ou initialisations ici si nécessaire
    run_GUI()
