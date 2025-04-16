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
    QSpacerItem, QSizePolicy
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
try:
    if platform.system() == "Windows":
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")
except AttributeError:
     print("Note: Could not set AppUserModelID (not on Windows or ctypes issue).")
     pass


# --- Classe MainWindow (Inchangée) ---
class MainWindow(QMainWindow):
    #... (Identique à votre version)...
    update_progress_signal = pyqtSignal(int)

    def __init__(self, main_function_ref):
        super().__init__()
        self.MAIN_function = main_function_ref
        self.setWindowTitle("FragHub 1.3.0")
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

        self.running = False
        self.progress_window = None
        self.thread = None # Renommé pour éviter conflit avec QThread plus bas

    def open_progress_window(self):
        if not self.running:
            # Cacher la fenêtre principale pendant l'exécution
            self.hide()
            # Créer et montrer la fenêtre de progression
            # Assurez-vous que ProgressWindow est bien un QWidget ou QDialog
            self.progress_window = ProgressWindow(parent=None) # parent=None pour la rendre indépendante
            self.progress_window.show()
            # Démarrer l'exécution dans un thread séparé
            self.start_execution()

    def start_execution(self):
        if self.MAIN_function is None:
            print("Erreur critique: La référence à la fonction MAIN n'a pas été passée à MainWindow.")
            QMessageBox.critical(self, "Erreur de Configuration", "La fonction principale n'a pas pu être chargée.")
            self.show() # Remontrer la fenêtre principale en cas d'erreur
            if self.progress_window: self.progress_window.close()
            self.progress_window = None
            return

        self.running = True
        # Utilisation de threading.Thread pour la fonction MAIN
        self.thread = Thread(target=self.run_main_function, daemon=True)
        self.thread.start()

    def run_main_function(self):
        try:
            # Appel de la fonction principale avec les callbacks
            self.MAIN_function(
                progress_callback=self.progress_window.update_progress_signal.emit,
                total_items_callback=self.progress_window.update_total_signal.emit,
                prefix_callback=self.progress_window.update_prefix_signal.emit,
                item_type_callback=self.progress_window.update_item_type_signal.emit,
                step_callback=self.progress_window.update_step_signal.emit,
                completion_callback=self.progress_window.completion_callback.emit,
                deletion_callback=self.progress_window.deletion_callback.emit,
            )
        except Exception as e:
            print(f"Erreur pendant l'exécution de MAIN dans le thread: {e}")
            traceback.print_exc()
            # Optionnel: Afficher une erreur à l'utilisateur via un signal vers le thread principal
            # self.error_signal.emit(f"Une erreur est survenue: {e}")
        finally:
            # Assurer la fermeture propre et le retour à la fenêtre principale
            # Utiliser QMetaObject.invokeMethod pour interagir avec l'UI depuis ce thread
            if self.progress_window:
                 # Fermer la fenêtre de progression depuis le thread principal
                 # Ceci nécessite un signal ou QMetaObject.invokeMethod
                 # Pour simplifier, on suppose que completion_callback gère la fermeture
                 # ou on peut utiliser un signal dédié pour fermer la fenêtre
                 pass # La fermeture est gérée par les callbacks ou signaux normalement

            # Remontrer la fenêtre principale (doit être fait sur le thread principal)
            # Utiliser un signal pour demander à la fenêtre principale de se remontrer
            # Exemple: self.show_main_window_signal.emit()
            # Pour l'instant, on suppose que la logique de rappel gère cela.

            self.running = False
            self.thread = None # Nettoyer la référence au thread

    # Ajouter un slot pour remonter la fenêtre principale si nécessaire
    # @pyqtSlot()
    # def handle_show_main_window(self):
    #     if self.progress_window:
    #         self.progress_window.close()
    #         self.progress_window = None
    #     self.show()

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
    """
    Capture toutes les exceptions non interceptées et affiche une QMessageBox.
    """
    # Format du traceback
    error_message = ''.join(traceback.format_exception(exctype, value, tb))

    # Écrire dans stderr pour le tracer même sans UI
    sys.stderr.write(f"Uncaught exception:\n{error_message}")

    # Afficher l'erreur dans une boîte de dialogue
    message_box = QMessageBox()
    message_box.setIcon(QMessageBox.Icon.Critical)
    message_box.setWindowTitle("Unhandled Exception")
    message_box.setText("An unexpected error occurred!")
    message_box.setDetailedText(error_message)
    message_box.exec()

    # Quitter l'application après l'erreur
    sys.exit(1)


# --- Fonction run_GUI (MODIFIÉE pour utiliser le threading) ---
def run_GUI():
    app_args = sys.argv if hasattr(sys, 'argv') else ['']
    app = QApplication(app_args)

    # Configuration de l'icône de l'application (inchangé)
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

    if splash_pixmap.isNull():
        print(f"Warning: Could not load splash image at {splash_pix_path}. Splash will lack main icon.")
        # Créer un pixmap vide pour éviter une erreur si LoadingSplashScreen en a besoin
        splash_pixmap = QPixmap(1, 1) # Petit pixmap valide mais invisible
        splash_pixmap.fill(Qt.GlobalColor.transparent)

    splash = LoadingSplashScreen(splash_pixmap)
    splash.showMessage("Loading FragHub, please wait...")
    splash.show()
    # Pas besoin de app.processEvents() ici, le threading gère la non-blocage

    # --- Threading Setup for Startup Tasks ---
    startup_thread = QThread()
    startup_worker = StartupWorker(BASE_DIR)
    startup_worker.moveToThread(startup_thread)

    # --- Variables pour stocker les résultats et références ---
    # Utilisation d'une liste mutable pour contourner les limitations de portée des nested functions
    shared_state = {'main_window': None, 'main_function': None, 'splash_screen': splash}

    # --- Slots (fonctions internes à run_GUI) ---
    def on_startup_complete(imported_main_function):
        """Appelé quand le worker a fini avec succès."""
        shared_state['main_function'] = imported_main_function

        # Créer la fenêtre principale MAINTENANT sur le thread principal
        try:
            shared_state['main_window'] = MainWindow(main_function_ref=shared_state['main_function'])
            shared_state['main_window'].show()
        except Exception as e:
             print(f"FATAL ERROR during MainWindow creation: {e}")
             traceback.print_exc()
             on_startup_error(f"Failed to create the main application window.\nError: {e}")
             return # Ne pas fermer le splash si l'erreur vient de MainWindow

        # Fermer le splash screen
        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None # Libérer la référence

    def on_startup_error(error_message):
        """Appelé si le worker rencontre une erreur."""
        print(f"Startup error: {error_message}")
        # Fermer le splash screen avant d'afficher l'erreur
        if shared_state['splash_screen']:
            shared_state['splash_screen'].close()
            shared_state['splash_screen'] = None

        # Afficher une boîte de dialogue d'erreur critique
        QMessageBox.critical(None, "Fatal Startup Error", error_message)
        app.quit() # Quitter l'application en cas d'erreur de démarrage

    def update_splash(message, font_size=28):
        """Met à jour le message sur le splash screen."""
        if shared_state['splash_screen']:
            shared_state['splash_screen'].showMessage(message, font_size=font_size)

    # --- Connect Signals and Slots ---
    startup_thread.started.connect(startup_worker.run_startup_tasks)
    startup_worker.finished.connect(on_startup_complete)
    startup_worker.error.connect(on_startup_error)
    startup_worker.update_splash_message.connect(update_splash)

    # Nettoyage du thread et du worker
    startup_worker.finished.connect(startup_thread.quit)
    startup_worker.finished.connect(startup_worker.deleteLater) # Planifie la suppression
    startup_thread.finished.connect(startup_thread.deleteLater) # Planifie la suppression

    # --- Démarrer le thread de démarrage ---
    startup_thread.start()

    # --- Démarrer la boucle d'événements principale ---
    # L'application attendra ici, traitant les événements (y compris l'animation du spinner)
    # jusqu'à ce que app.quit() soit appelé (par ex. en cas d'erreur ou fermeture normale)
    exit_code = app.exec()
    print(f"Application exiting with code {exit_code}")
    sys.exit(exit_code)

# --- Bloc d'exécution principal (Inchangé) ---
if __name__ == "__main__":
    run_GUI()