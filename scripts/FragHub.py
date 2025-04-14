import sys
import os

# Si le fichier est exécuté comme un exécutable PyInstaller
if getattr(sys, 'frozen', False):
    BASE_DIR = os.path.dirname(sys.executable)
else:
    # Si le fichier est exécuté comme un script Python
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QProgressBar, QTabWidget
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal
from threading import Thread
import time
from scripts.GUI.tabs.tab_input import InputTab
from scripts.GUI.tabs.tab_output import OutputTab
from scripts.GUI.tabs.tab_filters import FiltersTab
from scripts.GUI.tabs.tab_output_settings import OutputSettingTab
from scripts.GUI.tabs.tab_projects import ProjectsTab
from scripts.progress_window import ProgressWindow

from scripts.MAIN import MAIN
import traceback
import ctypes

ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FragHub")


class MainWindow(QMainWindow):
    # Signal pour mettre à jour la barre de progression
    update_progress_signal = pyqtSignal(int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub 1.3.0")
        self.setWindowIcon(QIcon("./GUI/assets/FragHub_icon.png"))
        self.setGeometry(100, 100, 1280, 720)

        # Création du layout principal
        main_layout = QVBoxLayout()

        # Ajouter le logo FragHub en haut
        banner = QLabel()
        pixmap = QPixmap("./GUI/assets/FragHub_icon.png").scaled(
            200, 200, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
        )
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)

        # Créez les onglets principaux
        self.tabs = QTabWidget()

        # Créer des instances explicites de OutputTab et ProjectsTab
        self.output_tab = OutputTab()
        self.projects_tab = ProjectsTab()

        # Ajouter les instances d'onglet à QTabWidget
        self.tabs.addTab(InputTab(), "INPUT")
        self.tabs.addTab(self.output_tab, "OUTPUT")
        self.tabs.addTab(FiltersTab(), "Filters settings")
        self.tabs.addTab(OutputSettingTab(), "Output settings")
        self.tabs.addTab(self.projects_tab, "Projects settings")
        main_layout.addWidget(self.tabs)

        # Connecter le signal (après avoir ajouté les onglets)
        self.output_tab.output_directory_changed.connect(self.projects_tab.output_directory_changed_signal)

        # Bouton START/STOP
        self.start_button = QPushButton("START")
        self.start_button.setFixedSize(140, 60)
        self.start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        self.start_button.clicked.connect(self.open_progress_window)
        main_layout.addWidget(self.start_button, alignment=Qt.AlignmentFlag.AlignCenter)

        # Ajouter le layout principal au widget central
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        # Variables de contrôle
        self.running = False
        self.progress_window = None
        self.thread = None

    def open_progress_window(self):
        """
        Affiche la fenêtre de progression et masque la fenêtre principale (indépendante).
        """
        if not self.running:
            self.hide()  # Cache simplement la fenêtre principale.

            # Crée une nouvelle fenêtre de progression sans parent
            self.progress_window = ProgressWindow(parent=None)
            self.progress_window.show()  # Affiche la fenêtre de progression

            # Commencer le processus d'exécution
            self.start_execution()

    def start_execution(self):
        """
        Lance l'exécution de MAIN() dans un thread séparé.
        """
        self.running = True
        self.thread = Thread(target=self.run_main_function, daemon=True)
        self.thread.start()

    def run_main_function(self):
        """
        Exécution de la tâche principale avec gestion de la progression.
        """
        try:
            # Lancement du processus principal
            MAIN(
                progress_callback=self.progress_window.update_progress_signal.emit,
                total_items_callback=self.progress_window.update_total_signal.emit,
                prefix_callback=self.progress_window.update_prefix_signal.emit,
                item_type_callback=self.progress_window.update_item_type_signal.emit,
                step_callback=self.progress_window.update_step_signal.emit,
                completion_callback=self.progress_window.completion_callback.emit,
                deletion_callback=self.progress_window.deletion_callback.emit,
            )
        except Exception as e:
            traceback.print_exc()
        finally:
            # Restaurer la fenêtre principale et terminer la progression
            self.show()  # Réaffiche la fenêtre principale
            if self.progress_window:
                self.progress_window.close()  # Ferme la fenêtre de progression
                self.progress_window = None
            self.running = False


def run_GUI():
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("./GUI/assets/FragHub_Python_icon.ico"))
    window = MainWindow()
    window.show()
    app.exec()


run_GUI()
