import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QLabel, QProgressBar, QTabWidget
)
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt, pyqtSignal
from threading import Thread
import time
from GUI.tabs.tab_input import InputTab
from GUI.tabs.tab_output import OutputTab
from GUI.tabs.tab_filters import FiltersTab
from GUI.tabs.tab_output_settings import OutputSettingTab
from GUI.tabs.tab_projects import ProjectsTab
from progress_window import ProgressWindow

from MAIN import MAIN


class MainWindow(QMainWindow):
    # Signal pour mettre à jour la barre de progression
    update_progress_signal = pyqtSignal(int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub 1.2.0")
        self.setWindowIcon(QIcon("./GUI/assets/FragHub_icon.png"))
        self.setGeometry(100, 100, 1280, 720)

        # Création du layout principal
        main_layout = QVBoxLayout()

        # Ajouter le logo FragHub en haut
        banner = QLabel()
        pixmap = QPixmap("./GUI/assets/FragHub_icon.png").scaled(
            130, 130, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
        )
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(banner)

        # Créez les onglets principaux
        self.tabs = QTabWidget()
        self.tabs.addTab(InputTab(), "INPUT")
        self.tabs.addTab(OutputTab(), "OUTPUT")
        self.tabs.addTab(FiltersTab(), "Filters settings")
        self.tabs.addTab(OutputSettingTab(), "Output settings")
        self.tabs.addTab(ProjectsTab(), "Projects settings")
        main_layout.addWidget(self.tabs)

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
        Affiche une nouvelle fenêtre pour la progression.
        """
        if not self.running:
            self.progress_window = ProgressWindow(self)
            self.progress_window.show()
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
        Exécute la fonction MAIN dans un thread.
        """
        try:
            MAIN(
                progress_callback=self.progress_window.update_progress_signal.emit,
                total_items_callback=self.progress_window.update_total_signal.emit,
                prefix_callback=self.progress_window.update_prefix_signal.emit,
                item_type_callback=self.progress_window.update_item_type_signal.emit,  # Nouveau callback
                step_callback=self.progress_window.update_step_signal.emit  # Ajout du step_callback
            )
        except Exception as e:
            print(f"Erreur : {e}")
        finally:
            self.running = False
            self.progress_window.close()
            self.progress_window = None


def run_GUI():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()


run_GUI()
