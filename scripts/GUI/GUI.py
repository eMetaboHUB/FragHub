import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget, QPushButton, QLabel
from PyQt6.QtGui import QFont, QPixmap, QIcon
from PyQt6.QtCore import Qt

from .tabs.tab_input import InputTab
from .tabs.tab_output import OutputTab
from .tabs.tab_filters import FiltersTab
from .tabs.tab_output_settings import OutputSettingTab
from .tabs.tab_projects import ProjectsTab
from .utils.global_vars import parameters_dict


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FragHub 1.2.0")
        self.setWindowIcon(QIcon('./GUI/assets/FragHub_icon.png'))
        self.setGeometry(100, 100, 900, 600)

        # Création du layout principal
        main_layout = QVBoxLayout()

        # Ajouter le logo FragHub en haut
        banner = QLabel()
        pixmap = QPixmap('./GUI/assets/FragHub_icon.png').scaled(130, 130, Qt.AspectRatioMode.KeepAspectRatio,
                                                           Qt.TransformationMode.SmoothTransformation)
        banner.setPixmap(pixmap)
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Centrer l'image
        main_layout.addWidget(banner)  # Ajouter au layout principal

        # Créez les onglets principaux
        self.tabs = QTabWidget()
        self.tabs.addTab(InputTab(), "INPUT")
        self.tabs.addTab(OutputTab(), "OUTPUT")
        self.tabs.addTab(FiltersTab(), "Filters settings")
        self.tabs.addTab(OutputSettingTab(), "Output settings")
        self.tabs.addTab(ProjectsTab(), "Projects settings")  # Onglet déjà existant
        main_layout.addWidget(self.tabs)

        # Ajouter le bouton START en bas
        start_button = QPushButton("START")
        start_button.setFixedSize(140, 60)
        start_button.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        start_button.clicked.connect(self.start_action)  # Connecter le clic à la méthode start_action
        main_layout.addWidget(start_button, alignment=Qt.AlignmentFlag.AlignCenter)

        # Ajouter ce layout principal au widget central
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

    def start_action(self):
        """
        Action exécutée lors du clic sur START.
        Ferme la fenêtre et continue le programme.
        """
        self.close()  # Fermer la fenêtre principale


def run_GUI():
    """Lancer l'application GUI."""
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()


# if __name__ == "__main__":
#     run_GUI()
#     print(parameters_dict)
