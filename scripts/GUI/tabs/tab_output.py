from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSpacerItem, QSizePolicy, QLabel, QHBoxLayout
from PyQt6.QtGui import QFont, QIcon
from PyQt6.QtCore import QSize, Qt, pyqtSignal
from PyQt6.QtWidgets import QFileDialog
from ..utils.global_vars import parameters_dict  # Importer le dictionnaire global
import sys
import os

# Si le fichier est exécuté comme un exécutable PyInstaller
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # Si le fichier est exécuté comme un script Python
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

class OutputTab(QWidget):
    output_directory_changed = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()

        # Ajouter des espaces pour centrer verticalement
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Créer le bouton de sélection de répertoires
        button = QPushButton()
        button.setIcon(QIcon(os.path.join(BASE_DIR,'../assets/directory.png')))
        button.setIconSize(QSize(128, 128))
        button.setFixedSize(140, 140)
        button.clicked.connect(self.browse_output_files)  # Connecter le bouton à la fonction

        # Centrer le bouton
        button_layout = QHBoxLayout()
        button_layout.addWidget(button, alignment=Qt.AlignmentFlag.AlignCenter)
        self.layout.addLayout(button_layout)

        # Ajouter l'étiquette sous le bouton
        label = QLabel("Select output directory")
        label.setFont(QFont("Arial", 14, QFont.Weight.Bold))
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(label)

        # Espacement et bouton d'infos
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        info_button_layout = QHBoxLayout()
        info_button_layout.addSpacerItem(QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum))

        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)
        info_button.setToolTip("Create a new empty directory or Select an existing directory where FragHub has already written files")
        info_button_layout.addWidget(info_button, alignment=Qt.AlignmentFlag.AlignRight)
        self.layout.addLayout(info_button_layout)

        self.setLayout(self.layout)

    def browse_output_files(self):
        """Gestionnaire de sélection de répertoire pour l'onglet OUTPUT."""
        directory = QFileDialog.getExistingDirectory(self, "Choisir un répertoire pour OUTPUT")
        if directory:
            parameters_dict["output_directory"] = directory
            self.output_directory_changed.emit(directory)  # Émettre le signal
