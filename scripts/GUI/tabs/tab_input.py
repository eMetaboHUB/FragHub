from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSpacerItem, QSizePolicy, QLabel, QHBoxLayout, QComboBox
from PyQt6.QtGui import QFont, QIcon, QPixmap
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtWidgets import QFileDialog
from ..utils.global_vars import parameters_dict  # Importer le dictionnaire global
import os

import sys

# Si le fichier est exécuté comme un exécutable PyInstaller
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # Si le fichier est exécuté comme un script Python
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class InputTab(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()

        # Ajouter des espaces pour centrer verticalement
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Créer le bouton de sélection de fichiers
        button = QPushButton()
        button.setIcon(QIcon(os.path.join(BASE_DIR,'../assets/files_icon.png')))
        button.setIconSize(QSize(128, 128))
        button.setFixedSize(140, 140)
        button.clicked.connect(self.browse_files)  # Connecter la fonction

        # Centrer le bouton
        button_layout = QHBoxLayout()
        button_layout.addWidget(button, alignment=Qt.AlignmentFlag.AlignCenter)
        self.layout.addLayout(button_layout)

        # Ajouter étiquette sous le bouton
        label = QLabel("Select input files")
        label.setFont(QFont("Arial", 14, QFont.Weight.Bold))
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(label)

        # Ajouter un menu déroulant sous le label
        self.file_menu = QComboBox()
        self.file_menu.setFixedWidth(200)  # Réduire la largeur à 200px
        self.file_menu.setPlaceholderText("No files selected")  # Texte initial
        self.layout.addWidget(self.file_menu, alignment=Qt.AlignmentFlag.AlignCenter)

        # Espacement et bouton d'informations
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        info_button_layout = QHBoxLayout()
        info_button_layout.addSpacerItem(QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum))

        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)
        info_button.setToolTip("Select single or multiple .json, .csv, .msp, or .mgf files")
        info_button_layout.addWidget(info_button, alignment=Qt.AlignmentFlag.AlignRight)
        self.layout.addLayout(info_button_layout)

        self.setLayout(self.layout)

    def browse_files(self):
        """Gestionnaire de sélection de fichiers pour l'onglet INPUT."""
        files, _ = QFileDialog.getOpenFileNames(self, "Choisir des fichiers")
        if files:
            # Mise à jour du dictionnaire global
            parameters_dict["input_directory"] = files

            # Ajouter les noms de fichiers au menu déroulant
            self.file_menu.clear()  # Nettoyer le menu existant
            for file_path in files:
                basename = os.path.basename(file_path)  # Récupérer juste le nom du fichier
                self.file_menu.addItem(basename)  # Ajouter au menu déroulant

            # Afficher directement le premier fichier si disponible
            if files:
                self.file_menu.setCurrentText(os.path.basename(files[0]))
