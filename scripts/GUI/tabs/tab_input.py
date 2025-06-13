from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QSpacerItem, QSizePolicy, QLabel, QHBoxLayout, QComboBox
from PyQt6.QtGui import QFont, QIcon, QPixmap
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtWidgets import QFileDialog
from scripts.GUI.utils.global_vars import parameters_dict  # Importer le dictionnaire global
import os

import sys

# If the file is executed as a PyInstaller executable
if getattr(sys, 'frozen', False):
    BASE_DIR = sys._MEIPASS
else:
    # If the file is executed as a Python script
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


class InputTab(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()

        # Add space to vertically center elements
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Create the file selection button
        button = QPushButton()
        button.setIcon(QIcon(os.path.join(BASE_DIR, './GUI/assets/files_icon.png')))
        button.setIconSize(QSize(128, 128))
        button.setFixedSize(140, 140)
        button.clicked.connect(self.browse_files)  # Connect the function

        # Center the button
        button_layout = QHBoxLayout()
        button_layout.addWidget(button, alignment=Qt.AlignmentFlag.AlignCenter)
        self.layout.addLayout(button_layout)

        # Add a label below the button
        label = QLabel("Select input files")
        label.setFont(QFont("Arial", 14, QFont.Weight.Bold))
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(label)

        # Add a dropdown menu below the label
        self.file_menu = QComboBox()
        self.file_menu.setFixedWidth(200)  # Set the width to 200px
        self.file_menu.setPlaceholderText("No files selected")  # Initial text
        self.layout.addWidget(self.file_menu, alignment=Qt.AlignmentFlag.AlignCenter)

        # Add spacing and an info button
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
        """File selection handler for the INPUT tab."""
        files, _ = QFileDialog.getOpenFileNames(self, "Choose files")
        if files:
            # Update the global dictionary
            parameters_dict["input_directory"] = files

            # Add file names to the dropdown menu
            self.file_menu.clear()  # Clear existing items
            for file_path in files:
                basename = os.path.basename(file_path)  # Get only the file name
                self.file_menu.addItem(basename)  # Add to the dropdown

            # Automatically show the first file if available
            if files:
                self.file_menu.setCurrentText(os.path.basename(files[0]))