import os
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QComboBox, QHBoxLayout, QStyledItemDelegate, QPushButton, QSpacerItem, \
    QSizePolicy, QLabel
from PyQt6.QtCore import Qt, pyqtSignal, QSize, QTimer
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen, QFont
from utils.global_vars import parameters_dict  # Importer le dictionnaire global


class QToggleSwitch(QWidget):
    stateChanged = pyqtSignal(bool)  # Signal émis quand l'état change

    def __init__(self, parent=None, initial_state=False):
        super().__init__(parent)
        self.setFixedSize(100, 50)  # Dimensions augmentées du toggle switch
        self._state = initial_state
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event):
        self._state = not self._state  # Inverser l'état
        self.update()  # Redessiner le switch
        self.stateChanged.emit(self._state)  # Émettre le signal

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        rect = self.rect()
        background_color = QColor("#00C853") if self._state else QColor("#f44336")
        painter.setBrush(QBrush(background_color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRoundedRect(rect, rect.height() // 2, rect.height() // 2)
        font = QFont("Arial", 12, QFont.Weight.Bold)  # Police plus grande pour le texte "YES/NO"
        painter.setFont(font)
        text_color = QColor("#FFFFFF")
        painter.setPen(text_color)
        text = "YES" if self._state else "NO"
        text_pos = rect.left() + 10 if self._state else rect.right() - 35
        painter.drawText(text_pos, int(rect.height() * 0.65), text)
        button_color = QColor("#FFFFFF")
        button_x = rect.width() - rect.height() + 5 if self._state else 5
        painter.setBrush(QBrush(button_color))
        painter.setPen(QPen(QColor("#E0E0E0"), 1))
        painter.drawEllipse(button_x, 5, rect.height() - 10, rect.height() - 10)

    def sizeHint(self):
        return QSize(100, 50)  # Taille adaptée

    def isChecked(self):
        return self._state

    def setChecked(self, state):
        self._state = state
        self.update()


class CenteredItemDelegate(QStyledItemDelegate):
    """Délégué personnalisé pour centrer les éléments dans le menu déroulant."""

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        option.displayAlignment = Qt.AlignmentFlag.AlignCenter


class ProjectsTab(QWidget):
    def __init__(self):
        super().__init__()

        # Ajouter un layout principal
        main_layout = QVBoxLayout()

        # Ajouter un espace extensible avant le contenu principal
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Ajouter une étiquette pour la section "Projects"
        project_label = QLabel("SELECT PROJECT")
        project_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        project_label.setStyleSheet("font-size: 24px; font-weight: bold; color: black;")
        main_layout.addWidget(project_label)

        # Section Menu déroulant
        center_layout = QHBoxLayout()
        self.dropdown = QComboBox()
        self.dropdown.addItem("Basic")  # Projet par défaut
        self.dropdown.setItemDelegate(CenteredItemDelegate(self.dropdown))
        self.dropdown.setEditable(True)
        self.dropdown.lineEdit().setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.dropdown.setFixedWidth(200)  # Largeur de 200px pour le menu déroulant
        self.populate_dropdown_with_json_files()
        self.dropdown.editTextChanged.connect(self.handle_new_project)
        self.dropdown.currentTextChanged.connect(self.save_selected_project)
        center_layout.addWidget(self.dropdown)
        main_layout.addLayout(center_layout)

        # Ajouter un séparateur
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Ajouter une étiquette pour la section "Reset Project"
        reset_label = QLabel("RESET PROJECT ?")
        reset_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        reset_label.setStyleSheet("font-size: 32px; font-weight: bold; color: red;")
        main_layout.addWidget(reset_label)

        # Démarrage du clignotement du texte
        self.color_state = True
        self.timer = QTimer(self)
        self.timer.timeout.connect(lambda: self.toggle_label_color(reset_label))
        self.timer.start(500)

        # Ajouter un toggle switch en dessous de l'étiquette
        initial_state = parameters_dict.get("reset_updates", 0.0) == 1.0
        self.toggle_switch = QToggleSwitch(initial_state=initial_state)
        self.toggle_switch.stateChanged.connect(self.handle_toggle)
        main_layout.addWidget(self.toggle_switch, alignment=Qt.AlignmentFlag.AlignCenter)

        # Ajouter un espace extensible après le contenu principal
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Ajouter un bouton d'information général en bas à droite
        info_button_layout = QHBoxLayout()
        info_button_layout.setAlignment(Qt.AlignmentFlag.AlignRight)
        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)
        info_button.setToolTip(
            "FragHub projects allow you to save the splash keys of previous spectra processed.\nSo that at the next update, only new spectra from the database are processed, and added to previous FragHub processes.\nEnter the name of your project to create a new one, or select an existing one.\n\nReseting project allows you to delete splash keys and output files from the selected project, in order to start a new project from scratch."
        )
        info_button_layout.addWidget(info_button)
        main_layout.addLayout(info_button_layout)

        # Appliquer le layout final
        self.setLayout(main_layout)

    def populate_dropdown_with_json_files(self):
        """Recherche et ajoute les fichiers .json depuis ../../datas/updates dans le menu déroulant."""
        json_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../datas/updates"))
        if os.path.exists(json_dir) and os.path.isdir(json_dir):
            json_files = [file for file in os.listdir(json_dir) if file.endswith(".json")]
            if json_files:
                self.dropdown.addItems([os.path.splitext(file)[0] for file in json_files])
            else:
                self.dropdown.addItem("Aucun fichier JSON trouvé")
        else:
            self.dropdown.addItem("Répertoire introuvable")

    def handle_new_project(self, text):
        """Ajoute un nouveau projet si non existant et le sélectionne."""
        if text and text not in [self.dropdown.itemText(i) for i in range(self.dropdown.count())]:
            self.dropdown.addItem(text)
            self.dropdown.setCurrentText(text)
        self.save_selected_project(text)

    def save_selected_project(self, text):
        """Enregistre le projet sélectionné."""
        parameters_dict["selected_profile"] = text

    def toggle_label_color(self, label):
        """Alterner la couleur de l'étiquette "RESET PROJECT ?"."""
        if self.color_state:
            label.setStyleSheet("font-size: 32px; font-weight: bold; color: #FF4040;")
        else:
            label.setStyleSheet("font-size: 32px; font-weight: bold; color: #8B0000;")
        self.color_state = not self.color_state

    def handle_toggle(self, state):
        """Traiter le changement d'état du toggle switch."""
        parameters_dict["reset_updates"] = 1.0 if state else 0.0
