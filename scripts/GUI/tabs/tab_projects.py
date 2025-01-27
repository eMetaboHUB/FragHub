import os
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QSpacerItem, QSizePolicy, QLabel
from PyQt6.QtCore import Qt, pyqtSignal, QSize, QTimer
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen, QFont
from ..utils.global_vars import parameters_dict  # Importer le dictionnaire global


class QToggleSwitch(QWidget):
    stateChanged = pyqtSignal(bool)  # Signal émis quand l'état change

    def __init__(self, parent=None, initial_state=False):
        super().__init__(parent)
        self.setFixedSize(100, 50)  # Dimensions du toggle
        self._state = initial_state
        self._enabled = False  # La capacité de changer l'état est désactivée par défaut
        self.setCursor(Qt.CursorShape.ForbiddenCursor)  # Curseur changé pour interdire l'interaction au départ

    def mousePressEvent(self, event):
        if self._enabled:  # Permettre de basculer uniquement si le bouton est activé
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
        font = QFont("Arial", 12, QFont.Weight.Bold)
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
        return QSize(100, 50)

    def isChecked(self):
        return self._state

    def setChecked(self, state):
        self._state = state
        self.update()

    def setEnabled(self, enabled):
        self._enabled = enabled  # Active ou désactive la capacité d'interaction
        self.setCursor(Qt.CursorShape.PointingHandCursor if enabled else Qt.CursorShape.ForbiddenCursor)
        self.update()


class ProjectsTab(QWidget):
    output_directory_changed_signal = pyqtSignal(str)

    def __init__(self):
        super().__init__()

        if "reset_updates" not in parameters_dict:
            parameters_dict["reset_updates"] = 0.0  # Initialisation par défaut

        # Layout principal
        main_layout = QVBoxLayout()

        # Ajouter un espace extensible
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Ajouter une étiquette pour l'option "RESET PROJECT"
        reset_label = QLabel("RESET PROJECT ?")
        reset_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        reset_label.setStyleSheet("font-size: 32px; font-weight: bold; color: red;")
        main_layout.addWidget(reset_label)

        # Clignotement du label
        self.color_state = True
        self.timer = QTimer(self)
        self.timer.timeout.connect(lambda: self.toggle_label_color(reset_label))
        self.timer.start(500)

        # Ajouter le toggle switch
        self.toggle_switch = QToggleSwitch(initial_state=False)  # Toujours "NO" au départ
        main_layout.addWidget(self.toggle_switch, alignment=Qt.AlignmentFlag.AlignCenter)

        # Connecter le signal du toggle switch au gestionnaire
        self.toggle_switch.stateChanged.connect(self.on_toggle_state_changed)

        # Connecter le signal pour surveiller les changements de répertoire
        self.output_directory_changed_signal.connect(self.check_output_directory)

        # Ajouter un espace extensible après
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Ajouter un bouton d'information en bas à droite
        info_button_layout = QHBoxLayout()
        info_button_layout.setAlignment(Qt.AlignmentFlag.AlignRight)
        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)
        info_button.setToolTip(
            "Reseting the project deletes splash keys and output files from the selected project, to start fresh."
        )
        info_button_layout.addWidget(info_button)
        main_layout.addLayout(info_button_layout)

        # Appliquer le layout final
        self.setLayout(main_layout)

    def toggle_label_color(self, label):
        """Clignotement de l'étiquette RESET PROJECT."""
        if self.color_state:
            label.setStyleSheet("font-size: 32px; font-weight: bold; color: #FF4040;")
        else:
            label.setStyleSheet("font-size: 32px; font-weight: bold; color: #8B0000;")
        self.color_state = not self.color_state

    def check_output_directory(self, directory):
        """Vérifier la présence du fichier `.fraghub` dans le répertoire sélectionné."""
        if directory:
            fraghub_file = os.path.join(directory, ".fraghub")
            if os.path.isfile(fraghub_file):
                self.toggle_switch.setEnabled(True)  # Débloquer le toggle switch
            else:
                self.toggle_switch.setEnabled(False)  # Bloquer le toggle switch
                self.toggle_switch.setChecked(False)  # Rester sur "NO" par défaut

    def on_toggle_state_changed(self, state):
        """
        Gère les changements d'état du toggle switch.
        Met à jour le dictionnaire global pour refléter l'état actuel.
        """
        parameters_dict["reset_updates"] = 1.0 if state else 0.0  # 1.0 pour YES, 0.0 pour NO
