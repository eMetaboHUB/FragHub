from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSpacerItem, QSizePolicy
from PyQt6.QtCore import Qt, QSize, pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen, QFont
from utils.global_vars import parameters_dict  # Importer le dictionnaire global


# Classe pour le Toggle Switch personnalisé avec style moderne
class QToggleSwitch(QWidget):
    # Signal émis lorsque l'état change
    stateChanged = pyqtSignal(bool)

    def __init__(self, parent=None, initial_state=True):
        super().__init__(parent)
        self.setFixedSize(60, 30)  # Taille du switch
        self._state = initial_state  # État initial
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event):
        # Inverser l'état au clic
        self._state = not self._state
        self.update()  # Redessiner le widget
        self.stateChanged.emit(self._state)  # Émettre le signal de changement d'état

    def paintEvent(self, event):
        # Peindre l'interrupteur
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        # Dimensions du rectangle
        rect = self.rect()

        # Couleur d'arrière-plan basée sur l'état
        background_color = QColor("#00C853") if self._state else QColor("#f44336")
        painter.setBrush(QBrush(background_color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRoundedRect(rect, rect.height() // 2, rect.height() // 2)

        # Dessiner le texte ('YES' ou 'NO') légèrement décalé
        font = QFont("Arial", 8, QFont.Weight.Bold)  # Réduction de la taille de la police
        painter.setFont(font)
        text_color = QColor("#FFFFFF")  # Blanc pour le texte

        # Texte et positions
        if self._state:
            text = "YES"
            text_pos = rect.left() + 5  # "YES" encore plus à gauche
        else:
            text = "NO"
            text_pos = rect.right() - 25  # "NO" encore plus à droite

        painter.setPen(text_color)
        painter.drawText(
            text_pos, rect.height() // 2 + 5, text
        )  # Placer le texte à la bonne position

        # Dessiner le bouton rond glissant
        button_color = QColor("#FFFFFF")  # Bouton blanc
        button_x = rect.width() - rect.height() + 2 if self._state else 2
        button_y = 2
        painter.setBrush(QBrush(button_color))
        painter.setPen(QPen(QColor("#E0E0E0"), 1))  # Bordure subtile
        painter.drawEllipse(button_x, button_y, rect.height() - 4, rect.height() - 4)

    def sizeHint(self):
        return QSize(60, 30)

    def isChecked(self):
        return self._state

    def setChecked(self, state):
        self._state = state
        self.update()


# Classe pour l'onglet "Output settings"
class OutputSettingTab(QWidget):
    def __init__(self):
        super().__init__()

        # Initialisation des paramètres par défaut dans le dictionnaire global avec des valeurs numériques
        parameters_dict['csv'] = 1.0
        parameters_dict['msp'] = 1.0
        parameters_dict['json'] = 1.0

        # Définir le layout principal (vertical global)
        main_layout = QVBoxLayout()

        # Ajouter un espace extensible avant le contenu principal
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Layout centré pour les éléments centraux (switches et labels)
        center_layout = QVBoxLayout()
        center_layout.setAlignment(
            Qt.AlignmentFlag.AlignCenter)  # Centrer le contenu principal horizontalement et verticalement

        # Ajouter les lignes pour chaque toggle switch avec leurs labels
        center_layout.addLayout(self._create_row("CSV", "csv"))
        center_layout.addLayout(self._create_row("MSP", "msp"))
        center_layout.addLayout(self._create_row("JSON", "json"))

        # Ajouter le layout centré au layout principal
        main_layout.addLayout(center_layout)

        # Ajouter un espace extensible après le contenu principal
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Ajouter le bouton d'information en bas à droite
        info_button_layout = QHBoxLayout()
        info_button_layout.setAlignment(Qt.AlignmentFlag.AlignRight)  # Aligner le bouton à droite
        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)

        # Définir une info-bulle pour le bouton
        info_button.setToolTip("This tab lets you chose the output formats to be written by FragHub at the end of processing.")

        # Ajouter le bouton au layout en bas à droite
        info_button_layout.addWidget(info_button)

        # Ajouter le layout du bouton au layout principal
        main_layout.addLayout(info_button_layout)

        # Appliquer le layout principal à l'onglet
        self.setLayout(main_layout)

    def _create_row(self, label_text, parameter_key):
        # Crée une ligne avec un toggle switch et un label
        row_layout = QHBoxLayout()  # Disposition horizontale
        row_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Centrer chaque ligne

        # Créer le toggle switch et connecter son signal
        toggle = QToggleSwitch(initial_state=parameters_dict[parameter_key] == 1.0)
        toggle.stateChanged.connect(lambda state: self.update_parameter(parameter_key, state))

        # Créer le label pour le texte
        label = QLabel(label_text)
        label.setStyleSheet("font-size: 14px;")  # Augmenter légèrement la taille du texte

        # Ajouter le switch et le texte dans la ligne
        row_layout.addWidget(toggle)
        row_layout.addWidget(label)

        return row_layout

    def update_parameter(self, key, state):
        # Mettre à jour le dictionnaire global avec 1.0 ou 0.0 selon l'état
        parameters_dict[key] = 1.0 if state else 0.0