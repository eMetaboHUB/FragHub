from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QSpacerItem,
    QSizePolicy, QPushButton
)
from PyQt6.QtGui import QFont, QPainter, QColor, QBrush, QPen
from PyQt6.QtCore import pyqtSignal, Qt

# Import the global parameters dictionary from your project
from scripts.GUI.utils.global_vars import parameters_dict


# --- QToggleSwitch Class (identical to the one in tab_filters.py) ---
# A custom ON/OFF switch for a clear user interface.
class QToggleSwitch(QWidget):
    """
    A custom ON/OFF switch widget.
    Emits a stateChanged(bool) signal when its state changes.
    """
    stateChanged = pyqtSignal(bool)

    def __init__(self, parent=None, initial_state=False):
        super().__init__(parent)
        self.setFixedSize(60, 30)  # Compact size
        self._state = initial_state
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event):
        self._state = not self._state
        self.update()
        self.stateChanged.emit(self._state)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        rect = self.rect()
        background_color = QColor("#00C853") if self._state else QColor("#f44336")
        painter.setBrush(QBrush(background_color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRoundedRect(rect, rect.height() // 2, rect.height() // 2)

        font = QFont("Arial", 8, QFont.Weight.Bold)
        painter.setFont(font)
        text_color = QColor("#FFFFFF")
        painter.setPen(text_color)

        if self._state:
            text = "ON"
            text_pos = rect.left() + 12
        else:
            text = "OFF"
            text_pos = rect.right() - 32

        painter.drawText(text_pos, rect.height() // 2 + 5, text)

        button_color = QColor("#FFFFFF")
        button_x = rect.width() - rect.height() + 2 if self._state else 2
        button_y = 2
        painter.setBrush(QBrush(button_color))
        painter.setPen(QPen(QColor("#E0E0E0"), 1))
        painter.drawEllipse(button_x, button_y, rect.height() - 4, rect.height() - 4)

    def isChecked(self):
        return self._state

    def setChecked(self, state):
        if self._state != state:
            self._state = state
            self.update()


# --- Main Class for the De Novo Tab ---
class DeNovoTab(QWidget):
    """
    GUI tab to configure de novo calculation options.
    Contains controls to enable the feature and adjust its parameters.
    """

    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout(self)

        self.initialize_parameters()
        self.create_de_novo_options()
        self.setLayout(self.layout)

    def initialize_parameters(self):
        """Initializes the default values in parameters_dict."""
        global parameters_dict
        parameters_dict.setdefault("calculate_de_novo", 0.0)
        parameters_dict.setdefault("de_novo_ppm_tolerance", 10.0)

    def create_de_novo_options(self):
        """Creates and arranges the tab's widgets."""

        # Spacer to push content to the vertical center
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # --- Row 1: De Novo Calculation Activation (horizontally centered) ---
        de_novo_layout = QHBoxLayout()
        de_novo_layout.addStretch()

        initial_state = parameters_dict.get("calculate_de_novo", 0.0) == 1.0
        self.toggle_de_novo = QToggleSwitch(initial_state=initial_state)
        self.toggle_de_novo.stateChanged.connect(self.on_toggle_changed)
        de_novo_layout.addWidget(self.toggle_de_novo)

        label = QLabel("Calculate de novo formulas")
        label.setFont(QFont("Arial", 12))
        de_novo_layout.addWidget(label)

        warning_label = QLabel("(warning not compatible with mzMine)")
        warning_font = QFont("Arial", 9)
        warning_font.setItalic(True)
        warning_label.setFont(warning_font)
        warning_label.setStyleSheet("color: #888;")
        de_novo_layout.addWidget(warning_label)

        de_novo_layout.addStretch()
        self.layout.addLayout(de_novo_layout)

        # --- Row 2: PPM Tolerance Parameter (horizontally centered) ---
        ppm_layout = QHBoxLayout()
        ppm_layout.addStretch()

        ppm_label = QLabel("ppm tolerance:")
        ppm_label.setFont(QFont("Arial", 10))
        ppm_layout.addWidget(ppm_label)

        initial_ppm = str(parameters_dict.get("de_novo_ppm_tolerance", 10.0))
        self.ppm_field = QLineEdit(initial_ppm)
        self.ppm_field.setFixedWidth(60)
        self.ppm_field.textChanged.connect(
            lambda value: self.handle_text_change("de_novo_ppm_tolerance", value)
        )
        ppm_layout.addWidget(self.ppm_field)

        ppm_layout.addStretch()
        self.layout.addLayout(ppm_layout)

        # Spacer that pushes the info button down
        self.layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # --- Info button in the bottom right ---
        info_button_layout = QHBoxLayout()
        info_button_layout.addStretch()  # Pushes the button to the right

        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)
        info_button.setToolTip(
            "Enables the de novo chemical formula calculation for each spectrum.\n"
            "The PPM tolerance is used for the precision of the formula matching.\n"
            "Warning, this option makes the output databases incompatible with most reprocessing software."
        )
        info_button_layout.addWidget(info_button)
        self.layout.addLayout(info_button_layout)

    def on_toggle_changed(self, state):
        """Updates the dictionary when the switch is toggled."""
        global parameters_dict
        parameters_dict["calculate_de_novo"] = 1.0 if state else 0.0

    def handle_text_change(self, key, value):
        """Updates the dictionary when a text field is changed."""
        global parameters_dict
        try:
            parameters_dict[key] = float(value)
        except ValueError:
            pass
