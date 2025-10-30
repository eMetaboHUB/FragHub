from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSpacerItem, QSizePolicy
from PyQt6.QtCore import Qt, QSize, pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen, QFont
from scripts.GUI.utils.global_vars import parameters_dict  # Importer le dictionnaire global


# Class for the custom Toggle Switch with modern styling
class QToggleSwitch(QWidget):
    # Signal emitted when the state changes
    stateChanged = pyqtSignal(bool)

    def __init__(self, parent=None, initial_state=True):
        super().__init__(parent)
        self.setFixedSize(60, 30)  # Switch size
        self._state = initial_state  # Initial state
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event):
        # Toggle state on click
        self._state = not self._state
        self.update()  # Redraw the widget
        self.stateChanged.emit(self._state)  # Emit the state change signal

    def paintEvent(self, event):
        # Paint the switch
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        # Rectangle dimensions
        rect = self.rect()

        # Background color based on the state
        background_color = QColor("#00C853") if self._state else QColor("#f44336")
        painter.setBrush(QBrush(background_color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRoundedRect(rect, rect.height() // 2, rect.height() // 2)

        # Draw the text ('YES' or 'NO') slightly adjusted
        font = QFont("Arial", 8, QFont.Weight.Bold)  # Smaller font size
        painter.setFont(font)
        text_color = QColor("#FFFFFF")  # White for the text

        # Text and positions
        if self._state:
            text = "YES"
            text_pos = rect.left() + 5  # "YES" slightly more to the left
        else:
            text = "NO"
            text_pos = rect.right() - 25  # "NO" slightly more to the right

        painter.setPen(text_color)
        painter.drawText(
            text_pos, rect.height() // 2 + 5, text
        )  # Place the text in the correct position

        # Draw the sliding round button
        button_color = QColor("#FFFFFF")  # White button
        button_x = rect.width() - rect.height() + 2 if self._state else 2
        button_y = 2
        painter.setBrush(QBrush(button_color))
        painter.setPen(QPen(QColor("#E0E0E0"), 1))  # Subtle border
        painter.drawEllipse(button_x, button_y, rect.height() - 4, rect.height() - 4)

    def sizeHint(self):
        return QSize(60, 30)

    def isChecked(self):
        return self._state

    def setChecked(self, state):
        self._state = state
        self.update()


# Class for the "Output settings" tab
class OutputSettingTab(QWidget):
    def __init__(self):
        super().__init__()

        # Initialize default parameters in the global dictionary with numeric values
        parameters_dict['csv'] = 1.0
        parameters_dict['msp'] = 1.0
        parameters_dict['json'] = 1.0

        # Define the main layout (global vertical layout)
        main_layout = QVBoxLayout()

        # Add expandable space before the main content
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Centered layout for central elements (switches and labels)
        center_layout = QVBoxLayout()
        center_layout.setAlignment(
            Qt.AlignmentFlag.AlignCenter)  # Center the main content horizontally and vertically

        # Add rows for each toggle switch with their labels
        center_layout.addLayout(self._create_row("CSV", "csv"))
        center_layout.addLayout(self._create_row("MSP", "msp"))
        center_layout.addLayout(self._create_row("JSON", "json"))

        # Add the centered layout to the main layout
        main_layout.addLayout(center_layout)

        # Add expandable space after the main content
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Add an information button at the bottom-right corner
        info_button_layout = QHBoxLayout()
        info_button_layout.setAlignment(Qt.AlignmentFlag.AlignRight)  # Align the button to the right
        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)

        # Set a tooltip for the button
        info_button.setToolTip("This tab lets you choose the output formats to be written by FragHub at the end of processing.")

        # Add the button to the bottom-right layout
        info_button_layout.addWidget(info_button)

        # Add the button layout to the main layout
        main_layout.addLayout(info_button_layout)

        # Apply the main layout to the tab
        self.setLayout(main_layout)

    def _create_row(self, label_text, parameter_key):
        # Create a row with a toggle switch and a label
        row_layout = QHBoxLayout()  # Horizontal layout
        row_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center each row

        # Create the toggle switch and connect its signal
        toggle = QToggleSwitch(initial_state=parameters_dict[parameter_key] == 1.0)
        toggle.stateChanged.connect(lambda state: self.update_parameter(parameter_key, state))

        # Create the label for the text
        label = QLabel(label_text)
        label.setStyleSheet("font-size: 14px;")  # Slightly increase the text size

        # Add the switch and the text to the row
        row_layout.addWidget(toggle)
        row_layout.addWidget(label)

        return row_layout

    def update_parameter(self, key, state):
        # Update the global dictionary with 1.0 or 0.0 depending on the state
        parameters_dict[key] = 1.0 if state else 0.0