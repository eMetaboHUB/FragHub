import os
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QSpacerItem, QSizePolicy, QLabel
from PyQt6.QtCore import Qt, pyqtSignal, QSize, QTimer
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen, QFont
from scripts.GUI.utils.global_vars import parameters_dict  # Importer le dictionnaire global


class QToggleSwitch(QWidget):
    stateChanged = pyqtSignal(bool)  # Signal emitted when the state changes

    def __init__(self, parent=None, initial_state=False):
        super().__init__(parent)
        self.setFixedSize(100, 50)  # Toggle dimensions
        self._state = initial_state
        self._enabled = False  # The ability to change the state is disabled by default
        self.setCursor(Qt.CursorShape.ForbiddenCursor)  # Cursor set to indicate interaction is forbidden by default

    def mousePressEvent(self, event):
        if self._enabled:  # Allow toggling only if the button is enabled
            self._state = not self._state  # Toggle the state
            self.update()  # Redraw the switch
            self.stateChanged.emit(self._state)  # Emit the signal

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
        self._enabled = enabled  # Enables or disables interaction capability
        self.setCursor(Qt.CursorShape.PointingHandCursor if enabled else Qt.CursorShape.ForbiddenCursor)
        self.update()


class ProjectsTab(QWidget):
    output_directory_changed_signal = pyqtSignal(str)

    def __init__(self):
        super().__init__()

        if "reset_updates" not in parameters_dict:
            parameters_dict["reset_updates"] = 0.0  # Default initialization

        # Main layout
        main_layout = QVBoxLayout()

        # Add an expandable space
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Add a label for the "RESET PROJECT" option
        reset_label = QLabel("RESET PROJECT ?")
        reset_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        reset_label.setStyleSheet("font-size: 32px; font-weight: bold; color: red;")
        main_layout.addWidget(reset_label)

        # Blinking effect for the label
        self.color_state = True
        self.timer = QTimer(self)
        self.timer.timeout.connect(lambda: self.toggle_label_color(reset_label))
        self.timer.start(500)

        # Add the toggle switch
        self.toggle_switch = QToggleSwitch(initial_state=False)  # Always starts with "NO"
        main_layout.addWidget(self.toggle_switch, alignment=Qt.AlignmentFlag.AlignCenter)

        # Connect the toggle switch's signal to the handler
        self.toggle_switch.stateChanged.connect(self.on_toggle_state_changed)

        # Connect the signal to monitor directory changes
        self.output_directory_changed_signal.connect(self.check_output_directory)

        # Add an expandable space after
        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        # Add an information button at the bottom-right
        info_button_layout = QHBoxLayout()
        info_button_layout.setAlignment(Qt.AlignmentFlag.AlignRight)
        info_button = QPushButton("🛈")
        info_button.setFixedSize(30, 30)
        info_button.setToolTip(
            "Resetting the project deletes splash keys and output files from the selected project, to start fresh."
        )
        info_button_layout.addWidget(info_button)
        main_layout.addLayout(info_button_layout)

        # Apply the final layout
        self.setLayout(main_layout)

    def toggle_label_color(self, label):
        """Blinking effect for the RESET PROJECT label."""
        if self.color_state:
            label.setStyleSheet("font-size: 32px; font-weight: bold; color: #FF4040;")
        else:
            label.setStyleSheet("font-size: 32px; font-weight: bold; color: #8B0000;")
        self.color_state = not self.color_state

    def check_output_directory(self, directory):
        """Check the presence of the `.fraghub` file in the selected directory."""
        if directory:
            fraghub_file = os.path.join(directory, ".fraghub")
            if os.path.isfile(fraghub_file):
                self.toggle_switch.setEnabled(True)  # Enable the toggle switch
            else:
                self.toggle_switch.setEnabled(False)  # Disable the toggle switch
                self.toggle_switch.setChecked(False)  # Default back to "NO"

    def on_toggle_state_changed(self, state):
        """
        Handles state changes of the toggle switch.
        Updates the global dictionary to reflect the current state.
        """
        parameters_dict["reset_updates"] = 1.0 if state else 0.0  # 1.0 for YES, 0.0 for NO
