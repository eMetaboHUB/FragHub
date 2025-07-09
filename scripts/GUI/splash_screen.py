from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSpacerItem, QSizePolicy, QApplication
)
from PyQt6.QtGui import QPainter, QColor, QPen, QPixmap
from PyQt6.QtCore import Qt, QTimer, QRectF

class InfiniteSpinnerWidget(QWidget):
    """Widget pour une animation de chargement infinie."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self._angle = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._update_angle)
        self._timer.setInterval(25)
        self._spinner_color = QColor(Qt.GlobalColor.white)
        self._pen_width = 3.0
        self.setFixedSize(40, 40)

    def _update_angle(self):
        self._angle = (self._angle - 10) % 360
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        padding = self._pen_width / 2.0
        rect = QRectF(padding, padding, float(self.width()) - self._pen_width, float(self.height()) - self._pen_width)
        if rect.width() <= 0 or rect.height() <= 0:
            return
        pen = QPen(self._spinner_color, self._pen_width, Qt.PenStyle.SolidLine)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawArc(rect, self._angle * 16, 90 * 16)

    def start_animation(self):
        if not self._timer.isActive():
            self._timer.start()

    def stop_animation(self):
        if self._timer.isActive():
            self._timer.stop()

    def set_color(self, color: QColor):
        self._spinner_color = color
        self.update()

    def hideEvent(self, event):
        self.stop_animation()
        super().hideEvent(event)

    def showEvent(self, event):
        self.start_animation()
        super().showEvent(event)


class LoadingSplashScreen(QWidget):
    """Écran de chargement personnalisé avec un spinner."""
    def __init__(self, icon_pixmap: QPixmap, parent=None):
        super().__init__(parent)
        self.setWindowFlags(
            Qt.WindowType.SplashScreen | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        self.icon_label = QLabel()
        if not icon_pixmap.isNull():
            self.icon_label.setPixmap(icon_pixmap)
        self.icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.icon_label)
        layout.addSpacerItem(QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed))
        self.spinner = InfiniteSpinnerWidget(self)
        spinner_layout = QHBoxLayout()
        spinner_layout.addStretch()
        spinner_layout.addWidget(self.spinner)
        spinner_layout.addStretch()
        layout.addLayout(spinner_layout)
        layout.addSpacerItem(QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed))
        self.message_label = QLabel("Loading FragHub, please wait...")
        self.message_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.message_label.setStyleSheet("QLabel { color : white; font-weight: bold; }")
        layout.addWidget(self.message_label)
        self.setLayout(layout)
        self.adjustSize()
        if QApplication.primaryScreen():
            center_point = QApplication.primaryScreen().availableGeometry().center()
            self.move(center_point.x() - self.width() // 2, center_point.y() - self.height() // 2)

    def showMessage(self, message, font_size=28, alignment=Qt.AlignmentFlag.AlignCenter, color=Qt.GlobalColor.white):
        self.message_label.setText(message)
        self.message_label.setAlignment(alignment)
        self.message_label.setStyleSheet(f"QLabel {{ color: white; font-weight: bold; font-size: {font_size}px; }}")

    def closeEvent(self, event):
        self.spinner.stop_animation()
        super().closeEvent(event)