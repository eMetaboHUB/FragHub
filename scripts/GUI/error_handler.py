import sys
import traceback
from PyQt6.QtWidgets import QMessageBox, QTextEdit, QApplication

def show_error_message(parent, title, message, on_close=None):
    """Affiche une erreur dans une QMessageBox élargie avec un QTextEdit."""
    msg_box = QMessageBox(parent)
    msg_box.setIcon(QMessageBox.Icon.Critical)
    msg_box.setWindowTitle(title)
    msg_box.setText("An error occurred during execution.")
    msg_box.setStandardButtons(QMessageBox.StandardButton.Ok)

    text_area = QTextEdit()
    text_area.setText(message)
    text_area.setReadOnly(True)
    text_area.setStyleSheet("font-family: Consolas, monospace; font-size: 10pt;")
    text_area.setMinimumSize(700, 350)

    try:
        grid_layout = msg_box.layout()
        grid_layout.addWidget(text_area, grid_layout.rowCount(), 0, 1, grid_layout.columnCount())
        msg_box.setMinimumSize(750, 450)
    except Exception as layout_e:
        print(f"Warning: Could not add QTextEdit to QMessageBox layout ({layout_e}). Using setDetailedText as fallback.")
        msg_box.setDetailedText(message)

    clicked_button = msg_box.exec()

    if clicked_button == QMessageBox.StandardButton.Ok and on_close:
        try:
            on_close()
        except Exception as callback_e:
            print(f"ERREUR: Exception dans le callback on_close de show_error_message: {callback_e}")
            traceback.print_exc()

def exception_hook(exctype, value, tb):
    """Hook global pour capturer les exceptions non gérées."""
    error_message = ''.join(traceback.format_exception(exctype, value, tb))
    sys.stderr.write(f"Uncaught exception:\n{error_message}")
    if QApplication.instance():
        # Utilise la fonction locale pour afficher le message
        show_error_message(
            parent=None,
            title="Unhandled Exception",
            message=error_message
        )
    sys.exit(1)