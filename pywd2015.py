import platform
import sys
from PyQt5 import QtWidgets
from src.interfaces import mainwindow_interface


def run():
    if platform.system() == "Windows":  # used for making app icon visible in windows taskbar
        import ctypes
        appid = "pywd2015-qt5"
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(appid)
    app = QtWidgets.QApplication(sys.argv)
    gui = mainwindow_interface.Widget(app)
    if gui.start_up():
        gui.show()
        sys.exit(app.exec_())
    else:
        return 0


if __name__ == "__main__":
    run()
