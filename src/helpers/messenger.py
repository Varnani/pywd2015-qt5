from PyQt5 import QtWidgets, QtGui
from src import constants


class Messenger:
    def __init__(self, icon_type, msg):
        self.msg_box = QtWidgets.QMessageBox()
        self.msg_box.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))
        self.msg_box.setText(msg)
        self.msg_box.setWindowTitle("PyWD2015")
        if icon_type == "info":
            self.msg_box.setIcon(QtWidgets.QMessageBox.Information)
        elif icon_type == "warning":
            self.msg_box.setIcon(QtWidgets.QMessageBox.Warning)
        elif icon_type == "question":
            self.msg_box.setIcon(QtWidgets.QMessageBox.Question)
        elif icon_type == "error":
            self.msg_box.setIcon(QtWidgets.QMessageBox.Critical)

    def show(self):
        return self.msg_box.exec_()

    def set_info(self, msg):
        self.msg_box.setInformativeText(msg)
