from PyQt5 import QtWidgets, QtGui
from gui import lcdcpicker_dialog
from src import constants
from src.helpers import methods, messenger
from functools import partial
from configparser import ConfigParser, NoOptionError
import os


class Widget(QtWidgets.QDialog, lcdcpicker_dialog.Ui_LCDCPickerDialog):
    def __init__(self):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        def _pick(label):
            path = methods.load_file(self)
            if path is not None:
                label.setText(path)
                label.setToolTip(path)

        self.picklc_btn.clicked.connect(partial(_pick, self.lcpath_label))
        self.pickdc_btn.clicked.connect(partial(_pick, self.dcpath_label))

        self.save_btn.clicked.connect(self.save)
        self.exit_btn.clicked.connect(self.exit)

    def check_config(self):
        if os.path.isfile(constants.PYWD_CONFIG_PATH):
            try:
                config = ConfigParser()
                config.read(constants.PYWD_CONFIG_PATH)
                lc_path = config.get("wd_paths", "lc_path")
                dc_path = config.get("wd_paths", "dc_path")

                if os.path.isfile(lc_path) and os.path.isfile(dc_path):
                    self.lcpath_label.setText(lc_path)
                    self.dcpath_label.setText(dc_path)
                    return True

                else:
                    msg = messenger.Messenger("warning", "Paths from the existing pywd_config file are not valid.")
                    msg.set_info("Please reselect wd binary paths.")
                    msg.show()
                    return self.exec_()
            except NoOptionError as ex:
                msg = messenger.Messenger("error", "A NoOptionError occured during config file read.")
                msg.set_info("Config file is malformed. Please reselect wd binary paths.\nAdditional info:\n" +
                             ex.args[0])
                msg.show()
                return self.exec_()

        else:
            return self.exec_()

    def save(self):
        error = ""

        if os.path.isfile(self.lcpath_label.text()) is not True:
            error = "Can't find LC file in path:\n" + self.lcpath_label.text()

        if os.path.isfile(self.dcpath_label.text()) is not True:
            error = error + "\nCan't DC find file in path:\n" + self.dcpath_label.text()

        if error == "":
            config = ConfigParser()
            config.add_section("wd_paths")
            config.set("wd_paths", "lc_path", self.lcpath_label.text())
            config.set("wd_paths", "dc_path", self.dcpath_label.text())
            with open(constants.PYWD_CONFIG_PATH, "w") as f:
                config.write(f)

            self.accept()

        else:
            msg = messenger.Messenger("error", "Error(s) occured:")
            msg.set_info(error)
            msg.show()

    def exit(self):
        self.reject()
