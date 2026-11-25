from PyQt5 import QtWidgets, QtGui
from gui import eclipsetimings_widget
from src.helpers import methods, messenger
from src.helpers.observations import read_timings
import numpy
from src import constants


class Widget(QtWidgets.QWidget, eclipsetimings_widget.Ui_EclipseWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent
        self.datawidget.setFont(parent.monoFont)
        self.datawidget.header().setSectionResizeMode(3)
        self.connect_signals()

    def connect_signals(self):
        self.load_btn.clicked.connect(self.load_file)
        self.clear_btn.clicked.connect(self.clear_timings)

    def load_file(self):
        filepath = methods.load_file(self)
        if filepath is not None:
            self.load_timings(filepath)

    def get_data(self):
        return read_timings(self.filepath_label.text(), self.constant_weight_checkbox.isChecked())

    def load_timings(self, filepath):
        try:
            data = read_timings(filepath, self.constant_weight_checkbox.isChecked())
        except (ValueError, OSError) as error:
            msg = messenger.Messenger("error", "Eclipse times could not be loaded:")
            msg.set_info(str(filepath) + "\n" + str(error))
            msg.show()
            return False
        self.datawidget.clear()
        for row in data.T:
            item = QtWidgets.QTreeWidgetItem(self.datawidget)
            for column, value in enumerate(row):
                item.setText(column, str(value))
        self.filepath_label.setText(filepath)
        return True

    def clear_timings(self):
        self.datawidget.clear()
        self.ksd_box.setValue(1)
        self.sigma_ipt.setValue(0.0)
        self.iftime_chk.setChecked(False)
        self.constant_weight_checkbox.setChecked(False)
        self.filepath_label.setText("None")

    def write_into_parser(self, parser):
        parser.add_section(constants.CONFIG_SECTION_ECLIPSE_TIMINGS)
        parser.set(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "iftime", str(self.iftime_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "cw", str(self.constant_weight_checkbox.isChecked()))
        parser.set(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "ksd", str(self.ksd_box.value()))
        parser.set(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "sigma", str(self.sigma_ipt.value()))
        parser.set(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "filepath", str(self.filepath_label.text()))

        return parser

    def read_from_parser(self, parser):
        self.clear_timings()

        self.iftime_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "iftime"))
        self.constant_weight_checkbox.setChecked(parser.getboolean(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "cw"))
        self.ksd_box.setValue(parser.getint(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "ksd"))
        self.sigma_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "sigma"))
        path = parser.get(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "filepath")

        if path != "None":
            self.load_timings(path)
