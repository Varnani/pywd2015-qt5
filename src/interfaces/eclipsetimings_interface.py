from PyQt5 import QtWidgets, QtGui
from gui import eclipsetimings_widget
from src.helpers import methods, messenger
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
        data = numpy.loadtxt(self.filepath_label.text(), unpack=True)

        if len(data) == 2:
            _ = list(data)
            _.append([])
            data = numpy.array(_)

        if self.constant_weight_checkbox.isChecked():
            data[2] = numpy.ones(len(data[1]))

        return data

    def load_timings(self, filepath):
        try:
            data = numpy.loadtxt(filepath, unpack=True)
            if 1 > len(data):
                raise IndexError("File has less than 2 columns.")

        except ValueError as e:
            msg = messenger.Messenger("error", "A ValueError has occured:")
            msg.set_info(e.args[0] + "\nEclipse times are not loaded.")
            msg.show()
            return False

        except IndexError as e:
            msg = messenger.Messenger("error", "An IndexError has occured:")
            msg.set_info(e.args[0] + "\nEclipse times are not loaded.")
            msg.show()
            return False

        self.datawidget.clear()

        for row_idx in range(len(data[0])):
            item = QtWidgets.QTreeWidgetItem(self.datawidget)
            for col_idx in range(len(data)):
                item.setText(col_idx, str(data[col_idx][row_idx]))

        self.filepath_label.setText(filepath)

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
