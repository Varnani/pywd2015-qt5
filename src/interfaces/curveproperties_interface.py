from PyQt5 import QtWidgets, QtGui
import numpy
from gui import curveproperties_dialog
from src.helpers import methods, messenger
from src import constants
from matplotlib import pyplot
import os


TYPE_LIGHT = "light"
TYPE_VELOCITY = "velocity"


class CurveProperty:
    def __init__(self, curve_type, main_window):
        self.curve_type = curve_type

        self.main_window = main_window

        self.unpopulated = True
        self.constant_weights = False

        self.filepath = None
        self.band_id = None

        self.l1 = None
        self.l2 = None
        self.l3 = None

        self.x1 = None
        self.x2 = None
        self.y1 = None
        self.y2 = None

        self.wla = None
        self.aextinc = None
        self.xunit = None
        self.calib = None

        self.ksd = None
        self.noise = None
        self.sigma = None
        self.opsf = None

        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.e4 = None

        self.col1_idx = None
        self.col2_idx = None
        self.col3_idx = None

    def show_dialog(self):
        dialog = Widget(self.curve_type, self.main_window)
        if dialog.execute() == 1:
            self.unpopulated = False
            self.populate_from_ui(dialog)

    def edit(self):
        dialog = self.populate_ui(Widget(self.curve_type, self.main_window))
        if dialog.execute() == 1:
            self.populate_from_ui(dialog)

    def populate_ui(self, ui):

        ui.filepath = self.filepath
        ui.band_box.setValue(self.band_id)

        ui.l1_ipt.setValue(self.l1)
        ui.l2_ipt.setValue(self.l2)
        ui.l3_ipt.setValue(self.l3)

        ui.x1_ipt.setValue(self.x1)
        ui.x2_ipt.setValue(self.x2)
        ui.y1_ipt.setValue(self.y1)
        ui.y2_ipt.setValue(self.y2)

        ui.wla_ipt.setValue(self.wla)
        ui.aextinc_ipt.setValue(self.aextinc)
        ui.xunit_ipt.setValue(self.xunit)
        ui.calib_ipt.setValue(self.calib)

        ui.ksd_spinbox.setValue(self.ksd)
        ui.noise_combobox.setCurrentIndex(self.noise)
        ui.sigma_ipt.setValue(self.sigma)
        ui.opsf_ipt.setValue(self.opsf)

        ui.e1_ipt.setValue(self.e1)
        ui.e2_ipt.setValue(self.e2)
        ui.e3_ipt.setValue(self.e3)
        ui.e4_ipt.setValue(self.e4)

        ui.populate_data_from_filepath(self.filepath)

        ui.time_combobox.setCurrentIndex(self.col1_idx)
        ui.obs_combobox.setCurrentIndex(self.col2_idx)
        ui.weight_combobox.setCurrentIndex(self.col3_idx)

        return ui

    def populate_from_ui(self, ui):

        self.filepath = ui.filepath
        self.band_id = ui.band_box.value()

        self.l1 = ui.l1_ipt.value()
        self.l2 = ui.l2_ipt.value()
        self.l3 = ui.l3_ipt.value()

        self.x1 = ui.x1_ipt.value()
        self.x2 = ui.x2_ipt.value()
        self.y1 = ui.y1_ipt.value()
        self.y2 = ui.y2_ipt.value()

        self.wla = ui.wla_ipt.value()
        self.aextinc = ui.aextinc_ipt.value()
        self.xunit = ui.xunit_ipt.value()
        self.calib = ui.calib_ipt.value()

        self.ksd = ui.ksd_spinbox.value()
        self.noise = ui.noise_combobox.currentIndex()
        self.sigma = ui.sigma_ipt.value()
        self.opsf = ui.opsf_ipt.value()

        self.e1 = ui.e1_ipt.value()
        self.e2 = ui.e2_ipt.value()
        self.e3 = ui.e3_ipt.value()
        self.e4 = ui.e4_ipt.value()

        self.col1_idx = ui.time_combobox.currentIndex()
        self.col2_idx = ui.obs_combobox.currentIndex()
        self.col3_idx = ui.weight_combobox.currentIndex()

        if str(ui.weight_combobox.currentText()) == "Constant (1.0)":
            self.constant_weights = True
        else:
            self.constant_weights = False


    def populate_parser(self, parser, section):

        parser.set(section, "filepath", self.filepath)
        parser.set(section, "band_id", str(self.band_id))
        parser.set(section, "constant_weights", str(self.constant_weights))

        parser.set(section, "l1", str(self.l1))
        parser.set(section, "l2", str(self.l2))
        parser.set(section, "l3", str(self.l3))

        parser.set(section, "x1", str(self.x1))
        parser.set(section, "x2", str(self.x2))
        parser.set(section, "y1", str(self.y1))
        parser.set(section, "y2", str(self.y2))

        parser.set(section, "wla", str(self.wla))
        parser.set(section, "aextinc", str(self.aextinc))
        parser.set(section, "xunit", str(self.xunit))
        parser.set(section, "calib", str(self.calib))

        parser.set(section, "ksd", str(self.ksd))
        parser.set(section, "noise", str(self.noise))
        parser.set(section, "sigma", str(self.sigma))
        parser.set(section, "opsf", str(self.opsf))

        parser.set(section, "e1", str(self.e1))
        parser.set(section, "e2", str(self.e2))
        parser.set(section, "e3", str(self.e3))
        parser.set(section, "e4", str(self.e4))

        parser.set(section, "col1_idx", str(self.col1_idx))
        parser.set(section, "col2_idx", str(self.col2_idx))
        parser.set(section, "col3_idx", str(self.col3_idx))

        return parser

    def populate_from_parser(self, parser, section):
        self.unpopulated = False

        self.filepath = parser.get(section, "filepath")
        self.band_id = parser.getint(section, "band_id")
        self.constant_weights = parser.getboolean(section, "constant_weights")

        self.l1 = parser.getfloat(section, "l1")
        self.l2 = parser.getfloat(section, "l2")
        self.l3 = parser.getfloat(section, "l3")

        self.x1 = parser.getfloat(section, "x1")
        self.x2 = parser.getfloat(section, "x2")
        self.y1 = parser.getfloat(section, "y1")
        self.y2 = parser.getfloat(section, "y2")

        self.wla = parser.getfloat(section, "wla")
        self.aextinc = parser.getfloat(section, "aextinc")
        self.xunit = parser.getfloat(section, "xunit")
        self.calib = parser.getfloat(section, "calib")

        self.ksd = parser.getfloat(section, "ksd")
        self.noise = parser.getint(section, "noise")
        self.sigma = parser.getfloat(section, "sigma")
        self.opsf = parser.getfloat(section, "opsf")

        self.e1 = parser.getfloat(section, "e1")
        self.e2 = parser.getfloat(section, "e2")
        self.e3 = parser.getfloat(section, "e3")
        self.e4 = parser.getfloat(section, "e4")

        self.col1_idx = parser.getint(section, "col1_idx")
        self.col2_idx = parser.getint(section, "col2_idx")
        self.col3_idx = parser.getint(section, "col3_idx")

    def get_data(self):
        if self.filepath is not None:
            data = numpy.loadtxt(self.filepath, unpack=True)

            x = data[self.col1_idx]
            y = data[self.col2_idx]
            z = []

            if self.constant_weights:
                z = numpy.ones(len(x))
            else:
                z = data[self.col3_idx]

            return x, y, z


class Widget(QtWidgets.QDialog, curveproperties_dialog.Ui_CurvePropertiesDialog):
    def __init__(self, curve_type, main_window):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.main_window = main_window

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.data_widget.setFont(main_window.monoFont)
        self.data_widget.header().setSectionResizeMode(3)

        self.filepath = None

        self.old_bandbox_value = 7

        self.apply_constraints()

        if curve_type == TYPE_LIGHT:
            self.title_label.setText("Create or modify a light curve")

        elif curve_type == TYPE_VELOCITY:
            self.title_label.setText("Create or modify a velocity curve")

            self.band_box.setDisabled(True)
            self.band_name_label.setText(" - ")
            self.bandpasscontextlist_btn.setDisabled(True)

            self.l1_ipt.setDisabled(True)
            self.l1_ipt.setValue(1.0)

            self.l2_ipt.setDisabled(True)
            self.l2_ipt.setValue(1.0)

            self.l3_ipt.setDisabled(True)
            self.l3_ipt.setValue(0.0)

            self.x1_ipt.setDisabled(True)
            self.x2_ipt.setDisabled(True)
            self.y1_ipt.setDisabled(True)
            self.y2_ipt.setDisabled(True)

            self.opsf_ipt.setDisabled(True)
            self.noise_combobox.setDisabled(True)

            self.e1_ipt.setDisabled(False)
            self.e2_ipt.setDisabled(False)
            self.e3_ipt.setDisabled(False)
            self.e4_ipt.setDisabled(False)

            self.aextinc_ipt.setDisabled(True)
            self.xunit_ipt.setDisabled(True)
            self.calib_ipt.setDisabled(True)

        self.connect_signals()

    def connect_signals(self):
        self.accept_btn.clicked.connect(self.accept_changes)
        self.discard_btn.clicked.connect(self.discard_changes)
        self.bandpasscontextlist_btn.clicked.connect(self.bandpass_menu)
        self.band_box.valueChanged.connect(self.band_box_value_change)
        self.whatsthis_btn.clicked.connect(QtWidgets.QWhatsThis.enterWhatsThisMode)
        self.plot_btn.clicked.connect(self.plot_preview_data)
        self.repick_btn.clicked.connect(self.repick_data)

    def accept_changes(self):
        self.done(1)

    def discard_changes(self):
        self.done(0)

    def repick_data(self):
        filepath = methods.load_file(self, name_filter="Observation data (*)")
        if filepath is not None and filepath != self.filepath:
            if self.populate_data_from_filepath(filepath):
                self.filepath = filepath
                self.time_combobox.setCurrentIndex(0)
                self.obs_combobox.setCurrentIndex(1)
                self.weight_combobox.setCurrentIndex(2)

    def plot_preview_data(self):
        data = numpy.loadtxt(self.filepath, unpack=True)

        x = data[self.time_combobox.currentIndex()]
        y = data[self.obs_combobox.currentIndex()]

        pyplot.plot(x, y, linestyle="", marker="o", markersize=2, color=constants.COLOR_BLUE)
        pyplot.get_current_fig_manager().set_window_title("Matplotlib - " + os.path.basename(self.filepath))
        pyplot.show()

    def bandpass_menu(self):
        selection = methods.create_bandpass_menu(self).exec_(QtGui.QCursor.pos())

        if selection is not None:
            band_id = int(constants.BANDPASS_ID_DICT[selection.objectName()])
            self.band_box.setValue(band_id)

    def band_box_value_change(self):
        try:
            self.band_name_label.setText(constants.ID_BANDPASS_DICT[str(self.band_box.value())])
            self.old_bandbox_value = self.band_box.value()
        except KeyError:
            self.band_box.setValue(self.old_bandbox_value)

    def populate_data_from_filepath(self, filepath):

        try:
            data = numpy.loadtxt(filepath, unpack=True)
            if 1 > len(data):
                raise IndexError("File has less than 2 columns.")

        except ValueError as e:
            msg = messenger.Messenger("error", "A ValueError has occured:")
            msg.set_info(e.args[0])
            msg.show()
            return False

        except IndexError as e:
            msg = messenger.Messenger("error", "An IndexError has occured:")
            msg.set_info(e.args[0])
            msg.show()
            return False

        self.clear_widgets()

        self.filepath_label.setText(filepath)
        self.filepath_label.setToolTip(filepath)

        header_item = QtWidgets.QTreeWidgetItem()
        for col_idx in range(len(data)):
            col_name = "Column " + str(col_idx)
            header_item.setText(col_idx, col_name)
            self.time_combobox.addItem(col_name)
            self.obs_combobox.addItem(col_name)
            self.weight_combobox.addItem(col_name)
        self.data_widget.setHeaderItem(header_item)

        self.weight_combobox.addItem("Constant (1.0)")

        for row_idx in range(len(data[0])):
            item = QtWidgets.QTreeWidgetItem(self.data_widget)
            for col_idx in range(len(data)):
                item.setText(col_idx, str(data[col_idx][row_idx]))

        self.data_widget.header().setSectionResizeMode(3)

        return True

    def execute(self):
        if self.filepath is None:
            self.filepath = methods.load_file(self, suffix="", name_filter="Observation Data (*)")
            if self.filepath is None:
                return 0

            if self.populate_data_from_filepath(self.filepath):
                self.time_combobox.setCurrentIndex(0)
                self.obs_combobox.setCurrentIndex(1)
                self.weight_combobox.setCurrentIndex(2)
            else:
                return 0

        return self.exec_()

    def apply_constraints(self):
        if self.main_window.ld1_chk.isChecked() is not True:
            self.x1_ipt.setDisabled(True)
            self.y1_ipt.setDisabled(True)
        if self.main_window.ld2_chk.isChecked() is not True:
            self.x2_ipt.setDisabled(True)
            self.y2_ipt.setDisabled(True)
        if self.main_window.ipb_chk.isChecked() is not True:
            self.l2_ipt.setDisabled(True)
            if self.l2_ipt.value() == 0.0:
                self.l2_ipt.setValue(1.0)
        if str(self.main_window.mode_combobox.currentText()) == "Mode 1":
            if self.main_window.ld1_chk.isChecked() and self.main_window.ld2_chk.isChecked():
                self.x2_ipt.setDisabled(True)
                self.y2_ipt.setDisabled(True)

                #def _lockx2y2():
                    #self.x2_ipt.setText(self.x1_ipt.text())
                    #self.y2_ipt.setText(self.y1_ipt.text())

                #self.x1_ipt.valueChanged.connect(_lockx2y2)
                #self.y1_ipt.valueChanged.connect(_lockx2y2)

    def clear_widgets(self):
        self.data_widget.clear()
        self.time_combobox.clear()
        self.obs_combobox.clear()
        self.weight_combobox.clear()
