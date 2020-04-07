from PyQt5 import QtWidgets, QtGui
from gui import dimensions_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src import constants
from src.helpers.wd_utils import wd_io


class Widget(QtWidgets.QWidget, dimensions_widget.Ui_DimensionWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.s1_chart = MatplotlibWidget(self.s1_plot_widget, 1, 1)
        self.s2_chart = MatplotlibWidget(self.s2_plot_widget, 1, 1)
        self.s1_chart.create_axis(0, 0, name="Primary Star Dimensions", labels=("", "Fract. Radius"))
        self.s2_chart.create_axis(0, 0, name="Secondary Star Dimensions", labels=("", "Fract. Radius"))

        self.s1_pole = []
        self.s1_point = []
        self.s1_side = []
        self.s1_back = []

        self.s2_pole = []
        self.s2_point = []
        self.s2_side = []
        self.s2_back = []

        self.x = []

        self.connect_signals()

    def connect_signals(self):
        self.plot_btn.clicked.connect(self.calculate_and_plot)

        self.s1_pole_chk.toggled.connect(self.update_plots)
        self.s1_point_chk.toggled.connect(self.update_plots)
        self.s1_side_chk.toggled.connect(self.update_plots)
        self.s1_back_chk.toggled.connect(self.update_plots)

        self.s2_pole_chk.toggled.connect(self.update_plots)
        self.s2_point_chk.toggled.connect(self.update_plots)
        self.s2_side_chk.toggled.connect(self.update_plots)
        self.s2_back_chk.toggled.connect(self.update_plots)

    def calculate_and_plot(self):
        self.calculate_component_radii()
        self.update_plots()

    def update_plots(self):
        self.s1_chart.clear_all()
        self.s2_chart.clear_all()

        if self.s1_pole_chk.isChecked():
            self.s1_chart.plot(self.x, self.s1_pole, clear=False, color=constants.COLOR_BLUE)
        if self.s1_point_chk.isChecked():
            self.s1_chart.plot(self.x, self.s1_point, clear=False, color="black")
        if self.s1_side_chk.isChecked():
            self.s1_chart.plot(self.x, self.s1_side, clear=False, color=constants.COLOR_RED)
        if self.s1_back_chk.isChecked():
            self.s1_chart.plot(self.x, self.s1_back, clear=False, color=constants.COLOR_GREEN)

        if self.s2_pole_chk.isChecked():
            self.s2_chart.plot(self.x, self.s2_pole, clear=False, color=constants.COLOR_BLUE)
        if self.s2_point_chk.isChecked():
            self.s2_chart.plot(self.x, self.s2_point, clear=False, color="black")
        if self.s2_side_chk.isChecked():
            self.s2_chart.plot(self.x, self.s2_side, clear=False, color=constants.COLOR_RED)
        if self.s2_back_chk.isChecked():
            self.s2_chart.plot(self.x, self.s2_back, clear=False, color=constants.COLOR_GREEN)

    def calculate_component_radii(self):
        lc_params = self.main_window.get_lc_params()

        lc_io = wd_io.LCIO(lc_params,
                           wd_path=self.main_window.lc_path,
                           lc_binary_name=self.main_window.lc_binary)

        results = lc_io.fill_for_component_dimensions().save().run().read_component_dimensions()

        x_index = None
        if self.main_window.jd_radiobtn.isChecked():
            x_index = 0
            self.s1_chart.set_labels("HJD", "Fract. Radius")
            self.s2_chart.set_labels("HJD", "Fract. Radius")

        elif self.main_window.phase_radiobtn.isChecked():
            x_index = 1
            self.s1_chart.set_labels("Phase", "Fract. Radius")
            self.s2_chart.set_labels("Phase", "Fract. Radius")

        self.x = results[x_index]

        self.s1_pole = results[2]
        self.s1_point = results[3]
        self.s1_side = results[4]
        self.s1_back = results[5]

        self.s2_pole = results[6]
        self.s2_point = results[7]
        self.s2_side = results[8]
        self.s2_back = results[9]

        self.update_plots()

    def clear(self):
        self.s1_chart.clear_all()
        self.s2_chart.clear_all()

        self.s1_pole = []
        self.s1_point = []
        self.s1_side = []
        self.s1_back = []

        self.s2_pole = []
        self.s2_point = []
        self.s2_side = []
        self.s2_back = []

        self.x = []
