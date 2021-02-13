from PyQt5 import QtWidgets, QtGui
from gui import oc_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src import constants
import os
from src.helpers import messenger
from src.helpers.wd_utils import wd_io
import numpy
from src.helpers import methods


delimiter = constants.EXPORT_DELIMITER


class Widget(QtWidgets.QWidget, oc_widget.Ui_OCWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.chart = MatplotlibWidget(self.plot_widget, 1, 1, exportable=False)
        self.chart.create_axis(0, 0, name="", labels=("E", "O - C (Day)"))

        self.data_treewidget.setFont(parent.monoFont)
        self.data_treewidget.header().setSectionResizeMode(3)

        self.cycles = []
        self.linear_resid = []
        self.dpdt_resid = []

        self.connect_signals()

    def connect_signals(self):
        self.compute_btn.clicked.connect(self.compute_oc)
        self.dpdt_chk.toggled.connect(self.update_plots)
        self.linear_chk.toggled.connect(self.update_plots)
        self.calculate_btn.clicked.connect(self.calculate_dp_and_dt)
        self.update_btn.clicked.connect(self.update_ephemeris_and_period)
        self.export_btn.clicked.connect(self.export_data)

    def export_data(self):
        root = self.data_treewidget.invisibleRootItem()
        if root.childCount() == 0:
            return 1

        filepath = methods.save_file(self, suffix="txt")
        if filepath is not None:
            lenght = 3
            output = "# HJD" + delimiter + "Linear Residuals" + delimiter + "Residuals with dP/dt\n"

            root = self.data_treewidget.invisibleRootItem()
            i = 0
            while i < root.childCount():
                idx = 0
                while idx < lenght:
                    output = output + root.child(i).text(idx) + delimiter
                    idx = idx + 1

                output = output + "\n"
                i = i + 1

            with open(filepath, "w") as destination:
                destination.write(output)

            msg = messenger.Messenger("info", "File saved:")
            msg.set_info(filepath)
            msg.show()

    def compute_oc(self):
        eclipse_timings_path = self.main_window.eclipsetimings_widget.filepath_label.text()
        if os.path.isfile(eclipse_timings_path):

            self.clear()

            lc_params = self.main_window.get_lc_params()

            ec_data = numpy.loadtxt(eclipse_timings_path, unpack=True)
            lc_params.add_eclipse_times(ec_data[0], ec_data[1])

            lc_io = wd_io.LCIO(lc_params,
                               wd_path=self.main_window.lc_path,
                               lc_binary_name=self.main_window.lc_binary)

            results = lc_io.fill_for_etv().save().run().read_etv()

            i = 0
            while i < len(results[0]):
                jd = results[0][i]
                lin_res = results[3][i]
                dpdt_res = results[5][i]

                item = QtWidgets.QTreeWidgetItem(self.data_treewidget)
                item.setText(0, str(jd))
                item.setText(1, str(lin_res))
                item.setText(2, str(dpdt_res))

                i = i + 1

            t0 = self.main_window.jd0_ipt.value()
            p = self.main_window.p0_ipt.value()

            for t in results[0]:
                cycle = (t - t0) / p
                e = numpy.around(cycle * 2.0, decimals=0) / 2.0
                self.cycles.append(e)

            self.linear_resid = results[3]
            self.dpdt_resid = results[5]

            self.update_plots()

        else:
            msg = messenger.Messenger("error", "An eclipse timings file must be provided for O - C calculation.")
            msg.set_info("You can load eclipse timings from the main tab.")
            msg.show()

    def update_plots(self):
        self.chart.clear_all()

        if self.linear_chk.isChecked():
            self.chart.plot(self.cycles, self.linear_resid, clear=False,
                            markersize=2, marker="o", linestyle="", color=constants.COLOR_BLUE)

        if self.dpdt_chk.isChecked():
            self.chart.plot(self.cycles, self.dpdt_resid, clear=False,
                            markersize=2, marker="o", linestyle="", color=constants.COLOR_RED)

    def calculate_dp_and_dt(self):
        if len(self.cycles) != 0 and self.linear_chk.isChecked():

            model_coeffs = numpy.polyfit(self.cycles, self.linear_resid, 1)
            model_y = [(y * model_coeffs[0] + model_coeffs[1]) for y in self.cycles]

            self.dp_otpt.setValue(model_coeffs[0])
            self.dt_otpt.setValue(model_coeffs[1])

            self.chart.plot(self.cycles, model_y, clear=False, color=constants.COLOR_RED)

    def update_ephemeris_and_period(self):
        if len(self.cycles) != 0:
            self.main_window.jd0_ipt.setValue(self.main_window.jd0_ipt.value() + self.dt_otpt.value())
            self.main_window.p0_ipt.setValue(self.main_window.p0_ipt.value() + self.dp_otpt.value())

            self.update_plots()

            self.dp_otpt.setValue(0.0)
            self.dt_otpt.setValue(0.0)

            self.compute_btn.click()

    def clear(self):
        self.cycles = []
        self.linear_resid = []
        self.dpdt_resid = []

        self.chart.clear_all()
        self.data_treewidget.clear()
        self.dp_otpt.setValue(0.0)
        self.dt_otpt.setValue(0.0)
