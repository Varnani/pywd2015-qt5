from PyQt5 import QtWidgets, QtGui
from gui import conjunction_widget
from src import constants
from src.helpers.wd_utils import wd_io
from src.helpers import methods
from src.helpers import messenger
import numpy


delimiter = constants.EXPORT_DELIMITER


class Widget(QtWidgets.QWidget, conjunction_widget.Ui_conjunctionwidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.data_treewidget.setFont(parent.monoFont)
        self.data_treewidget.header().setSectionResizeMode(3)

        self.connect_signals()

    def connect_signals(self):
        self.compute_btn.clicked.connect(self.compute_conjunction)
        self.export_btn.clicked.connect(self.export_data)

    def export_data(self):
        root = self.data_treewidget.invisibleRootItem()
        if root.childCount() == 0:
            return 1

        filepath = methods.save_file(self, suffix="txt")
        if filepath is not None:
            lenght = None
            output = ""

            if self.ut_groupbox.isChecked():
                lenght = 3
                if self.dt_groupbox.isChecked():
                    output = "# JD" + delimiter + "MinType" + delimiter + "Time\n"
                else:
                    output = "# HJD" + delimiter + "MinType" + delimiter + "Time\n"

            else:
                lenght = 2
                output = "# HJD" + delimiter + "MinType\n"

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

    def compute_conjunction(self):
        self.data_treewidget.clear()

        lc_params = self.main_window.get_lc_params()

        lc_io = wd_io.LCIO(lc_params,
                           wd_path=self.main_window.lc_path,
                           lc_binary_name=self.main_window.lc_binary)

        results = lc_io.fill_for_conjunction(self.kstep_spinbox.value()).save().run().read_conjunction()

        i = 0
        while i < len(results[0]):
            jd = results[0][i]
            mn = results[1][i]

            item = QtWidgets.QTreeWidgetItem(self.data_treewidget)

            if self.ut_groupbox.isChecked():
                if self.dt_groupbox.isChecked():
                    jd = methods.convert_hjd_to_jd(jd,
                                                   self.ra_h_spinbox.value(),
                                                   self.ra_m_spinbox.value(),
                                                   self.ra_s_spinbox.value(),
                                                   self.dec_d_spinbox.value(),
                                                   self.dec_m_spinbox.value(),
                                                   self.dec_s_spinbox.value())

                year, month, day, hour, minute, second = methods.convert_jd_to_ut(jd)
                item.setText(2, str(day) + "/" + str(month) + "/" + str(year) + " - " +
                             str(hour) + "/" + str(minute) + "/" + str(numpy.round(second, decimals=3)))

            item.setText(0, str(numpy.round(jd, decimals=5)))
            item.setText(1, str(int(mn)))
            i = i + 1

        header_item = QtWidgets.QTreeWidgetItem()
        if self.ut_groupbox.isChecked() and self.dt_groupbox.isChecked():
            header_item.setText(0, "JD")
        else:
            header_item.setText(0, "HJD")
        header_item.setText(1, "Min")
        header_item.setText(2, "D/M/Y - H:M:S")
        self.data_treewidget.setHeaderItem(header_item)
