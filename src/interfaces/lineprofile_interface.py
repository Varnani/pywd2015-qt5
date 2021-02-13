from PyQt5 import QtWidgets, QtGui
from gui import lineprofile_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src import constants
from src.helpers.wd_utils import wd_io


class Widget(QtWidgets.QWidget, lineprofile_widget.Ui_LineProfileWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.chart = MatplotlibWidget(self.plot_widget, 1, 1)
        self.chart.create_axis(0, 0, name="Spectral Line Profiles", labels=("Microns", "Norm. Flux"))

        self.s1_treewidget.header().setSectionResizeMode(1)
        self.s1_treewidget.header().setSectionsMovable(False)
        self.s2_treewidget.header().setSectionResizeMode(1)
        self.s2_treewidget.header().setSectionsMovable(False)

        self.connect_signals()

    def connect_signals(self):
        self.plot_btn.clicked.connect(self.plot_lines)
        self.s1_add_btn.clicked.connect(self.add_row_to_primary)
        self.s2_add_btn.clicked.connect(self.add_row_to_secondary)
        self.s1_remove_btn.clicked.connect(self.remove_row_from_primary)
        self.s2_remove_btn.clicked.connect(self.remove_row_from_secondary)

    def add_row_to_primary(self):
        self.add_row(self.s1_treewidget)

    def add_row_to_secondary(self):
        self.add_row(self.s2_treewidget)

    def remove_row_from_primary(self):
        self.remove_row(self.s1_treewidget)

    def remove_row_from_secondary(self):
        self.remove_row(self.s2_treewidget)

    def get_selected(self, widget):
        selected_items = widget.selectedItems()
        if len(selected_items) > 0:
            return selected_items[0]
        else:
            return None

    def add_row(self, widget):
        item = QtWidgets.QTreeWidgetItem(widget)
        i = 0
        values = [0.65627, 0.00001, 0.5]
        while i < 3:
            spinbox = QtWidgets.QDoubleSpinBox(widget)
            spinbox.setButtonSymbols(2)
            spinbox.setDecimals(7)
            spinbox.setMaximum(999)
            spinbox.setMinimum(0)
            spinbox.setValue(values[i])
            widget.setItemWidget(item, i, spinbox)
            i = i + 1

        integer_spinbox = QtWidgets.QSpinBox(widget)
        integer_spinbox.setButtonSymbols(2)
        integer_spinbox.setMinimum(-999)
        integer_spinbox.setMaximum(999)
        integer_spinbox.setValue(0)
        widget.setItemWidget(item, 3, integer_spinbox)

    def remove_row(self, widget):
        selected_item = self.get_selected(widget)
        if selected_item is not None:
            widget.headerItem().removeChild(selected_item)

    def plot_lines(self):
        self.chart.clear_all()

        lc_params = self.main_window.get_lc_params()

        lc_params["binwm1"] = self.s1_binwidth_spinbox.value()
        lc_params["sc1"] = self.s1_contscale_spinbox.value()
        lc_params["sl1"] = self.s1_contslope_spinbox.value()
        lc_params["nf1"] = self.s1_subgrid_spinbox.value()

        lc_params["binwm2"] = self.s2_binwidth_spinbox.value()
        lc_params["sc2"] = self.s2_contscale_spinbox.value()
        lc_params["sl2"] = self.s2_contslope_spinbox.value()
        lc_params["nf2"] = self.s2_subgrid_spinbox.value()

        # def add_spectral_line(self, star, wll, ewid, depth, kks):

        s1_child_count = self.s1_treewidget.invisibleRootItem().childCount()
        index = 0
        while index < s1_child_count:
            item = self.s1_treewidget.invisibleRootItem().child(index)
            lc_params.add_spectral_line(
                1,
                self.s1_treewidget.itemWidget(item, 0).value(),
                self.s1_treewidget.itemWidget(item, 1).value(),
                self.s1_treewidget.itemWidget(item, 2).value(),
                self.s1_treewidget.itemWidget(item, 3).value()
            )

            index = index + 1

        s2_child_count = self.s2_treewidget.invisibleRootItem().childCount()
        index = 0
        while index < s2_child_count:
            item = self.s2_treewidget.invisibleRootItem().child(index)
            lc_params.add_spectral_line(
                2,
                self.s2_treewidget.itemWidget(item, 0).value(),
                self.s2_treewidget.itemWidget(item, 1).value(),
                self.s2_treewidget.itemWidget(item, 2).value(),
                self.s2_treewidget.itemWidget(item, 3).value()
            )

            index = index + 1

        lc_params["jdphs"] = 2
        lc_params["phstrt"] = self.phase_spinbox.value()
        lc_params["phstop"] = self.phase_spinbox.value()
        lc_params["phin"] = 0.1

        lc_io = wd_io.LCIO(lc_params,
                           wd_path=self.main_window.lc_path,
                           lc_binary_name=self.main_window.lc_binary)

        s1_results, s2_results = lc_io.fill_for_spectral_lines().save().run().read_spectral_lines()

        if len(s1_results[0]) > 0:
            self.chart.plot(s1_results[0][2], s1_results[0][4], color=constants.COLOR_BLUE)

        if len(s2_results[0]) > 0:
            self.chart.plot(s2_results[0][2], s2_results[0][4], clear=False, color=constants.COLOR_RED)

    def clear(self):
        self.chart.clear_all()
