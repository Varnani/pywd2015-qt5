from PyQt5 import QtWidgets, QtGui
from gui import syntheticcurve_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src.helpers import methods
from functools import partial
from src import constants
import os
from src.helpers.wd_utils import wd_io
from src.helpers import messenger


class Widget(QtWidgets.QWidget, syntheticcurve_widget.Ui_SyntheticCurveWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent
        self.obs_widget = self.main_window.loadobservations_widget

        self.nlc = 0

        # setup light_chart area
        self.light_chart = MatplotlibWidget(self.light_plot_widget, 2, 1, h_ratios=[1.5, 1])
        self.light_chart.create_axis(0, 0, name="Synthetic Light Curve")
        self.light_chart.create_axis(1, 0, sharex=self.light_chart[0])
        self.light_chart.figure.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95, hspace=0, wspace=0)
        self.light_chart[0].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True,
                                        labelright=False, labelleft=True, labelsize=11)

        # setup velocity_chart area
        self.velocity_chart = MatplotlibWidget(self.vel_plot_widget, 2, 1, h_ratios=[1.5, 1])
        self.velocity_chart.create_axis(0, 0, name="Synthetic Velocity Curve")
        self.velocity_chart.create_axis(1, 0, sharex=self.velocity_chart[0])
        self.velocity_chart.figure.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95, hspace=0, wspace=0)
        self.velocity_chart[0].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True,
                                           labelright=False, labelleft=True, labelsize=11)

        self.light_treewidget.header().setSectionResizeMode(3)
        self.light_treewidget.header().setSectionsMovable(False)

        self.insert_synthetic_curves()

        self.connect_signals()

    def connect_signals(self):
        self.light_plot_btn.clicked.connect(self.plot_selected_light_curve)
        self.vel_plot_btn.clicked.connect(self.plot_selected_velocity_curve)

    def add_div_to_light(self):
        item = QtWidgets.QTreeWidgetItem(self.light_treewidget)
        i = 0
        while i < 15:
            separator = QtWidgets.QFrame(self.light_treewidget)
            separator.setFrameShape(QtWidgets.QFrame.HLine)
            separator.setFrameShadow(QtWidgets.QFrame.Sunken)
            self.light_treewidget.setItemWidget(item, i, separator)
            i = i + 1

    def add_row_to_light(self, filename):

        # internal functions
        def _bandpass_button_pressed(btn):
            selection = methods.create_bandpass_menu(btn).exec_(QtGui.QCursor.pos())
            if selection is not None:
                btn.setText(selection.objectName())

        def _bandpass_button_pressed_edit_curve(btn, rw, mw):
            _bandpass_button_pressed(btn)
            mw.loadobservations_widget.light_curves[rw].band_id = int(constants.BANDPASS_ID_DICT[str(btn.text())])
            mw.loadobservations_widget.update_curve_list()

        def _spinbox_edited(spnbx, rw, cl, mw):
            crv = mw.loadobservations_widget.light_curves[rw]
            if cl == 2:
                crv.l1 = spnbx.value()
            elif cl == 3:
                crv.l2 = spnbx.value()
            elif cl == 4:
                crv.l3 = spnbx.value()
            elif cl == 5:
                crv.x1 = spnbx.value()
            elif cl == 6:
                crv.x2 = spnbx.value()
            elif cl == 7:
                crv.y1 = spnbx.value()
            elif cl == 8:
                crv.y2 = spnbx.value()
            elif cl == 9:
                crv.opsf = spnbx.value()
            elif cl == 10:
                crv.aextinc = spnbx.value()
            elif cl == 11:
                crv.calib = spnbx.value()

        item = QtWidgets.QTreeWidgetItem(self.light_treewidget)
        item.setText(0, filename)

        bandpass_button = QtWidgets.QPushButton(self.light_treewidget)
        bandpass_button.setText("Johnson V")

        if filename != "[Synthetic]":
            bandpass_button.clicked.connect(partial(_bandpass_button_pressed_edit_curve,
                                                    bandpass_button, self.nlc, self.main_window))
        else:
            bandpass_button.clicked.connect(partial(_bandpass_button_pressed,
                                                    bandpass_button))

        self.light_treewidget.setItemWidget(item, 1, bandpass_button)

        i = 2
        while i < 14:
            spinbox = QtWidgets.QDoubleSpinBox(self.light_treewidget)
            spinbox.setMaximumWidth(90)
            spinbox.setButtonSymbols(2)
            spinbox.setDecimals(7)
            spinbox.setMaximum(999)
            spinbox.setMinimum(-999)
            if i == 12:
                spinbox.setValue(1.0)
            if i == 13:
                spinbox.setValue(8.0)
            if filename != "[Synthetic]":
                spinbox.editingFinished.connect(partial(_spinbox_edited, spinbox, self.nlc, i, self.main_window))
            self.light_treewidget.setItemWidget(item, i, spinbox)
            i = i + 1

        if filename != "[Synthetic]":
            self.nlc = self.nlc + 1

        return item

    def insert_synthetic_curves(self):
        self.add_row_to_light("[Synthetic]")
        self.add_div_to_light()
        self.vel_pri_label.setText("[Synthetic]")
        self.vel_sec_label.setText("[Synthetic]")

    def populate_from_loaded_curves(self):

        if self.obs_widget.velocity_curves[0] is not None:
            self.vel_pri_label.setText(os.path.basename(self.obs_widget.velocity_curves[0].filepath))

        if self.obs_widget.velocity_curves[1] is not None:
            self.vel_sec_label.setText(os.path.basename(self.obs_widget.velocity_curves[1].filepath))

        for lc in self.obs_widget.light_curves:
            item = self.add_row_to_light(os.path.basename(lc.filepath))
            self.light_treewidget.itemWidget(item, 1).setText(constants.ID_BANDPASS_DICT[str(lc.band_id)])
            self.light_treewidget.itemWidget(item, 2).setValue(lc.l1)
            self.light_treewidget.itemWidget(item, 3).setValue(lc.l2)
            self.light_treewidget.itemWidget(item, 4).setValue(lc.l3)
            self.light_treewidget.itemWidget(item, 5).setValue(lc.x1)
            self.light_treewidget.itemWidget(item, 6).setValue(lc.x2)
            self.light_treewidget.itemWidget(item, 7).setValue(lc.y1)
            self.light_treewidget.itemWidget(item, 8).setValue(lc.y2)
            self.light_treewidget.itemWidget(item, 9).setValue(lc.opsf)
            self.light_treewidget.itemWidget(item, 10).setValue(lc.aextinc)
            self.light_treewidget.itemWidget(item, 11).setValue(lc.calib)
            self.light_treewidget.itemWidget(item, 12).setValue(1.0)
            self.light_treewidget.itemWidget(item, 13).setValue(8.0)

    def reset_light_treewidget(self):
        self.nlc = 0
        self.light_treewidget.clear()
        #self.light_treewidget_2.clear()
        self.insert_synthetic_curves()

    def reset_widget(self):
        self.reset_light_treewidget()
        #self.light_chart.clear_all()

    def reset_and_repopulate(self):
        self.reset_widget()
        self.populate_from_loaded_curves()

    def selected_light_item(self):
        selecteditem = self.light_treewidget.selectedItems()
        if len(selecteditem) > 0:
            return selecteditem[0]
        else:
            return None

    def alias(self, obs, model_start, model_end):
        obs_x = obs[0]
        obs_y = obs[1]
        t0 = self.main_window.jd0_ipt.value()
        p = self.main_window.p0_ipt.value()
        n_x = None
        n_y = None

        if self.main_window.phase_radiobtn.isChecked():
            if self.main_window.jdphs_combobox.currentText() == "Phase":
                n_x, n_y = methods.alias_phased_obs_with_phase(obs_x, obs_y, model_start, model_end)

            elif self.main_window.jdphs_combobox.currentText() == "Time":
                n_x, n_y = methods.alias_jd_obs_with_phase(obs_x, obs_y, model_start, model_end, t0, p)

        elif self.main_window.jd_radiobtn.isChecked():
            if self.main_window.jdphs_combobox.currentText() == "Phase":
                n_x, n_y = methods.alias_phased_obs_with_jd(obs_x, obs_y, model_start, model_end, t0, p)

            elif self.main_window.jdphs_combobox.currentText() == "Time":
                n_x = obs_x
                n_y = obs_y

        return n_x, n_y

    def plot_selected_light_curve(self):
        selected_item = self.selected_light_item()
        if selected_item is not None:
            if self.light_treewidget.invisibleRootItem().indexOfChild(selected_item) == 1:
                return 0
            self.light_chart.clear_all()
            lc_params = self.main_window.get_lc_params()
            model_start, model_end = self.main_window.get_lc_boundaries()

            x_index = None
            if self.main_window.jd_radiobtn.isChecked():
                x_index = 0

            elif self.main_window.phase_radiobtn.isChecked():
                x_index = 1

            y_index = None
            if self.main_window.maglite_combobox.currentText() == "Magnitude":
                y_index = 8

            elif self.main_window.maglite_combobox.currentText() == "Flux":
                y_index = 4

            lc_params.set_synthetic_curve(
                int(constants.BANDPASS_ID_DICT[self.light_treewidget.itemWidget(selected_item, 1).text()]),
                self.light_treewidget.itemWidget(selected_item, 2).value(),  # l1
                self.light_treewidget.itemWidget(selected_item, 3).value(),  # l2
                self.light_treewidget.itemWidget(selected_item, 5).value(),  # x1
                self.light_treewidget.itemWidget(selected_item, 6).value(),  # x2
                self.light_treewidget.itemWidget(selected_item, 7).value(),  # y1
                self.light_treewidget.itemWidget(selected_item, 8).value(),  # y2
                self.light_treewidget.itemWidget(selected_item, 4).value(),  # l3
                self.light_treewidget.itemWidget(selected_item, 9).value(),  # opsf
                self.light_treewidget.itemWidget(selected_item, 13).value(),  # zero
                self.light_treewidget.itemWidget(selected_item, 12).value(),  # factor
                0.55,  # wl, dummy
                self.light_treewidget.itemWidget(selected_item, 10).value(),  # aextinc
                self.light_treewidget.itemWidget(selected_item, 11).value()  # calib
            )

            if self.light_plotobs_chk.isChecked() and selected_item.text(0) != "[Synthetic]":
                index = self.light_treewidget.invisibleRootItem().indexOfChild(selected_item)
                data = self.main_window.loadobservations_widget.light_curves[index - 2].get_data()
                lc_obs = data[0], data[1]

                if self.light_alias_chk.isChecked():
                    try:
                        lc_obs = self.alias(lc_obs, model_start, model_end)
                    except ValueError as e:
                        msg = messenger.Messenger("error", "A ValueError has occured:")
                        msg.set_info(e.args[0])
                        msg.show()

                        return 0

                self.light_chart.plot(lc_obs[0], lc_obs[1],
                                      linestyle="", marker="o",
                                      markersize=constants.MARKER_SIZE, color=constants.COLOR_BLUE)

            if self.light_plotmodel_chk.isChecked():
                lc_io = wd_io.LCIO(lc_params,
                                   wd_path=self.main_window.lc_path,
                                   lc_binary_name=self.main_window.lc_binary)

                results = lc_io.fill_for_synthetic_light_curve().save().run().read_synthetic_light_curve()

                absolute_params, teffs, sma, lds, lums = lc_io.read_abs_params()
                teffs = float(teffs[0][0])*10000, float(teffs[1][0])*10000
                sma = float(sma[1][0])
                L1, L2, logL1, logL2 = methods.compute_luminosity(teffs[0],teffs[1],absolute_params[2][0],absolute_params[2][1])
                #self.light_treewidget_2.clear()
                item = self.light_treewidget_2

                sma = str(sma)
                a = str(absolute_params[1][0])
                b = str(absolute_params[1][1])
                c = str(absolute_params[2][0])
                d = str(absolute_params[2][1])
                e = str(int(teffs[0]))
                f = str(int(teffs[1]))

                aa = str(absolute_params[4][0])
                bb = str(absolute_params[4][1])
                cc = str(absolute_params[3][0])
                dd = str(absolute_params[3][1])
                ee = str(logL1)
                ff = str(logL2)


                for index, val in enumerate((sma, a, b, c, d, e, f, aa, bb, cc, dd, ee, ff)):
                    if val == "nan":
                        item.topLevelItem(index).setBackground(index, QtGui.QBrush(QtGui.QColor("red")))
                        val = "NaN"
                    item.topLevelItem(index).setText(1, val)

                self.light_chart.plot(results[x_index], results[y_index], clear=False, color=constants.COLOR_RED)

            if self.light_plotobs_chk.isChecked() and self.light_plotmodel_chk.isChecked() and \
                    selected_item.text(0) != "[Synthetic]":
                residuals = methods.compute_residuals(lc_obs[0], lc_obs[1], results[x_index], results[y_index])
                self.light_chart.plot(lc_obs[0], residuals, index=1, linestyle="", marker="o",
                                      markersize=constants.MARKER_SIZE, color=constants.COLOR_BLUE)

                self.light_chart.axes[1].axhline([0], color=constants.COLOR_RED)
                self.light_chart.redraw()

            x_label = ""
            y_label = self.main_window.maglite_combobox.currentText()

            if y_label == "Magnitude":
                self.light_chart.axes[0].invert_yaxis()
                self.light_chart.axes[1].invert_yaxis()

            if x_index == 0:
                x_label = "HJD"

            elif x_index == 1:
                x_label = "Phase"

            self.light_chart.set_labels("", y_label, index=0)
            self.light_chart.set_labels(x_label, "Residuals", index=1)

    def plot_selected_velocity_curve(self):
        self.velocity_chart.clear_all()
        lc_params = self.main_window.get_lc_params()
        model_start, model_end = self.main_window.get_lc_boundaries()

        x_index = None
        if self.main_window.jd_radiobtn.isChecked():
            x_index = 0

        elif self.main_window.phase_radiobtn.isChecked():
            x_index = 1

        lc_params.set_dummy_synthetic_curve()

        plot_vc1 = False
        plot_vc2 = False

        if self.vel_pri_radiobtn.isChecked():
            plot_vc1 = True

        elif self.vel_sec_radiobtn.isChecked():
            plot_vc2 = True

        if self.vel_both_radiobtn.isChecked():
            plot_vc1 = True
            plot_vc2 = True

        vunit = float(self.main_window.vunit_ipt.value())
        if self.vel_plotobs_chk.isChecked():
            if plot_vc1 and self.vel_pri_label.text() != "[Synthetic]":
                data = self.main_window.loadobservations_widget.velocity_curves[0].get_data()
                s1_obs = data[0], [i * vunit for i in data[1]]
                if self.vel_alias_chk.isChecked():
                    try:
                        s1_obs = self.alias(s1_obs, model_start, model_end)
                    except ValueError as e:
                        msg = messenger.Messenger("error", "A ValueError has occured:")
                        msg.set_info(e.args[0])
                        msg.show()

                        return 0

                self.velocity_chart.plot(s1_obs[0], s1_obs[1],
                                         linestyle="", marker="o",
                                         markersize=constants.MARKER_SIZE, color=constants.COLOR_BLUE)

            if plot_vc2 and self.vel_sec_label.text() != "[Synthetic]":
                data = self.main_window.loadobservations_widget.velocity_curves[1].get_data()
                s2_obs = data[0], [i * vunit for i in data[1]]
                if self.vel_alias_chk.isChecked():
                    try:
                        s2_obs = self.alias(s2_obs, model_start, model_end)
                    except ValueError as e:
                        msg = messenger.Messenger("error", "A ValueError has occured:")
                        msg.set_info(e.args[0])
                        msg.show()

                        return 0

                self.velocity_chart.plot(s2_obs[0], s2_obs[1], clear=False,
                                         linestyle="", marker="o",
                                         markersize=constants.MARKER_SIZE, color=constants.COLOR_GREEN)

        if self.vel_plotmodel_chk.isChecked():
            lc_io = wd_io.LCIO(lc_params,
                               wd_path=self.main_window.lc_path,
                               lc_binary_name=self.main_window.lc_binary)

            results = lc_io.fill_for_synthetic_velocity_curve().save().run().read_synthetic_velocity_curve()
            s1_index = 6
            s2_index = 7
            results[s1_index] = [i * vunit for i in results[s1_index]]
            results[s2_index] = [i * vunit for i in results[s2_index]]

            if plot_vc1:
                self.velocity_chart.plot(results[x_index], results[s1_index], clear=False, color=constants.COLOR_RED)

                s1_model = results[x_index], results[s1_index]

            if plot_vc2:
                self.velocity_chart.plot(results[x_index], results[s2_index], clear=False, color=constants.COLOR_ORANGE)

                s2_model = results[x_index], results[s2_index]

        if self.vel_plotobs_chk.isChecked() and self.vel_plotmodel_chk.isChecked():
            if plot_vc1 and self.vel_pri_label.text() != "[Synthetic]":
                residuals = methods.compute_residuals(s1_obs[0], s1_obs[1], s1_model[0], s1_model[1])
                self.velocity_chart.plot(s1_obs[0], residuals, index=1,
                                         linestyle="", marker="o",
                                         markersize=constants.MARKER_SIZE, color=constants.COLOR_BLUE)

            if plot_vc2 and self.vel_sec_label.text() != "[Synthetic]":
                residuals = methods.compute_residuals(s2_obs[0], s2_obs[1], s2_model[0], s2_model[1])
                self.velocity_chart.plot(s2_obs[0], residuals, index=1, clear=False,
                                         linestyle="", marker="o",
                                         markersize=constants.MARKER_SIZE, color=constants.COLOR_GREEN)

            self.velocity_chart.axes[1].axhline([0], color=constants.COLOR_RED)
            self.velocity_chart.redraw()

        x_label = ""
        y_label = "Km/s"

        if x_index == 0:
            x_label = "HJD"

        elif x_index == 1:
            x_label = "Phase"

        self.velocity_chart.set_labels("", y_label, index=0)
        self.velocity_chart.set_labels(x_label, "Residuals", index=1)

        self.velocity_chart.axes[0].axhline([self.main_window.vgamma_ipt.value()*vunit], color="black", linestyle="--")
        self.velocity_chart.redraw()
