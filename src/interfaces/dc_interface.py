from PyQt5 import QtWidgets, QtGui, QtCore
from gui import dc_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src import constants
from src.helpers.wd_utils import wd_io
from src.helpers import methods
from src.helpers import messenger
import os
import math
import numpy
import webbrowser  # sending open file commands with this might not be portable
import timeit

delimiter = constants.EXPORT_DELIMITER


class Widget(QtWidgets.QWidget, dc_widget.Ui_DCWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.ptm_font = QtGui.QFont(parent.monoFont)
        self.ptm_font.setPointSize(constants.DC_RESULTS_FONTSIZE)
        self.component_treewidget.setFont(self.ptm_font)
        self.curvestat_treewidget.setFont(self.ptm_font)
        self.residual_treewidget.setFont(self.ptm_font)
        self.result_treewidget.setFont(self.ptm_font)

        self.component_treewidget.header().setSectionResizeMode(3)
        self.curvestat_treewidget.header().setSectionResizeMode(3)
        self.residual_treewidget.header().setSectionResizeMode(3)
        self.result_treewidget.header().setSectionResizeMode(3)

        self.component_treewidget.header().setSectionsMovable(False)
        self.curvestat_treewidget.header().setSectionsMovable(False)
        self.residual_treewidget.header().setSectionsMovable(False)
        self.result_treewidget.header().setSectionsMovable(False)

        # setup chart area
        self.chart = MatplotlibWidget(self.plotwidget, 2, 1, h_ratios=[1.5, 1])
        self.chart.create_axis(0, 0, name="Solution Curve")
        self.chart.create_axis(1, 0, sharex=self.chart[0])
        self.chart.figure.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95, hspace=0, wspace=0)
        self.chart[0].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True,
                                  labelright=False, labelleft=True, labelsize=11)

        self.current_iteration = 1
        self.iterator = None

        self.results = None
        self.dimensions = None
        self.observations = None
        self.solution_stats = None
        self.curve_stats = None

        self.set_default_dels()

        self.connect_signals()

    def connect_signals(self):
        self.rundc2015_btn.clicked.connect(self.start_differential_correction)
        self.plot_btn.clicked.connect(self.plot)
        self.viewlastdcin_btn.clicked.connect(self.open_dcin)
        self.viewlaastdcout_btn.clicked.connect(self.open_dcout)
        self.updateinputs_btn.clicked.connect(self.update_inputs)
        self.clearbaseset_btn.clicked.connect(self.clear_keeps)
        self.setdeldefaults_btn.clicked.connect(self.set_default_dels)
        self.exportresults_btn.clicked.connect(self.export_data)

    def open_dcin(self):
        webbrowser.open(os.path.join(self.main_window.dc_path, "dcin.active"))

    def open_dcout(self):
        webbrowser.open(os.path.join(self.main_window.dc_path, "dcout.active"))

    def export_data(self):
        if self.results is not None:
            menu = QtWidgets.QMenu(self)
            plaintext = menu.addAction("Plaintext")
            plaintext.setObjectName("plaintext")
            latex = menu.addAction("Latex")
            latex.setObjectName("latex")
            selection = menu.exec_(QtGui.QCursor.pos())

            if selection is not None:
                filepath = methods.save_file(self, suffix="txt")

                if filepath is not None:
                    with open(filepath, "w") as f:
                        if selection.objectName() == "plaintext":
                            f.write("#Parameter" + delimiter + "Value" + delimiter + "Sigma\n")
                        if selection.objectName() == "latex":
                            f.write("\\begin{table}\n\\begin{center}\n\\begin{tabular}{c|c}\n")
                            f.write("Parameter & Value" + "\\" + "\\" + "\n")
                            f.write("\hline")

                        for result in self.results:

                            output = result[2]
                            stderr = result[5]
                            par_id = result[0]
                            c_id = int(result[1])

                            if par_id in (19.0, 20.0):
                                output = int(output * 10000.0)
                                stderr = int(stderr * 10000.0)

                            if selection.objectName() == "plaintext":
                                name = constants.KEEPS_ID_NAME_DICT[par_id]
                                if c_id != 0:
                                    band = constants.ID_BANDPASS_DICT[
                                        str(self.main_window.loadobservations_widget.light_curves[c_id - 1].band_id)
                                    ]
                                    name = name + " (" + band + ")"
                                f.write(name + delimiter + str(output) + delimiter + str(stderr) + "\n")

                            if selection.objectName() == "latex":
                                name = constants.ID_LATEX_DICT[str(int(par_id))]
                                if name != "nan":
                                    if par_id in (1.0, 2.0, 3.0, 5.0, 6.0, 7.0, 11.0):
                                        output = output * 180.0 / numpy.pi
                                        stderr = stderr * 180.0 / numpy.pi
                                    if c_id != 0:
                                        band = constants.ID_BANDPASS_DICT[
                                            str(self.main_window.loadobservations_widget.light_curves[c_id - 1].band_id)
                                        ]
                                        name = name.format(band=band.replace(" ", "~"))
                                    if par_id == 23.0 and self.main_window.mode_combobox.currentText() in ("Mode 1", "Mode 3"):
                                        name = "$\Omega_{1}$ = $\Omega_{2}$"
                                    f.write(name + " & " + str(output) + " $\pm$ " + str(stderr) + " \\" + "\\" + "\n")

                        f.write("\n")

                        s1_dimensions = list(numpy.array(self.dimensions[0]).transpose())
                        s2_dimensions = list(numpy.array(self.dimensions[1]).transpose())

                        for i, component in enumerate((s1_dimensions, s2_dimensions)):
                            for result in component:
                                if selection.objectName() == "plaintext":
                                    value = "r" + str(i + 1) + "_" + str(result[1])
                                    f.write(value + delimiter + str(result[2]) + delimiter + str(result[3]) + "\n")
                                elif selection.objectName() == "latex":
                                    value = "$r_{" + str(i + 1) + "~" + str(result[1]) + "}$"
                                    f.write(value + " & " + str(result[2]) + " $\pm$ " +
                                            str(result[3]) + " \\" + "\\" + "\n")

                        if selection.objectName() == "plaintext":
                            f.write("\n#Mean Residual for Input Values\n" + str(self.solution_stats[0][0]) + "\n")

                        elif selection.objectName() == "latex":
                            f.write("\n$Mean~residual~for~input~values$ & " + str(self.solution_stats[0][0])
                                    + "\\" + "\\" + "\n")
                            f.write("\\end{tabular}\n\\end{center}\n\\end{table}")

                    msg = messenger.Messenger("info", "File saved:")
                    info = filepath
                    if selection.objectName() == "latex":
                        info = info + "\n\n" + "Exported values are not truncated and exported as-is. " \
                                             "Please inspect your solutions carefully and " \
                                             "truncate your results before using the Latex file."
                    msg.set_info(info)
                    msg.show()

    def clear(self):
        self.result_treewidget.clear()
        self.residual_treewidget.clear()
        self.curvestat_treewidget.clear()
        self.component_treewidget.clear()

        self.results = None
        self.dimensions = None
        self.observations = None
        self.solution_stats = None
        self.curve_stats = None

        self.update_curve_list()

        self.chart.clear_all()

        self.clear_iterator()

    def clear_keeps(self):
        self.s1lat_chk.setChecked(False)
        self.s1long_chk.setChecked(False)
        self.s1rad_chk.setChecked(False)
        self.s1temp_chk.setChecked(False)

        self.s2lat_chk.setChecked(False)
        self.s2long_chk.setChecked(False)
        self.s2rad_chk.setChecked(False)
        self.s2temp_chk.setChecked(False)

        self.a_chk.setChecked(False)
        self.e_chk.setChecked(False)
        self.perr0_chk.setChecked(False)
        self.f1_chk.setChecked(False)
        self.f2_chk.setChecked(False)
        self.pshift_chk.setChecked(False)
        self.vgam_chk.setChecked(False)
        self.incl_chk.setChecked(False)
        self.g1_chk.setChecked(False)
        self.g2_chk.setChecked(False)
        self.t1_chk.setChecked(False)
        self.t2_chk.setChecked(False)
        self.alb1_chk.setChecked(False)
        self.alb2_chk.setChecked(False)
        self.pot1_chk.setChecked(False)
        self.pot2_chk.setChecked(False)
        self.q_chk.setChecked(False)
        self.jd0_chk.setChecked(False)
        self.p0_chk.setChecked(False)
        self.dpdt_chk.setChecked(False)
        self.dperdt_chk.setChecked(False)
        self.a3b_chk.setChecked(False)
        self.p3b_chk.setChecked(False)
        self.xinc3b_chk.setChecked(False)
        self.e3b_chk.setChecked(False)
        self.perr3b_chk.setChecked(False)
        self.jd0_chk.setChecked(False)
        self.logd_chk.setChecked(False)
        self.desextinc_chk.setChecked(False)

        self.s1tstart_chk.setChecked(False)
        self.s1tmax1_chk.setChecked(False)
        self.s1tmax2_chk.setChecked(False)
        self.s1tend_chk.setChecked(False)

        self.s2tstart_chk.setChecked(False)
        self.s2tmax1_chk.setChecked(False)
        self.s2tmax2_chk.setChecked(False)
        self.s2tend_chk.setChecked(False)

        self.l1_chk.setChecked(False)
        self.l2_chk.setChecked(False)
        self.x1_chk.setChecked(False)
        self.x2_chk.setChecked(False)
        self.el3_chk.setChecked(False)

        # iteration params
        self.marqmul_ipt.setValue(0.00001)
        self.vlr_ipt.setValue(1.0)

    def set_default_dels(self):
        self.del_s1lat_ipt.setValue(0.02)
        self.del_s1lng_ipt.setValue(0.02)
        self.del_s1agrad_ipt.setValue(0.001)
        self.del_s1tmpf_ipt.setValue(0.02)

        self.del_s2lat_ipt.setValue(0.02)
        self.del_s2lng_ipt.setValue(0.02)
        self.del_s2agrad_ipt.setValue(0.001)
        self.del_s2tmpf_ipt.setValue(0.02)

        self.del_a_ipt.setValue(0.0500)
        self.del_e_ipt.setValue(0.0010)
        self.del_perr0_ipt.setValue(0.0100)
        self.del_f1_ipt.setValue(0.0100)
        self.del_f2_ipt.setValue(0.0100)
        self.del_pshift_ipt.setValue(0.0020)
        self.del_i_ipt.setValue(0.2000)
        self.del_g1_ipt.setValue(0.0100)
        self.del_g2_ipt.setValue(0.0100)
        self.del_t1_ipt.setValue(0.0200)
        self.del_t2_ipt.setValue(0.0200)

        self.del_alb1_ipt.setValue(0.0500)
        self.del_alb2_ipt.setValue(0.0500)
        self.del_pot1_ipt.setValue(0.0200)
        self.del_pot2_ipt.setValue(0.0200)
        self.del_q_ipt.setValue(0.0030)
        self.del_l1_ipt.setValue(0.0100)
        self.del_l2_ipt.setValue(0.0100)
        self.del_x1_ipt.setValue(0.0100)
        self.del_x2_ipt.setValue(0.0100)

    def get_dc_params(self):
        def _eval_checkbox(checkbox):
            if checkbox.isChecked():
                return 0
            else:
                return 1

        dc_params = self.main_window.get_dc_params()

        dc_params.keeps["spot_a_lat"] = _eval_checkbox(self.s1lat_chk)
        dc_params.keeps["spot_a_long"] = _eval_checkbox(self.s1long_chk)
        dc_params.keeps["spot_a_rad"] = _eval_checkbox(self.s1rad_chk)
        dc_params.keeps["spot_a_tempf"] = _eval_checkbox(self.s1temp_chk)

        dc_params.keeps["spot_b_lat"] = _eval_checkbox(self.s2lat_chk)
        dc_params.keeps["spot_b_long"] = _eval_checkbox(self.s2long_chk)
        dc_params.keeps["spot_b_rad"] = _eval_checkbox(self.s2rad_chk)
        dc_params.keeps["spot_b_tempf"] = _eval_checkbox(self.s2temp_chk)

        dc_params.keeps["a"] = _eval_checkbox(self.a_chk)
        dc_params.keeps["e"] = _eval_checkbox(self.e_chk)
        dc_params.keeps["perr"] = _eval_checkbox(self.perr0_chk)
        dc_params.keeps["f1"] = _eval_checkbox(self.f1_chk)
        dc_params.keeps["f2"] = _eval_checkbox(self.f2_chk)
        dc_params.keeps["pshift"] = _eval_checkbox(self.pshift_chk)
        dc_params.keeps["vga"] = _eval_checkbox(self.vgam_chk)
        dc_params.keeps["xincl"] = _eval_checkbox(self.incl_chk)
        dc_params.keeps["g1"] = _eval_checkbox(self.g1_chk)
        dc_params.keeps["g2"] = _eval_checkbox(self.g2_chk)
        dc_params.keeps["tavh"] = _eval_checkbox(self.t1_chk)
        dc_params.keeps["tavc"] = _eval_checkbox(self.t2_chk)
        dc_params.keeps["alb1"] = _eval_checkbox(self.alb1_chk)
        dc_params.keeps["alb2"] = _eval_checkbox(self.alb2_chk)
        dc_params.keeps["phsv"] = _eval_checkbox(self.pot1_chk)
        dc_params.keeps["pcsv"] = _eval_checkbox(self.pot2_chk)
        dc_params.keeps["rm"] = _eval_checkbox(self.q_chk)
        dc_params.keeps["hjd0"] = _eval_checkbox(self.jd0_chk)
        dc_params.keeps["pzero"] = _eval_checkbox(self.p0_chk)
        dc_params.keeps["dpdt"] = _eval_checkbox(self.dpdt_chk)
        dc_params.keeps["dperdt"] = _eval_checkbox(self.dperdt_chk)
        dc_params.keeps["a3b"] = _eval_checkbox(self.a3b_chk)
        dc_params.keeps["p3b"] = _eval_checkbox(self.p3b_chk)
        dc_params.keeps["xincl3b"] = _eval_checkbox(self.xinc3b_chk)
        dc_params.keeps["e3b"] = _eval_checkbox(self.e3b_chk)
        dc_params.keeps["perr3b"] = _eval_checkbox(self.perr3b_chk)
        dc_params.keeps["t03b"] = _eval_checkbox(self.tc3b_chk)
        dc_params.keeps["dpclog"] = _eval_checkbox(self.logd_chk)
        dc_params.keeps["desextinc"] = _eval_checkbox(self.desextinc_chk)

        dc_params.keeps["spot_a_tstart"] = _eval_checkbox(self.s1tstart_chk)
        dc_params.keeps["spot_a_tmax1"] = _eval_checkbox(self.s1tmax1_chk)
        dc_params.keeps["spot_a_tmax2"] = _eval_checkbox(self.s1tmax2_chk)
        dc_params.keeps["spot_a_tend"] = _eval_checkbox(self.s1tend_chk)

        dc_params.keeps["spot_b_tstart"] = _eval_checkbox(self.s2tstart_chk)
        dc_params.keeps["spot_b_tmax1"] = _eval_checkbox(self.s2tmax1_chk)
        dc_params.keeps["spot_b_tmax2"] = _eval_checkbox(self.s2tmax2_chk)
        dc_params.keeps["spot_b_tend"] = _eval_checkbox(self.s2tend_chk)

        dc_params.keeps["hla"] = _eval_checkbox(self.l1_chk)
        dc_params.keeps["cla"] = _eval_checkbox(self.l2_chk)
        dc_params.keeps["x1a"] = _eval_checkbox(self.x1_chk)
        dc_params.keeps["x2a"] = _eval_checkbox(self.x2_chk)
        dc_params.keeps["el3a"] = _eval_checkbox(self.el3_chk)

        # iteration params
        dc_params.keeps["niter"] = self.niter_spinbox.value()
        dc_params.keeps["xlamda"] = self.marqmul_ipt.value()
        dc_params.keeps["vlr"] = self.vlr_ipt.value()

        # dels
        dc_params.dels["spot_a_lat"] = self.del_s1lat_ipt.value()
        dc_params.dels["spot_a_long"] = self.del_s1lng_ipt.value()
        dc_params.dels["spot_a_rad"] = self.del_s1agrad_ipt.value()
        dc_params.dels["spot_a_tempf"] = self.del_s1tmpf_ipt.value()

        dc_params.dels["spot_b_lat"] = self.del_s2lat_ipt.value()
        dc_params.dels["spot_b_long"] = self.del_s2lng_ipt.value()
        dc_params.dels["spot_b_rad"] = self.del_s2agrad_ipt.value()
        dc_params.dels["spot_b_tempf"] = self.del_s2tmpf_ipt.value()

        dc_params.dels["a"] = self.del_a_ipt.value()
        dc_params.dels["e"] = self.del_e_ipt.value()
        dc_params.dels["perr"] = self.del_perr0_ipt.value()
        dc_params.dels["f1"] = self.del_f1_ipt.value()
        dc_params.dels["f2"] = self.del_f2_ipt.value()
        dc_params.dels["pshift"] = self.del_pshift_ipt.value()
        dc_params.dels["xincl"] = self.del_i_ipt.value()
        dc_params.dels["g1"] = self.del_g1_ipt.value()
        dc_params.dels["g2"] = self.del_g2_ipt.value()
        dc_params.dels["tavh"] = self.del_t1_ipt.value()
        dc_params.dels["tavc"] = self.del_t2_ipt.value()

        dc_params.dels["alb1"] = self.del_alb1_ipt.value()
        dc_params.dels["alb2"] = self.del_alb2_ipt.value()
        dc_params.dels["phsv"] = self.del_pot1_ipt.value()
        dc_params.dels["pcsv"] = self.del_pot2_ipt.value()
        dc_params.dels["rm"] = self.del_q_ipt.value()
        dc_params.dels["hla"] = self.del_l1_ipt.value()
        dc_params.dels["cla"] = self.del_l2_ipt.value()
        dc_params.dels["x1a"] = self.del_x1_ipt.value()
        dc_params.dels["x2a"] = self.del_x2_ipt.value()

        return dc_params

    def update_inputs(self):
        if self.results is None:
            msg = messenger.Messenger("error", "Can't update inputs:")
            msg.set_info("No result found.")
            msg.show()
        elif type(self.results) == float:
            if math.isnan(self.results):
                msg = messenger.Messenger("error", "Can't update inputs:")
                msg.set_info("Last iteration returned NaN's.")
                msg.show()
        else:
            def _update_interface(par_id, value, main_window):
                paramdict = {
                    9: main_window.a_ipt,
                    10: main_window.e_ipt,
                    11: main_window.perr0_ipt,
                    12: main_window.f1_ipt,
                    13: main_window.f2_ipt,
                    14: main_window.pshift_ipt,
                    15: main_window.vgamma_ipt,
                    16: main_window.incl_ipt,
                    17: main_window.gr1_ipt,
                    18: main_window.gr2_ipt,
                    19: main_window.t1_ipt,
                    20: main_window.t2_ipt,
                    21: main_window.alb1_ipt,
                    22: main_window.alb2_ipt,
                    23: main_window.pot1_ipt,
                    24: main_window.pot2_ipt,
                    25: main_window.q_ipt,
                    26: main_window.jd0_ipt,
                    27: main_window.p0_ipt,
                    28: main_window.dpdt_ipt,
                    29: main_window.dperdt_ipt,
                    30: main_window.a3b_ipt,
                    31: main_window.p3b_ipt,
                    32: main_window.incl3b_ipt,
                    33: main_window.e3b_ipt,
                    34: main_window.perr03b_ipt,
                    35: main_window.conj3b_ipt,
                    41: main_window.dpclog_ipt,
                    42: main_window.dc_desextinc_ipt
                }

                if par_id in (19.0, 20.0):
                    value = value * 10000.0

                paramdict[par_id].setValue(value)

            def _update_curve(par_id, c_id, value, main_window):
                curve = main_window.loadobservations_widget.light_curves[c_id - 1]
                if par_id == 56:
                    curve.l1 = value

                elif par_id == 57:
                    curve.l2 = value

                elif par_id == 58:
                    curve.x1 = value

                elif par_id == 59:
                    curve.x2 = value

                elif par_id == 60:
                    curve.l3 = value

            def _update_spot(par_id, value, main_window):

                spot_idx_paramdict = {
                    1: 2,
                    2: 3,
                    3: 4,
                    4: 5,
                    5: 2,
                    6: 3,
                    7: 4,
                    8: 5,
                    43: 6,
                    44: 7,
                    45: 8,
                    46: 9,
                    47: 6,
                    48: 7,
                    49: 8,
                    50: 9,
                }

                spot_a = (1, 2, 3, 4, 43, 44, 45, 46)
                spot_b = (5, 6, 7, 8, 47, 48, 49, 50)

                radio_index = None

                if par_id in spot_a:
                    radio_index = 0

                elif par_id in spot_b:
                    radio_index = 1

                s1_parent = main_window.configurespot_widget.star1_treewidget
                s2_parent = main_window.configurespot_widget.star2_treewidget

                s1_spot_items = s1_parent.get_all_items()
                s2_spot_items = s2_parent.get_all_items()

                for item in s1_spot_items:
                    radio_btn = s1_parent.tree_widget.itemWidget(item, radio_index)
                    if radio_btn.isChecked():
                        ipt = s1_parent.tree_widget.itemWidget(item, spot_idx_paramdict[par_id])
                        ipt.setValue(value)

                for item in s2_spot_items:
                    radio_btn = s2_parent.tree_widget.itemWidget(item, radio_index)
                    if radio_btn.isChecked():
                        ipt = s2_parent.tree_widget.itemWidget(item, spot_idx_paramdict[par_id])
                        ipt.setValue(value)

            for result in self.results:
                if result[1] == 0.0:
                    if result[0] in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                                     43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0):
                        _update_spot(int(result[0]), result[4], self.main_window)
                    else:
                        _update_interface(int(result[0]), result[4], self.main_window)

                elif result[1] != 0.0:
                    _update_curve(int(result[0]), int(result[1]), result[4], self.main_window)

            self.main_window.lc_synthetic_curve_widget.reset_and_repopulate()

    def plot(self):
        if self.observations is not None:
            if self.data_combobox.currentText() == "Velocity Curve #1 + #2":
                self.plot_both_vc()

            elif self.data_combobox.currentText() != "":
                self.plot_curve()

    def plot_both_vc(self):
        self.chart.clear_all()

        dcout_obs_x_index = 0
        dcout_obs_y_index = 1
        dcout_mdl_y_index = 2

        if self.main_window.jdphs_combobox.currentText() == "Time":
            dcout_obs_y_index = 2
            dcout_mdl_y_index = 3

            if self.forcephase_chk.isChecked():
                dcout_obs_x_index = 1

        vc1_data = self.observations[0]
        vc2_data = self.observations[1]

        vc1_obs_x = vc1_data[dcout_obs_x_index]
        vc2_obs_x = vc2_data[dcout_obs_x_index]

        vc1_obs_y = vc1_data[dcout_obs_y_index]
        vc2_obs_y = vc2_data[dcout_obs_y_index]

        vc1_resd = vc1_data[-1]
        vc2_resd = vc2_data[-1]

        vc1_mdl_x = None
        vc2_mdl_x = None
        vc1_mdl_y = None
        vc2_mdl_y = None

        if self.uselc_chk.isChecked():
            x_index = None
            y_index = None

            data = numpy.append(vc1_obs_x, vc2_obs_x)

            mdl_start = min(data)
            mdl_end = max(data)

            lc_params = self.main_window.get_lc_params()

            if self.main_window.jdphs_combobox.currentText() == "Time":
                x_index = 0

                lc_params["jdphs"] = 1
                lc_params["hjdst"] = float(mdl_start)
                lc_params["hjdsp"] = float(mdl_end)
                lc_params["hjdin"] = float(self.main_window.p0_ipt.value()) / 500.0

            if self.forcephase_chk.isChecked() or self.main_window.jdphs_combobox.currentText() == "Phase":
                x_index = 1

                lc_params["jdphs"] = 2
                lc_params["phstrt"] = 0.0
                lc_params["phstop"] = 1.0
                lc_params["phin"] = 0.001

            lc_params.set_dummy_synthetic_curve()

            lc_io = wd_io.LCIO(lc_params,
                               wd_path=self.main_window.lc_path,
                               lc_binary_name=self.main_window.lc_binary)

            results = lc_io.fill_for_synthetic_velocity_curve().save().run().read_synthetic_velocity_curve()

            vc1_mdl_x = results[x_index]
            vc2_mdl_x = results[x_index]
            vc1_mdl_y = results[6]
            vc2_mdl_y = results[7]

        else:
            vc1_mdl_x = vc1_obs_x
            vc2_mdl_x = vc2_obs_x
            vc1_mdl_y = vc1_data[dcout_mdl_y_index]
            vc2_mdl_y = vc2_data[dcout_mdl_y_index]

        self.chart.plot(vc1_obs_x, vc1_obs_y, clear=False,
                        markersize=constants.MARKER_SIZE, marker="o", linestyle="", color=constants.COLOR_BLUE)

        self.chart.plot(vc2_obs_x, vc2_obs_y, clear=False,
                        markersize=constants.MARKER_SIZE, marker="o", linestyle="", color=constants.COLOR_GREEN)

        if self.uselc_chk.isChecked():
            self.chart.plot(vc1_mdl_x, vc1_mdl_y, clear=False, color=constants.COLOR_RED)
            self.chart.plot(vc2_mdl_x, vc2_mdl_y, clear=False, color=constants.COLOR_ORANGE)

        else:
            self.chart.plot(vc1_mdl_x, vc1_mdl_y, clear=False,
                            markersize=constants.MARKER_SIZE - 1, marker="o",
                            linestyle="", color=constants.COLOR_RED)

            self.chart.plot(vc2_mdl_x, vc2_mdl_y, clear=False,
                            markersize=constants.MARKER_SIZE - 1, marker="o",
                            linestyle="", color=constants.COLOR_RED)

        self.chart.plot(vc1_obs_x, vc1_resd, index=1, clear=False,
                        markersize=constants.MARKER_SIZE, marker="o", linestyle="", color=constants.COLOR_BLUE)

        self.chart.plot(vc2_obs_x, vc2_resd, index=1, clear=False,
                        markersize=constants.MARKER_SIZE, marker="o", linestyle="", color=constants.COLOR_GREEN)

        self.chart.axes[1].axhline([0], color=constants.COLOR_RED)
        self.chart.axes[0].axhline([self.main_window.vgamma_ipt.value()], color="black", linestyle="--")

        self.chart.set_labels("", "Km/s", index=0)

        if self.forcephase_chk.isChecked():
            self.chart.set_labels("Phase", "Residuals", index=1)

        else:
            if self.main_window.jdphs_combobox.currentText() == "Time":
                self.chart.set_labels("HJD", "Residuals", index=1)
            else:
                self.chart.set_labels("Phase", "Residuals", index=1)

        self.chart.redraw()

    def plot_curve(self):
        self.chart.clear_all()

        curves = self.main_window.loadobservations_widget.get_all_curves()
        curve = curves[self.data_combobox.currentIndex()]
        mdl_data = self.observations[self.data_combobox.currentIndex()]

        dcout_obs_x_index = 0
        dcout_obs_y_index = 1
        dcout_mdl_y_index = 2
        dcout_residual_index = -1

        if self.main_window.jdphs_combobox.currentText() == "Time":
            dcout_obs_y_index = 2
            dcout_mdl_y_index = 3

            if self.forcephase_chk.isChecked():
                dcout_obs_x_index = 1

        obs_x = mdl_data[dcout_obs_x_index]
        obs_y = mdl_data[dcout_obs_y_index]
        residuals = mdl_data[dcout_residual_index]

        mdl_x = None
        mdl_y = None

        obs_color = constants.COLOR_BLUE
        mdl_color = constants.COLOR_RED

        if (self.data_combobox.currentIndex() == 0 and
            self.main_window.loadobservations_widget.velocity_curves[1] is not None and
            self.main_window.loadobservations_widget.velocity_curves[0] is None) or \
                self.data_combobox.currentIndex() == 1:

            obs_color = constants.COLOR_GREEN
            mdl_color = constants.COLOR_ORANGE

        if self.uselc_chk.isChecked():
            x_index = None
            y_index = None

            data = curve.get_data()

            mdl_start = min(data[0])
            mdl_end = max(data[0])

            lc_params = self.main_window.get_lc_params()

            if self.main_window.jdphs_combobox.currentText() == "Time":
                x_index = 0

                lc_params["jdphs"] = 1
                lc_params["hjdst"] = float(mdl_start)
                lc_params["hjdsp"] = float(mdl_end)
                lc_params["hjdin"] = float(self.main_window.p0_ipt.value()) / 500.0

            if self.forcephase_chk.isChecked() or self.main_window.jdphs_combobox.currentText() == "Phase":
                x_index = 1

                lc_params["jdphs"] = 2
                lc_params["phstrt"] = 0.0
                lc_params["phstop"] = 1.0
                lc_params["phin"] = 0.001

            if self.main_window.maglite_combobox.currentText() == "Magnitude":
                y_index = 8

            elif self.main_window.maglite_combobox.currentText() == "Flux":
                y_index = 4

            results = None

            if curve.curve_type == "light":
                lc_params.set_synthetic_curve(
                    curve.band_id,
                    curve.l1,  # l1
                    curve.l2,  # l2
                    curve.x1,  # x1
                    curve.x2,  # x2
                    curve.y1,  # y1
                    curve.y2,  # y2
                    curve.l3,  # l3
                    curve.opsf,  # opsf
                    8.0,  # zero
                    1.0,  # factor
                    0.55,  # wl, dummy
                    curve.aextinc,  # aextinc
                    curve.calib  # calib
                )

                lc_io = wd_io.LCIO(lc_params,
                                   wd_path=self.main_window.lc_path,
                                   lc_binary_name=self.main_window.lc_binary)

                results = lc_io.fill_for_synthetic_light_curve().save().run().read_synthetic_light_curve()

            elif curve.curve_type == "velocity":
                lc_params.set_dummy_synthetic_curve()

                lc_io = wd_io.LCIO(lc_params,
                                   wd_path=self.main_window.lc_path,
                                   lc_binary_name=self.main_window.lc_binary)

                results = lc_io.fill_for_synthetic_velocity_curve().save().run().read_synthetic_velocity_curve()

                y_index = 6

                if (self.data_combobox.currentIndex() == 0 and
                    self.main_window.loadobservations_widget.velocity_curves[1] is not None and
                    self.main_window.loadobservations_widget.velocity_curves[0] is None) or \
                        self.data_combobox.currentIndex() == 1:

                    y_index = 7

            mdl_x = results[x_index]
            mdl_y = results[y_index]

        else:
            mdl_x = obs_x
            mdl_y = mdl_data[dcout_mdl_y_index]

        if curve.curve_type == "light" and self.main_window.maglite_combobox.currentText() == "Magnitude":
            _ = [-2.5 * numpy.log10(x) for x in obs_y]
            obs_y = _
            if not self.uselc_chk.isChecked():
                _ = [-2.5 * numpy.log10(x) for x in mdl_y]
                mdl_y = _
                residuals = numpy.array(obs_y) - numpy.array(mdl_y)
            else:
                _ = [-2.5 * numpy.log10(x) for x in mdl_data[dcout_mdl_y_index]]
                residuals = numpy.array(obs_y) - numpy.array(_)

        self.chart.plot(obs_x, obs_y, clear=False,
                        markersize=constants.MARKER_SIZE, marker="o", linestyle="", color=obs_color)

        if self.uselc_chk.isChecked():
            self.chart.plot(mdl_x, mdl_y, clear=False, color=mdl_color)

        else:
            self.chart.plot(obs_x, mdl_y, clear=False,
                            markersize=constants.MARKER_SIZE - 1, marker="o",
                            linestyle="", color=mdl_color)

        self.chart.plot(obs_x, residuals, index=1,
                        markersize=constants.MARKER_SIZE, marker="o", linestyle="", color=obs_color)

        self.chart.axes[1].axhline([0], color=mdl_color)

        if curve.curve_type == "light" and self.main_window.maglite_combobox.currentText() == "Magnitude":
            self.chart.axes[0].invert_yaxis()
            self.chart.axes[1].invert_yaxis()

        if curve.curve_type == "light":
            self.chart.set_labels("", self.main_window.maglite_combobox.currentText(), index=0)

        elif curve.curve_type == "velocity":
            self.chart.set_labels("", "Km/s", index=0)
            self.chart.axes[0].axhline([self.main_window.vgamma_ipt.value()], color="black", linestyle="--")

        if self.forcephase_chk.isChecked():
            self.chart.set_labels("Phase", "Residuals", index=1)

        else:
            if self.main_window.jdphs_combobox.currentText() == "Time":
                self.chart.set_labels("HJD", "Residuals", index=1)
            else:
                self.chart.set_labels("Phase", "Residuals", index=1)

        self.chart.redraw()

    def update_curve_list(self):
        self.data_combobox.clear()

        curves = self.main_window.loadobservations_widget.get_all_curves()
        for curve in curves:
            self.data_combobox.addItem(os.path.basename(curve.filepath))

        vc1 = self.main_window.loadobservations_widget.velocity_curves[0]
        vc2 = self.main_window.loadobservations_widget.velocity_curves[1]
        if vc1 is not None and vc2 is not None:
            self.data_combobox.addItem("Velocity Curve #1 + #2")

        if len(curves) > 0:
            self.data_combobox.setCurrentIndex(0)

    def update_abort_button_label(self):
        self.rundc2015_btn.setText("Abort (Iteration {0} of {1})".format(self.current_iteration,
                                                                         self.iteration_spinbox.value()))

    def clear_iterator(self):
        if self.iterator is not None:
            self.iterator.stop()
            self.iterator.iteration_done.disconnect()
            self.iterator = None

    def prepare_iterator(self):
        self.clear_iterator()

        dc_params = self.get_dc_params()

        dc_io = wd_io.DCIO(dc_params,
                           wd_path=self.main_window.dc_path,
                           dc_binary_name=self.main_window.dc_binary)

        self.iterator = IterationWorker(dc_io)
        self.iterator.iteration_done.connect(self.receive_iteration)

    def start_differential_correction(self):
        self.lock_ui()

        self.update_abort_button_label()

        self.rundc2015_btn.disconnect()
        self.rundc2015_btn.clicked.connect(self.abort_differential_correction)

        # self.result_treewidget.clear()
        # self.residual_treewidget.clear()
        # self.curvestat_treewidget.clear()
        # self.component_treewidget.clear()

        # self.results = None
        # self.dimensions = None
        # self.observations = None
        # self.solution_stats = None
        # self.curve_stats = None

        # self.data_combobox.clear()

        # self.chart.clear_all()

        self.prepare_iterator()
        self.iterator.start()

    def abort_differential_correction(self):
        self.clear_iterator()

        self.rundc2015_btn.disconnect()
        self.rundc2015_btn.clicked.connect(self.start_differential_correction)
        self.rundc2015_btn.setText("Run DC")

        self.current_iteration = 1

        self.release_ui()

    def update_result_treewidget(self):
        def _format(idx, itm, val, par_id):
            if math.isnan(val) is not True:
                if par_id in (19.0, 20.0):
                    val = val * 10000.0
                itm.setText(idx, str(val))

            else:
                itm.setText(idx, "NaN")
                itm.setBackground(idx, QtGui.QBrush(QtGui.QColor("red")))

        def _populate_item(itm, rslt):
            par_id = rslt[0]
            ipt = rslt[2]
            correction = rslt[3]
            output = rslt[4]
            stderr = rslt[5]

            _format(1, itm, ipt, par_id)
            _format(2, itm, correction, par_id)
            _format(3, itm, output, par_id)
            _format(4, itm, stderr, par_id)

            if numpy.absolute(float(stderr)) > numpy.absolute(float(correction)) or \
                    (numpy.absolute(float(stderr)) == 0.0 and numpy.absolute(float(correction)) == 0.0):
                itm.setBackground(3, QtGui.QBrush(QtGui.QColor("green")))

        if self.results is not None:
            self.result_treewidget.clear()
            for result in self.results:
                param_id = result[0]
                item = QtWidgets.QTreeWidgetItem()
                root = self.result_treewidget.invisibleRootItem()
                if param_id in (56.0, 57.0, 58.0, 59.0, 60.0):
                    curve_name = constants.ID_BANDPASS_DICT[
                        str(self.main_window.loadobservations_widget.light_curves[int(result[1]) - 1].band_id)
                    ]

                    curve_tooltip = os.path.basename(
                        self.main_window.loadobservations_widget.light_curves[int(result[1]) - 1].filepath
                    )

                    found_parent = False
                    i = 0
                    while i < root.childCount():
                        child = root.child(i)
                        if child.toolTip(0) == curve_tooltip:
                            child.addChild(item)
                            item.setText(0, constants.KEEPS_ID_NAME_DICT[result[0]])
                            found_parent = True
                            break
                        i = i + 1

                    if found_parent is False:
                        child = QtWidgets.QTreeWidgetItem(self.result_treewidget)
                        child.setText(0, curve_name)
                        child.setToolTip(0, curve_tooltip)
                        child.addChild(item)
                        child.setExpanded(True)
                        item.setText(0, constants.KEEPS_ID_NAME_DICT[result[0]])
                else:
                    root.addChild(item)
                    item.setText(0, constants.KEEPS_ID_NAME_DICT[result[0]])

                _populate_item(item, result)

    def update_dimension_treewidget(self):
        if self.dimensions is not None:
            self.component_treewidget.clear()
            s1_root = QtWidgets.QTreeWidgetItem(self.component_treewidget)
            s1_root.setText(0, "Star 1")

            s2_root = QtWidgets.QTreeWidgetItem(self.component_treewidget)
            s2_root.setText(0, "Star 2")

            temp_dimensions = list(self.dimensions)

            temp_dimensions[0].pop(3)
            temp_dimensions[0].pop(3)

            temp_dimensions[1].pop(3)
            temp_dimensions[1].pop(3)

            s1_dimensions = list(numpy.array(temp_dimensions[0][1:]).transpose())
            s2_dimensions = list(numpy.array(temp_dimensions[1][1:]).transpose())

            for row in s1_dimensions:
                item = QtWidgets.QTreeWidgetItem(s1_root)
                for index, val in enumerate(row):
                    item.setText(index, val)
                    if val == "nan":
                        item.setText(index, "NaN")
                        item.setBackground(index, QtGui.QBrush(QtGui.QColor("red")))

            for row in s2_dimensions:
                item = QtWidgets.QTreeWidgetItem(s2_root)
                for index, val in enumerate(row):
                    item.setText(index, val)
                    if val == "nan":
                        item.setText(index, "NaN")
                        item.setBackground(index, QtGui.QBrush(QtGui.QColor("red")))

            s1_root.setExpanded(True)
            s2_root.setExpanded(True)

    def update_solution_stat_treewidget(self):
        if self.solution_stats is not None:
            self.residual_treewidget.clear()
            item = QtWidgets.QTreeWidgetItem(self.residual_treewidget)
            a = str(self.solution_stats[0][0])
            b = str(self.solution_stats[1][0])
            c = str(self.solution_stats[2][0])
            d = str(runtime)

            for index, val in enumerate((a, b, c, d)):
                if val == "nan":
                    item.setBackground(index, QtGui.QBrush(QtGui.QColor("red")))
                    val = "NaN"
                item.setText(index, val)

    def update_curve_stat_treewidget(self):
        if self.curve_stats is not None:
            self.curvestat_treewidget.clear()
            curves = self.main_window.loadobservations_widget.get_all_curves()
            for index, row in enumerate(self.curve_stats):
                item = QtWidgets.QTreeWidgetItem(self.curvestat_treewidget)
                item.setText(0, os.path.basename(curves[index].filepath))
                item.setText(1, str(int(row[1])))
                item.setText(2, str(row[2]))
                item.setText(3, str(int(row[3])))
                item.setText(4, str(row[4]))

    def read_dcout(self):
        self.results = self.iterator.dc_io.read_results()
        self.dimensions = self.iterator.dc_io.read_component_dimensions()
        self.observations = self.iterator.dc_io.read_unweighted_observations(split_by_observation=True)
        self.solution_stats = self.iterator.dc_io.read_solution_stats()
        self.curve_stats = self.iterator.dc_io._read_table(self.iterator.dc_io._get_output_path(),
                                                           "    Curve   No. of obs.        Std. dev.",
                                                           tidy=False)

    def receive_iteration(self):
        _template_msg = "\n\nThis is likely caused by malformed dcout file. Check your inputs, KEEP's and observations."
        error = False
        try:
            self.read_dcout()

        except ValueError as e:
            msg = messenger.Messenger("error", "An ValueError has occured:")
            msg.set_info(e.args[0] + _template_msg)
            msg.show()

            error = True

        except IndexError as e:
            msg = messenger.Messenger("error", "An IndexError has occured:")
            msg.set_info(e.args[0] + _template_msg)
            msg.show()

            error = True

        if error:
            self.abort_differential_correction()

            self.results = None
            self.dimensions = None
            self.observations = None
            self.solution_stats = None
            self.curve_stats = None

            self.data_combobox.clear()

            return 0

        self.update_result_treewidget()
        self.update_solution_stat_treewidget()
        self.update_dimension_treewidget()
        self.update_curve_stat_treewidget()

        for result in self.results:
            for value in result:
                if math.isnan(value):
                    msg = messenger.Messenger("warning", "NaN's are encountered while parsing output.")
                    msg.set_info("Please check your inputs and KEEP's.")
                    msg.show()

                    self.abort_differential_correction()

                    self.results = float("nan")
                    self.dimensions = None
                    self.observations = None
                    self.solution_stats = None
                    self.curve_stats = None

                    self.data_combobox.clear()

                    return 0

        if self.autoupdate_chk.isChecked():
            self.plot()

        self.current_iteration = self.current_iteration + 1

        self.main_window.dc_history_widget.receive_solution(self.results, self.solution_stats)

        if self.current_iteration > self.iteration_spinbox.value():
            self.abort_differential_correction()

        else:
            self.update_abort_button_label()
            self.update_inputs()
            self.prepare_iterator()
            self.iterator.start()

    def lock_ui(self):
        self.niter_spinbox.setDisabled(True)
        self.updateinputs_btn.setDisabled(True)
        self.exportresults_btn.setDisabled(True)

        self.keepdel_tab_widget.setDisabled(True)
        self.chart.export_button.setDisabled(True)

        self.main_window.setDisabled(True)

    def release_ui(self):
        self.niter_spinbox.setDisabled(False)
        self.updateinputs_btn.setDisabled(False)
        self.exportresults_btn.setDisabled(False)

        self.keepdel_tab_widget.setDisabled(False)
        self.chart.export_button.setDisabled(False)

        self.main_window.setDisabled(False)

    def write_into_parser(self, parser):
        parser.add_section(constants.CONFIG_SECTION_KEEPS)
        parser.set(constants.CONFIG_SECTION_KEEPS, "xlamda", str(self.marqmul_ipt.value()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "vlr", str(self.vlr_ipt.value()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "t0", str(self.jd0_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "p", str(self.p0_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "dpdt", str(self.dpdt_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "omega", str(self.perr0_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "dperdt", str(self.dperdt_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "pshift", str(self.pshift_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "a", str(self.a_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "e", str(self.e_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "incl", str(self.incl_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "q", str(self.q_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "logd", str(self.logd_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "vgam", str(self.vgam_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "t1", str(self.t1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "t2", str(self.t2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "g1", str(self.g1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "g2", str(self.g2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "alb1", str(self.alb1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "alb2", str(self.alb2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "f1", str(self.f1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "f2", str(self.f2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "l1", str(self.l1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "l2", str(self.l2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "x1", str(self.x1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "x2", str(self.x2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "pot1", str(self.pot1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "pot2", str(self.pot2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "a3b", str(self.a3b_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "p3b", str(self.p3b_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "i3b", str(self.xinc3b_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "e3b", str(self.e3b_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "t03b", str(self.tc3b_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "omega3b", str(self.perr3b_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "l3", str(self.el3_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "desextinc", str(self.desextinc_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1lat", str(self.s1lat_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1lon", str(self.s1long_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1rad", str(self.s1rad_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1temp", str(self.s1temp_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2lat", str(self.s2lat_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2lon", str(self.s2long_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2rad", str(self.s2rad_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2temp", str(self.s2temp_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1tstart", str(self.s1tstart_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1tmax1", str(self.s1tmax1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1tmax2", str(self.s1tmax2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s1tend", str(self.s1tend_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2tstart", str(self.s2tstart_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2tmax1", str(self.s2tmax1_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2tmax2", str(self.s2tmax2_chk.isChecked()))
        parser.set(constants.CONFIG_SECTION_KEEPS, "s2tend", str(self.s2tend_chk.isChecked()))

        parser.add_section(constants.CONFIG_SECTION_DELS)
        parser.set(constants.CONFIG_SECTION_DELS, "a", str(self.del_a_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "e", str(self.del_e_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "incl", str(self.del_i_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "q", str(self.del_q_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "pshift", str(self.del_pshift_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "omega", str(self.del_perr0_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "pot1", str(self.del_pot1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "pot2", str(self.del_pot2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "l1", str(self.del_l1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "l2", str(self.del_l2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "g1", str(self.del_g1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "g2", str(self.del_g2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "t1", str(self.del_t1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "t2", str(self.del_t2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "alb1", str(self.del_alb1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "alb2", str(self.del_alb2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "f1", str(self.del_f1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "f2", str(self.del_f2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "x1", str(self.del_x1_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "x2", str(self.del_x2_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s1lat", str(self.del_s1lat_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s1lon", str(self.del_s1lng_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s1temp", str(self.del_s1tmpf_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s1rad", str(self.del_s1agrad_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s2lat", str(self.del_s2lat_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s2lon", str(self.del_s2lng_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s2temp", str(self.del_s2tmpf_ipt.value()))
        parser.set(constants.CONFIG_SECTION_DELS, "s2rad", str(self.del_s2agrad_ipt.value()))

        return parser

    def read_from_parser(self, parser):
        self.marqmul_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_KEEPS, "xlamda"))
        self.vlr_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_KEEPS, "vlr"))
        self.jd0_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "t0"))
        self.p0_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "p"))
        self.dpdt_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "dpdt"))
        self.perr0_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "omega"))
        self.dperdt_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "dperdt"))
        self.pshift_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "pshift"))
        self.a_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "a"))
        self.e_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "e"))
        self.incl_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "incl"))
        self.q_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "q"))
        self.logd_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "logd"))
        self.vgam_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "vgam"))
        self.t1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "t1"))
        self.t2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "t2"))
        self.g1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "g1"))
        self.g2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "g2"))
        self.alb1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "alb1"))
        self.alb2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "alb2"))
        self.f1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "f1"))
        self.f2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "f2"))
        self.l1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "l1"))
        self.l2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "l2"))
        self.x1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "x1"))
        self.x2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "x2"))
        self.pot1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "pot1"))
        self.pot2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "pot2"))
        self.a3b_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "a3b"))
        self.p3b_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "p3b"))
        self.xinc3b_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "i3b"))
        self.e3b_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "e3b"))
        self.tc3b_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "t03b"))
        self.perr3b_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "omega3b"))
        self.el3_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "l3"))
        self.desextinc_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "desextinc"))
        self.s1lat_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1lat"))
        self.s1long_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1lon"))
        self.s1rad_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1rad"))
        self.s1temp_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1temp"))
        self.s2lat_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2lat"))
        self.s2long_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2lon"))
        self.s2rad_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2rad"))
        self.s2temp_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2temp"))
        self.s1tstart_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1tstart"))
        self.s1tmax1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1tmax1"))
        self.s1tmax2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1tmax2"))
        self.s1tend_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s1tend"))
        self.s2tstart_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2tstart"))
        self.s2tmax1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2tmax1"))
        self.s2tmax2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2tmax2"))
        self.s2tend_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_KEEPS, "s2tend"))

        self.del_a_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "a"))
        self.del_e_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "e"))
        self.del_i_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "incl"))
        self.del_q_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "q"))
        self.del_pshift_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "pshift"))
        self.del_perr0_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "omega"))
        self.del_pot1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "pot1"))
        self.del_pot2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "pot2"))
        self.del_l1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "l1"))
        self.del_l2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "l2"))
        self.del_g1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "g1"))
        self.del_g2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "g2"))
        self.del_t1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "t1"))
        self.del_t2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "t2"))
        self.del_alb1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "alb1"))
        self.del_alb2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "alb2"))
        self.del_f1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "f1"))
        self.del_f2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "f2"))
        self.del_x1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "x1"))
        self.del_x2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "x2"))
        self.del_s1lat_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s1lat"))
        self.del_s1lng_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s1lon"))
        self.del_s1tmpf_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s1temp"))
        self.del_s1agrad_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s1rad"))
        self.del_s2lat_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s2lat"))
        self.del_s2lng_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s2lon"))
        self.del_s2tmpf_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s2temp"))
        self.del_s2agrad_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DELS, "s2rad"))


class IterationWorker(QtCore.QThread):
    iteration_done = QtCore.pyqtSignal(name="done")

    def __init__(self, dc_io):
        super(IterationWorker, self).__init__()

        self.dc_io = dc_io

    def run(self):
        global runtime
        starttime = timeit.default_timer()
        self.dc_io.fill_for_solution().save().run()
        self.iteration_done.emit()
        runtime = (timeit.default_timer() - starttime)
        if runtime < 120.0:
            runtime = ("- %0.2f second -" % (timeit.default_timer() - starttime))
        else:
            runtime = (timeit.default_timer() - starttime)/60.0
            runtime = ("- %0.2f minutes -" % ( runtime ) )


    def stop(self):
        if self.dc_io.process is not None:
            self.blockSignals(True)
            self.dc_io.process.kill()
