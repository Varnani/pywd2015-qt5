from PyQt5 import QtWidgets, QtGui
from gui import mainwindow_widget
from src import constants
from src.helpers import methods, messenger
from src.helpers.wd_utils import wd_containers
from configparser import ConfigParser
import os
import webbrowser  # sending open file commands with this might not be portable
import numpy
from . import lcdcpicker_interface
from . import loadobservations_interface
from . import syntheticcurve_interface
from . import configurespots_interface
from . import eclipsetimings_interface
from . import lineprofile_interface
from . import dimension_interface
from . import starpositions_interface
from . import conjunction_interface
from . import oc_interface
from . import dc_interface
from . import history_interface
from . import single_conjunction


class Widget(QtWidgets.QMainWindow, mainwindow_widget.Ui_MainWindow):
    def __init__(self, app):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.app = app

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.lc_path = None
        self.dc_path = None
        self.lc_binary = None
        self.dc_binary = None

        # create mono font
        if os.path.isfile(constants.MONO_FONT_PATH):
            fid = QtGui.QFontDatabase.addApplicationFont(constants.MONO_FONT_PATH)
            self.monoFont = QtGui.QFont(QtGui.QFontDatabase.applicationFontFamilies(fid)[0])

        else:
            self.monoFont = QtGui.QFont()  # this creates a copy of app's default font

        # interfaces
        self.loadobservations_widget = loadobservations_interface.Widget(self)
        self.configurespot_widget = configurespots_interface.Widget(self)
        self.eclipsetimings_widget = eclipsetimings_interface.Widget(self)
        self.lc_synthetic_curve_widget = syntheticcurve_interface.Widget(self)
        self.lc_lineprofile_widget = lineprofile_interface.Widget(self)
        self.lc_dimension_widget = dimension_interface.Widget(self)
        self.lc_starpositions_widget = starpositions_interface.Widget(self)
        self.lc_conjunction_widget = conjunction_interface.Widget(self)
        self.lc_oc_widget = oc_interface.Widget(self)
        self.single_conjunction_widget = single_conjunction.Widget(self)

        self.dc_widget = dc_interface.Widget(self)
        self.dc_history_widget = history_interface.Widget(self)

        self.populate_styles()
        self.apply_constraints()
        self.check_jdphs()

        self.connect_signals()

    def connect_signals(self):
        self.whatsthis_btn.clicked.connect(QtWidgets.QWhatsThis.enterWhatsThisMode)
        self.systemname_tipt.textChanged.connect(self.system_name_changed)
        self.theme_combobox.currentIndexChanged.connect(self.change_style)
        self.saveproject_btn.clicked.connect(self.save_project)
        self.loadproject_btn.clicked.connect(self.load_project)
        self.fill_btn.clicked.connect(self.fill_button_clicked)
        self.conjunctions_btn.clicked.connect(self.single_conjunction_widget.show)
        self.lc_defaultseed_btn.clicked.connect(self.set_default_seed)
        self.lc_randomizeseed_btn.clicked.connect(self.get_random_seed)

        # tools
        self.t_calc_pot_btn.clicked.connect(self.compute_omega)
        self.t_color_calc_hot_btn.clicked.connect(self.compute_optical_color_temp)
        self.t_color_calc_cool_btn.clicked.connect(self.compute_infrared_color_temp)
        self.t_time_jd_convert_btn.clicked.connect(self.compute_ut)
        self.t_time_ut_convert_btn.clicked.connect(self.compute_jd)

        # show's
        self.loadwidget_btn.clicked.connect(self.loadobservations_widget.show)
        self.spotwidget_btn.clicked.connect(self.configurespot_widget.show)
        self.eclipsewidget_btn.clicked.connect(self.eclipsetimings_widget.show)
        self.lc_lightcurve_btn.clicked.connect(self.lc_synthetic_curve_widget.show)
        self.lc_speclineprof_btn.clicked.connect(self.lc_lineprofile_widget.show)
        self.lc_stardimphase_btn.clicked.connect(self.lc_dimension_widget.show)
        self.lc_coordinates_btn.clicked.connect(self.lc_starpositions_widget.show)
        self.lc_conjunction_btn.clicked.connect(self.lc_conjunction_widget.show)
        self.lc_oc_btn.clicked.connect(self.lc_oc_widget.show)

        self.dc_rundc_btn.clicked.connect(self.dc_widget.show)
        self.dc_history_btn.clicked.connect(self.dc_history_widget.show)

        self.lc_lcin_btn.clicked.connect(self.open_lcin)
        self.lc_lcout_btn.clicked.connect(self.open_lcout)

        # constraints
        self.mode_combobox.currentIndexChanged.connect(self.apply_constraints)
        self.pot1_ipt.valueChanged.connect(self.update_input_pairs)
        self.t1_ipt.valueChanged.connect(self.update_input_pairs)
        self.gr1_ipt.valueChanged.connect(self.update_input_pairs)
        self.alb1_ipt.valueChanged.connect(self.update_input_pairs)
        self.jdphs_combobox.currentIndexChanged.connect(self.check_jdphs)
        self.ld1_chk.stateChanged.connect(self.check_ld)
        self.ld2_chk.stateChanged.connect(self.check_ld)
        self.ipb_chk.stateChanged.connect(self.check_ipb)
        self.e_ipt.valueChanged.connect(self.update_potentials)
        self.q_ipt.valueChanged.connect(self.update_potentials)

    def closeEvent(self, *args, **kwargs):
        self.loadobservations_widget.close()
        self.configurespot_widget.close()
        self.eclipsetimings_widget.close()
        self.single_conjunction_widget.close()
        self.lc_synthetic_curve_widget.close()
        self.lc_lineprofile_widget.close()
        self.lc_dimension_widget.close()
        self.lc_starpositions_widget.close()
        self.lc_conjunction_widget.close()
        self.lc_oc_widget.close()

        self.dc_widget.close()
        self.dc_history_widget.close()

    def open_lcin(self):
        webbrowser.open(os.path.join(self.lc_path, "lcin.active"))

    def open_lcout(self):
        webbrowser.open(os.path.join(self.lc_path, "lcout.active"))

    def set_default_seed(self):
        self.lc_seed_ipt.setValue(138472375.0)

    def get_random_seed(self):
        self.lc_seed_ipt.setValue(float(numpy.random.randint(100000000, high=999999999, dtype=int)))

    def clear_children(self):
        self.lc_synthetic_curve_widget.reset_and_repopulate()
        self.lc_lineprofile_widget.clear()
        self.lc_dimension_widget.clear()
        self.lc_oc_widget.clear()
        self.dc_widget.clear()
        self.dc_history_widget.clear()

    def start_up(self):
        lcdcpicker = lcdcpicker_interface.Widget()
        result = lcdcpicker.check_config()

        if result == QtWidgets.QDialog.Accepted:
            self.lc_path = os.path.dirname(lcdcpicker.lcpath_label.text())
            self.dc_path = os.path.dirname(lcdcpicker.dcpath_label.text())
            self.lc_binary = os.path.basename(lcdcpicker.lcpath_label.text())
            self.dc_binary = os.path.basename(lcdcpicker.dcpath_label.text())

            return True

        else:
            return False

    def system_name_changed(self):
        if self.systemname_tipt.text() == "":
            self.setWindowTitle("PyWD2015")
        else:
            self.setWindowTitle("PyWD2015 - " + self.systemname_tipt.text())

    def fill_button_clicked(self):
        menu = QtWidgets.QMenu()
        curves = self.loadobservations_widget.get_all_curves()

        if len(curves) != 0:

            for index, curve in enumerate(curves):
                action = menu.addAction(os.path.basename(curve.filepath))
                action.idx = index

            selection = menu.exec_(QtGui.QCursor.pos())
            if selection is not None:
                curve = curves[selection.idx]
                data = curve.get_data()
                start = min(data[0])
                end = max(data[0])

                if self.jd_radiobtn.isChecked():
                    self.lc_jd_start_ipt.setValue(start)
                    self.lc_jd_end_ipt.setValue(end)
                    self.lc_jd_incr_ipt.setValue(self.p0_ipt.value() / 1000.0)

                if self.phase_radiobtn.isChecked():
                    self.lc_phs_start_ipt.setValue(start)
                    self.lc_phs_stop_ipt.setValue(end)
                    self.lc_phs_incr_ipt.setValue(0.001)

    def populate_styles(self):
        style_factory = QtWidgets.QStyleFactory()
        style_keys = list(style_factory.keys())

        for style in style_keys:
            self.theme_combobox.addItem(str(style))

    def change_style(self):
        font = self.font()
        style_factory = QtWidgets.QStyleFactory()
        style = style_factory.create(self.theme_combobox.currentText())
        self.app.setStyle(style)
        self.app.setFont(font)

    def curve_list_changed(self):
        self.lc_synthetic_curve_widget.reset_and_repopulate()
        self.dc_widget.clear()

    def compute_jd(self):
        year = self.t_time_year_ipt.value()
        month = self.t_time_month_ipt.value()
        day = self.t_time_day_ipt.value()
        hour = self.t_time_hour_ipt.value()
        minute = self.t_time_minute_ipt.value()
        second = self.t_time_second_ipt.value()

        jd = methods.convert_ut_to_jd(year, month, day, hour, minute, second)

        self.t_time_jd_otpt.setValue(jd)

    def compute_ut(self):
        year, month, day, hour, minute, second = methods.convert_jd_to_ut(self.t_time_jd_ipt.value(), add_24=False)
        self.t_time_dmy_otpt.setText(str(day) + "/" + str(month) + "/" + str(year) + " - " +
                                     str(hour) + "/" + str(minute) + "/" + str(numpy.round(second, decimals=3)))

    def compute_optical_color_temp(self):
        color = self.t_color_bv_ipt.value()
        error = self.t_color_bv_err_ipt.value()

        try:

            gray, gray_error = methods.compute_temp_from_color("gray", color, error)
            flower, flower_error = methods.compute_temp_from_color("flower", color, error)
            dl, dl_error = methods.compute_temp_from_color("drilling_landolt", color, error)
            popper, popper_error = methods.compute_temp_from_color("popper", color, error)

            self.t_color_gray_otpt.setValue(gray)
            self.t_color_gray_err_otpt.setValue(gray_error)

            self.t_color_flower_otpt.setValue(flower)
            self.t_color_flower_err_otpt.setValue(flower_error)

            self.t_color_dl_otpt.setValue(dl)
            self.t_color_dl_err_otpt.setValue(dl_error)

            self.t_color_popper_otpt.setValue(popper)
            self.t_color_popper_err_otpt.setValue(popper_error)

        except OverflowError as e:
            msg = messenger.Messenger("error", "An error encountered. Check your inputs.")
            msg.set_info("Details:\n" + e.args[0])
            msg.show()

    def compute_infrared_color_temp(self):
        color = self.t_color_cool_ipt.value()
        error = self.t_color_cool_err_ipt.value()

        ref_dict = {"V - K": "tokunaga_vk",
                    "J - H": "tokunaga_jh",
                    "H - K": "tokunaga_hk"}

        try:

            temp, error = methods.compute_temp_from_color(ref_dict[self.t_color_cool_combobox.currentText()],
                                                          color, error)

            self.t_color_tokunaga_otpt.setValue(temp)
            self.t_color_tokunaga_err_otpt.setValue(error)

        except OverflowError as e:
            msg = messenger.Messenger("error", "An error encountered. Check your inputs.")
            msg.set_info("Details:\n" + e.args[0])
            msg.show()

    def compute_omega(self):
        try:
            pot = methods.compute_omega_potential(self.t_pot_q_ipt.value(), self.t_pot_rad_ipt.value(),
                                                  self.t_pot_f_ipt.value(), self.t_pot_d_ipt.value())
            self.t_pot_pot_otpt.setValue(pot)

        except:
            self.t_pot_pot_otpt.setValue(0.0)

    def save_project(self):
        save_path = methods.save_file(self, suffix=".pywdproject", name_filter="PyWD2015 project file (*.pywdproject)")
        if save_path is not None:
            parser = ConfigParser()

            parser.add_section(constants.CONFIG_SECTION_INFO)
            parser.set(constants.CONFIG_SECTION_INFO, "version", constants.MAIN_VERSION)
            parser.set(constants.CONFIG_SECTION_INFO, "config", constants.CONFIG_VERSION)

            parser.add_section(constants.CONFIG_SECTION_MAIN)
            parser.set(constants.CONFIG_SECTION_MAIN, "system_name", self.systemname_tipt.text())
            parser.set(constants.CONFIG_SECTION_MAIN, "operation_mode", str(self.mode_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_MAIN, "jdphs", str(self.jdphs_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_MAIN, "maglite", str(self.maglite_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_MAIN, "rv_corr_1", str(self.icor1_chk.isChecked()))
            parser.set(constants.CONFIG_SECTION_MAIN, "rv_corr_2", str(self.icor2_chk.isChecked()))
            parser.set(constants.CONFIG_SECTION_MAIN, "ifcgs", str(self.ifcgs_chk.isChecked()))

            parser.add_section(constants.CONFIG_SECTION_SYSTEM)
            parser.set(constants.CONFIG_SECTION_SYSTEM, "jd0", str(self.jd0_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "period", str(self.p0_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "dpdt", str(self.dpdt_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "ecc", str(self.e_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "omega", str(self.perr0_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "dperdt", str(self.dperdt_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "incl", str(self.incl_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "a", str(self.a_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "pshift", str(self.pshift_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "vgam", str(self.vgamma_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "vunit", str(self.vunit_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "logd", str(self.dpclog_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "q", str(self.q_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "abunin", str(self.abunin_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "t1", str(self.t1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "t2", str(self.t2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "pot1", str(self.pot1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "pot2", str(self.pot2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "f1", str(self.f1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "f2", str(self.f2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "the", str(self.the_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "delph", str(self.delph_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SYSTEM, "nga", str(self.nga_ipt.value()))

            parser.add_section(constants.CONFIG_SECTION_SURFACE)
            parser.set(constants.CONFIG_SECTION_SURFACE, "ifat1", str(self.ifat1_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ifat2", str(self.ifat2_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "alb1", str(self.alb1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "alb2", str(self.alb2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "gr1", str(self.gr1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "gr2", str(self.gr2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "n1", str(self.n1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "n2", str(self.n2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "n1l", str(self.n1l_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "n2l", str(self.n2l_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ld1", str(self.ld1_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ld2", str(self.ld2_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ld1_fixed", str(self.ld1_chk.isChecked()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ld2_fixed", str(self.ld2_chk.isChecked()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ipb", str(self.ipb_chk.isChecked()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "mref", str(self.mref_chk.isChecked()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "nref", str(self.nref_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "xbol_s1", str(self.xbol1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "xbol_s2", str(self.xbol2_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ybol_s1", str(self.ybol1_ipt.value()))
            parser.set(constants.CONFIG_SECTION_SURFACE, "ybol_s2", str(self.ybol2_ipt.value()))

            parser.add_section(constants.CONFIG_SECTION_3B)
            parser.set(constants.CONFIG_SECTION_3B, "if3b", str(self.if3b_groupbox.isChecked()))
            parser.set(constants.CONFIG_SECTION_3B, "a_3b", str(self.a3b_ipt.value()))
            parser.set(constants.CONFIG_SECTION_3B, "p_3b", str(self.p3b_ipt.value()))
            parser.set(constants.CONFIG_SECTION_3B, "incl_3b", str(self.incl3b_ipt.value()))
            parser.set(constants.CONFIG_SECTION_3B, "ecc_3b", str(self.e3b_ipt.value()))
            parser.set(constants.CONFIG_SECTION_3B, "omega_3b", str(self.perr03b_ipt.value()))
            parser.set(constants.CONFIG_SECTION_3B, "conj_time", str(self.conj3b_ipt.value()))

            parser.add_section(constants.CONFIG_SECTION_LC2015)
            parser.set(constants.CONFIG_SECTION_LC2015, "std_dev", str(self.lc_fract_std_dev_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "noise", str(self.lc_noise_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_LC2015, "seed", str(self.lc_seed_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "jd_start", str(self.lc_jd_start_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "jd_stop", str(self.lc_jd_end_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "jd_incr", str(self.lc_jd_incr_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "phase_start", str(self.lc_phs_start_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "phase_stop", str(self.lc_phs_stop_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "phase_incr", str(self.lc_phs_incr_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "lsp", str(self.lc_lsp_spinbox.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "tobs", str(self.lc_tobs_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "phobs", str(self.lc_phobs_ipt.value()))
            parser.set(constants.CONFIG_SECTION_LC2015, "phn", str(self.lc_phn_ipt.value()))

            parser.add_section(constants.CONFIG_SECTION_DC2015)
            parser.set(constants.CONFIG_SECTION_DC2015, "deriv_type", str(self.dc_isym_combobox.currentIndex()))
            parser.set(constants.CONFIG_SECTION_DC2015, "desextinc", str(self.dc_desextinc_ipt.value()))
            parser.set(constants.CONFIG_SECTION_DC2015, "linkext", str(self.dc_ext_band_spinbox.value()))

            parser = self.dc_widget.write_into_parser(parser)
            parser = self.configurespot_widget.write_into_parser(parser)
            parser = self.eclipsetimings_widget.write_into_parser(parser)
            parser = self.loadobservations_widget.write_into_parser(parser)

            with open(save_path, "w") as destination:
                parser.write(destination)

            msg = messenger.Messenger("info", "Save completed:")
            msg.set_info(save_path)
            msg.show()

    def load_project(self):
        msg = messenger.Messenger("warning", "Loading a project file will reset all widgets and clear all plots.")
        msg.set_info("Are you sure?")
        msg.msg_box.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)

        if msg.show() == QtWidgets.QMessageBox.Cancel:
            return 1

        load_path = methods.load_file(self, suffix=".pywdproject", name_filter="PyWD2015 project file (*.pywdproject)")
        if load_path is not None:
            parser = ConfigParser()
            with open(load_path, "r") as source:
                parser.read_file(source)

                if constants.MAIN_VERSION != parser.get(constants.CONFIG_SECTION_INFO, "version"):
                    msg = messenger.Messenger("error", "This file is saved with a different version of the program "
                                                       "and cannot be loaded.")
                    msg.set_info("Program Version: " + constants.MAIN_VERSION + "\nFile Version: " +
                                 parser.get(constants.CONFIG_SECTION_INFO, "version"))
                    msg.show()

                    return 1

                if constants.CONFIG_VERSION != parser.get(constants.CONFIG_SECTION_INFO, "config"):
                    msg = messenger.Messenger("error", "This file is saved with a different config version "
                                                       "and cannot be loaded.")
                    msg.set_info("Program Config Version: " + constants.CONFIG_VERSION + "\nFile Config Version: " +
                                 parser.get(constants.CONFIG_SECTION_INFO, "config"))
                    msg.show()

                    return 1

                path_errors = ""

                def _check_for_file(pth, second_dir):
                    if os.path.isfile(pth) is not True:
                        second_path = os.path.join(os.path.dirname(second_dir), os.path.basename(pth))
                        if os.path.isfile(second_path):
                            return second_path
                        else:
                            return None
                    else:
                        return pth

                lc_count = parser.getint(constants.CONFIG_SECTION_CURVE_COUNTS, "light_curves")

                if lc_count > 0:
                    i = 1
                    while i <= lc_count:
                        sct = constants.CONFIG_SECTION_LIGHT_CURVE_STUB + str(i)
                        path = parser.get(sct, "filepath")
                        n_path = _check_for_file(path, load_path)
                        if n_path is None:
                            path_errors = path_errors + path + "\n"

                        else:
                            parser.set(sct, "filepath", n_path)

                        i = i + 1

                vc1 = parser.getboolean(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_1")
                vc2 = parser.getboolean(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_2")

                if vc1:
                    sct = constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + "1"
                    path = parser.get(sct, "filepath")
                    n_path = _check_for_file(path, load_path)
                    if n_path is None:
                        path_errors = path_errors + path + "\n"

                    else:
                        parser.set(sct, "filepath", n_path)

                if vc2:
                    sct = constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + "2"
                    path = parser.get(sct, "filepath")
                    n_path = _check_for_file(path, load_path)
                    if n_path is None:
                        path_errors = path_errors + path + "\n"

                    else:
                        parser.set(sct, "filepath", n_path)

                eclipse_path = parser.get(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "filepath")

                if eclipse_path != "None":
                    n_path = _check_for_file(eclipse_path, load_path)
                    if n_path is None:
                        path_errors = path_errors + eclipse_path + "\n"

                    else:
                        parser.set(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "filepath", n_path)

                if path_errors != "":
                    msg = messenger.Messenger("error", "Some files can not be found:")
                    msg.set_info(path_errors)
                    msg.show()

                    return 1

                self.systemname_tipt.setText(parser.get(constants.CONFIG_SECTION_MAIN, "system_name"))
                self.mode_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_MAIN, "operation_mode"))
                self.jdphs_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_MAIN, "jdphs"))
                self.maglite_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_MAIN, "maglite"))
                self.icor1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_MAIN, "rv_corr_1"))
                self.icor2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_MAIN, "rv_corr_2"))
                self.ifcgs_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_MAIN, "ifcgs"))

                self.jd0_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "jd0"))
                self.p0_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "period"))
                self.dpdt_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "dpdt"))
                self.e_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "ecc"))
                self.perr0_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "omega"))
                self.dperdt_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "dperdt"))
                self.incl_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "incl"))
                self.a_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "a"))
                self.pshift_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "pshift"))
                self.vgamma_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "vgam"))
                self.vunit_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "vunit"))
                self.dpclog_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "logd"))
                self.q_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "q"))
                self.abunin_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "abunin"))
                self.t1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "t1"))
                self.t2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "t2"))
                self.pot1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "pot1"))
                self.pot2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "pot2"))
                self.f1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "f1"))
                self.f2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "f2"))
                self.the_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "the"))
                self.delph_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SYSTEM, "delph"))
                self.nga_ipt.setValue(parser.getint(constants.CONFIG_SECTION_SYSTEM, "nga"))

                self.ifat1_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_SURFACE, "ifat1"))
                self.ifat2_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_SURFACE, "ifat2"))
                self.alb1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "alb1"))
                self.alb2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "alb2"))
                self.gr1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "gr1"))
                self.gr2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "gr2"))
                self.n1_ipt.setValue(parser.getint(constants.CONFIG_SECTION_SURFACE, "n1"))
                self.n2_ipt.setValue(parser.getint(constants.CONFIG_SECTION_SURFACE, "n2"))
                self.n1l_ipt.setValue(parser.getint(constants.CONFIG_SECTION_SURFACE, "n1l"))
                self.n2l_ipt.setValue(parser.getint(constants.CONFIG_SECTION_SURFACE, "n2l"))
                self.ld1_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_SURFACE, "ld1"))
                self.ld2_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_SURFACE, "ld2"))
                self.ld1_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SURFACE, "ld1_fixed"))
                self.ld2_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SURFACE, "ld2_fixed"))
                self.ipb_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SURFACE, "ipb"))
                self.mref_chk.setChecked(parser.getboolean(constants.CONFIG_SECTION_SURFACE, "mref"))
                self.nref_ipt.setValue(parser.getint(constants.CONFIG_SECTION_SURFACE, "nref"))
                self.xbol1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "xbol_s1"))
                self.xbol2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "xbol_s2"))
                self.ybol1_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "ybol_s1"))
                self.ybol2_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_SURFACE, "ybol_s2"))

                self.if3b_groupbox.setChecked(parser.getboolean(constants.CONFIG_SECTION_3B, "if3b"))
                self.a3b_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_3B, "a_3b"))
                self.p3b_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_3B, "p_3b"))
                self.incl3b_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_3B, "incl_3b"))
                self.e3b_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_3B, "ecc_3b"))
                self.perr03b_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_3B, "omega_3b"))
                self.conj3b_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_3B, "conj_time"))

                self.lc_fract_std_dev_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "std_dev"))
                self.lc_noise_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_LC2015, "noise"))
                self.lc_seed_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "seed"))
                self.lc_jd_start_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "jd_start"))
                self.lc_jd_end_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "jd_stop"))
                self.lc_jd_incr_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "jd_incr"))
                self.lc_phs_start_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "phase_start"))
                self.lc_phs_stop_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "phase_stop"))
                self.lc_phs_incr_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "phase_incr"))
                self.lc_lsp_spinbox.setValue(parser.getint(constants.CONFIG_SECTION_LC2015, "lsp"))
                self.lc_tobs_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "tobs"))
                self.lc_phobs_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "phobs"))
                self.lc_phn_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_LC2015, "phn"))

                self.dc_isym_combobox.setCurrentIndex(parser.getint(constants.CONFIG_SECTION_DC2015, "deriv_type"))
                self.dc_desextinc_ipt.setValue(parser.getfloat(constants.CONFIG_SECTION_DC2015, "desextinc"))
                self.dc_ext_band_spinbox.setValue(parser.getint(constants.CONFIG_SECTION_DC2015, "linkext"))

                self.dc_widget.read_from_parser(parser)
                self.configurespot_widget.read_from_parser(parser)
                self.eclipsetimings_widget.read_from_parser(parser)
                self.loadobservations_widget.read_from_parser(parser)

                self.clear_children()

                msg = messenger.Messenger("info", "Load completed:")
                msg.set_info(load_path)
                msg.show()

    @staticmethod
    def evaluate_checkbox(checkbox, reverse=False):
        if checkbox.isChecked():
            if reverse:
                return 0
            else:
                return 1

        else:
            if reverse:
                return 1
            else:
                return 0

    def get_lc_boundaries(self):
        start, end = None, None

        if self.jd_radiobtn.isChecked():
            start = self.lc_jd_start_ipt.value()
            end = self.lc_jd_end_ipt.value()
        elif self.phase_radiobtn.isChecked():
            start = self.lc_phs_start_ipt.value()
            end = self.lc_phs_stop_ipt.value()

        return start, end

    def get_lc_params(self):
        lc_params = wd_containers.LCParameterContainer()

        # general params
        jdphs = 0

        if self.jd_radiobtn.isChecked():
            jdphs = 1

        elif self.phase_radiobtn.isChecked():
            jdphs = 2

        lc_params["jdphs"] = jdphs
        lc_params["ifcgs"] = self.evaluate_checkbox(self.ifcgs_chk)
        lc_params["mode"] = constants.MODE_DICT[self.mode_combobox.currentText()]
        lc_params["icor1"] = self.evaluate_checkbox(self.icor1_chk)
        lc_params["icor2"] = self.evaluate_checkbox(self.icor2_chk)

        # system params
        lc_params["hjd0"] = self.jd0_ipt.value()
        lc_params["pzero"] = self.p0_ipt.value()
        lc_params["dpdt"] = self.dpdt_ipt.value()
        lc_params["pshift"] = self.pshift_ipt.value()
        lc_params["delph"] = self.delph_ipt.value()
        lc_params["nga"] = self.nga_ipt.value()
        lc_params["e"] = self.e_ipt.value()
        lc_params["a"] = self.a_ipt.value()
        lc_params["f1"] = self.f1_ipt.value()
        lc_params["f2"] = self.f2_ipt.value()
        lc_params["vga"] = self.vgamma_ipt.value()
        lc_params["xincl"] = self.incl_ipt.value()
        lc_params["tavh"] = self.t1_ipt.value()
        lc_params["tavc"] = self.t2_ipt.value()
        lc_params["phsv"] = self.pot1_ipt.value()
        lc_params["pcsv"] = self.pot2_ipt.value()
        lc_params["rm"] = self.q_ipt.value()
        lc_params["perr"] = self.perr0_ipt.value()
        lc_params["dperdt"] = self.dperdt_ipt.value()
        lc_params["the"] = self.the_ipt.value()
        lc_params["vunit"] = self.vunit_ipt.value()
        lc_params["abunin"] = self.abunin_ipt.value()
        lc_params["dpclog"] = self.dpclog_ipt.value()

        # surface params
        lc_params["ifat1"] = constants.ATM_DICT[self.ifat1_combobox.currentText()]
        lc_params["ifat2"] = constants.ATM_DICT[self.ifat2_combobox.currentText()]
        lc_params["gr1"] = self.gr1_ipt.value()
        lc_params["gr2"] = self.gr2_ipt.value()
        lc_params["ipb"] = self.evaluate_checkbox(self.ipb_chk)
        lc_params["mref"] = self.evaluate_checkbox(self.mref_chk) + 1
        lc_params["nref"] = self.nref_ipt.value()
        lc_params["n1"] = self.n1_ipt.value()
        lc_params["n2"] = self.n2_ipt.value()
        lc_params["alb1"] = self.alb1_ipt.value()
        lc_params["alb2"] = self.alb2_ipt.value()
        lc_params["xbol1"] = self.xbol1_ipt.value()
        lc_params["xbol2"] = self.xbol2_ipt.value()
        lc_params["ybol1"] = self.ybol1_ipt.value()
        lc_params["ybol2"] = self.ybol2_ipt.value()

        if not self.ld1_chk.isChecked():
            lc_params["ld1"] = constants.LIMB_DICT[self.ld1_combobox.currentText()] * -1
        else:
            lc_params["ld1"] = constants.LIMB_DICT[self.ld1_combobox.currentText()]

        if not self.ld2_chk.isChecked():
            lc_params["ld2"] = constants.LIMB_DICT[self.ld2_combobox.currentText()] * -1
        else:
            lc_params["ld2"] = constants.LIMB_DICT[self.ld2_combobox.currentText()]

        # third body
        lc_params["if3b"] = self.evaluate_checkbox(self.if3b_groupbox)
        lc_params["a3b"] = self.a3b_ipt.value()
        lc_params["p3b"] = self.p3b_ipt.value()
        lc_params["xincl3b"] = self.incl3b_ipt.value()
        lc_params["e3b"] = self.e3b_ipt.value()
        lc_params["perr3b"] = self.perr03b_ipt.value()
        lc_params["tc3b"] = self.conj3b_ipt.value()

        # spot params
        lc_params["nomax"] = constants.NOMAX_DICT[self.configurespot_widget.nomax_combobox.currentText()]
        lc_params["kspev"] = self.evaluate_checkbox(self.configurespot_widget.kspev_groupbox)
        lc_params["kspot"] = self.evaluate_checkbox(self.configurespot_widget.kspot_chk) + 1
        lc_params["fspot1"] = self.configurespot_widget.fspot1_ipt.value()
        lc_params["fspot2"] = self.configurespot_widget.fspot2_ipt.value()
        lc_params["ifsmv1"] = self.evaluate_checkbox(self.configurespot_widget.ifsmv1_chk)
        lc_params["ifsmv2"] = self.evaluate_checkbox(self.configurespot_widget.ifsmv2_chk)

        # noise
        lc_params["stdev"] = self.lc_fract_std_dev_ipt.value()
        lc_params["noise"] = constants.NOISE_DICT[self.lc_noise_combobox.currentText()]
        lc_params["seed"] = int(self.lc_seed_ipt.value())

        # steps
        lc_params["hjdst"] = self.lc_jd_start_ipt.value()
        lc_params["hjdsp"] = self.lc_jd_end_ipt.value()
        lc_params["hjdin"] = self.lc_jd_incr_ipt.value()
        lc_params["phstrt"] = self.lc_phs_start_ipt.value()
        lc_params["phstop"] = self.lc_phs_stop_ipt.value()
        lc_params["phin"] = self.lc_phs_incr_ipt.value()

        # normalization
        lc_params["phn"] = self.lc_phn_ipt.value()

        # temperature estimation params
        lc_params["phobs"] = self.lc_phobs_ipt.value()
        lc_params["lsp"] = self.lc_lsp_spinbox.value()
        lc_params["tobs"] = self.lc_tobs_ipt.value()

        # line profile parameters, default values
        lc_params["binwm1"] = 0.00001
        lc_params["sc1"] = 1.0
        lc_params["sl1"] = 0.0
        lc_params["nf1"] = 1

        lc_params["binwm2"] = 0.00001
        lc_params["sc2"] = 1.0
        lc_params["sl2"] = 0.0
        lc_params["nf2"] = 1

        s1_spots, s2_spots = self.configurespot_widget.get_all_spots()

        for spot in s1_spots:
            lc_params.add_spot(1,
                               spot[0], spot[1], spot[2], spot[3],
                               spot[4], spot[5], spot[6], spot[7])

        for spot in s2_spots:
            lc_params.add_spot(2,
                               spot[0], spot[1], spot[2], spot[3],
                               spot[4], spot[5], spot[6], spot[7])

        return lc_params

    def get_dc_params(self):
        dc_params = wd_containers.DCParameterContainer()

        # general params
        dc_params["jdphs"] = constants.JDPHS_DICT[self.jdphs_combobox.currentText()]
        dc_params["ifcgs"] = self.evaluate_checkbox(self.ifcgs_chk)
        dc_params["mode"] = constants.MODE_DICT[self.mode_combobox.currentText()]
        dc_params["icor1"] = self.evaluate_checkbox(self.icor1_chk)
        dc_params["icor2"] = self.evaluate_checkbox(self.icor2_chk)

        # system params
        dc_params["hjd0"] = self.jd0_ipt.value()
        dc_params["pzero"] = self.p0_ipt.value()
        dc_params["dpdt"] = self.dpdt_ipt.value()
        dc_params["pshift"] = self.pshift_ipt.value()
        dc_params["delph"] = self.delph_ipt.value()
        dc_params["nga"] = self.nga_ipt.value()
        dc_params["e"] = self.e_ipt.value()
        dc_params["a"] = self.a_ipt.value()
        dc_params["f1"] = self.f1_ipt.value()
        dc_params["f2"] = self.f2_ipt.value()
        dc_params["vga"] = self.vgamma_ipt.value()
        dc_params["xincl"] = self.incl_ipt.value()
        dc_params["tavh"] = self.t1_ipt.value()
        dc_params["tavc"] = self.t2_ipt.value()
        dc_params["phsv"] = self.pot1_ipt.value()
        dc_params["pcsv"] = self.pot2_ipt.value()
        dc_params["rm"] = self.q_ipt.value()
        dc_params["perr"] = self.perr0_ipt.value()
        dc_params["dperdt"] = self.dperdt_ipt.value()
        dc_params["the"] = self.the_ipt.value()
        dc_params["vunit"] = self.vunit_ipt.value()
        dc_params["abunin"] = self.abunin_ipt.value()
        dc_params["dpclog"] = self.dpclog_ipt.value()

        # surface params
        dc_params["ifat1"] = constants.ATM_DICT[self.ifat1_combobox.currentText()]
        dc_params["ifat2"] = constants.ATM_DICT[self.ifat2_combobox.currentText()]
        dc_params["gr1"] = self.gr1_ipt.value()
        dc_params["gr2"] = self.gr2_ipt.value()
        dc_params["ipb"] = self.evaluate_checkbox(self.ipb_chk)
        dc_params["mref"] = self.evaluate_checkbox(self.mref_chk) + 1
        dc_params["nref"] = self.nref_ipt.value()
        dc_params["n1"] = self.n1_ipt.value()
        dc_params["n2"] = self.n2_ipt.value()
        dc_params["alb1"] = self.alb1_ipt.value()
        dc_params["alb2"] = self.alb2_ipt.value()
        dc_params["xbol1"] = self.xbol1_ipt.value()
        dc_params["xbol2"] = self.xbol2_ipt.value()
        dc_params["ybol1"] = self.ybol1_ipt.value()
        dc_params["ybol2"] = self.ybol2_ipt.value()

        if not self.ld1_chk.isChecked():
            dc_params["ld1"] = constants.LIMB_DICT[self.ld1_combobox.currentText()] * -1
        else:
            dc_params["ld1"] = constants.LIMB_DICT[self.ld1_combobox.currentText()]

        if not self.ld2_chk.isChecked():
            dc_params["ld2"] = constants.LIMB_DICT[self.ld2_combobox.currentText()] * -1
        else:
            dc_params["ld2"] = constants.LIMB_DICT[self.ld2_combobox.currentText()]

        # third body
        dc_params["if3b"] = self.evaluate_checkbox(self.if3b_groupbox)
        dc_params["a3b"] = self.a3b_ipt.value()
        dc_params["p3b"] = self.p3b_ipt.value()
        dc_params["xincl3b"] = self.incl3b_ipt.value()
        dc_params["e3b"] = self.e3b_ipt.value()
        dc_params["perr3b"] = self.perr03b_ipt.value()
        dc_params["tc3b"] = self.conj3b_ipt.value()

        # general dc params
        dc_params["isym"] = constants.DERIV_DICT[self.dc_isym_combobox.currentText()]
        dc_params["maglite"] = constants.MAGLITE_DICT[self.maglite_combobox.currentText()]
        dc_params["linkext"] = self.dc_ext_band_spinbox.value()
        dc_params["desextinc"] = self.dc_desextinc_ipt.value()
        dc_params["n1l"] = self.n1l_ipt.value()
        dc_params["n2l"] = self.n2l_ipt.value()

        # spot params
        dc_params["nomax"] = constants.NOMAX_DICT[self.configurespot_widget.nomax_combobox.currentText()]
        dc_params["kspev"] = self.evaluate_checkbox(self.configurespot_widget.kspev_groupbox)
        dc_params["kspot"] = self.evaluate_checkbox(self.configurespot_widget.kspot_chk) + 1
        dc_params["fspot1"] = self.configurespot_widget.fspot1_ipt.value()
        dc_params["fspot2"] = self.configurespot_widget.fspot2_ipt.value()
        dc_params["ifsmv1"] = self.evaluate_checkbox(self.configurespot_widget.ifsmv1_chk)
        dc_params["ifsmv2"] = self.evaluate_checkbox(self.configurespot_widget.ifsmv2_chk)

        s1_spots, s2_spots = self.configurespot_widget.get_all_spots(get_ab=True)

        def _parse_spots(spots, params, star):
            for i, spot in enumerate(spots):
                a = spot[8]
                b = spot[9]

                if a:
                    params["kspa"] = star
                    params["nspa"] = i + 1
                if b:
                    params["kspb"] = star
                    params["nspb"] = i + 1

                params.add_spot(star,
                                spot[0], spot[1], spot[2], spot[3],
                                spot[4], spot[5], spot[6], spot[7])

        _parse_spots(s1_spots, dc_params, 1)
        _parse_spots(s2_spots, dc_params, 2)

        for lc in self.loadobservations_widget.light_curves:
            if lc.l2 == 0.0:
                lc.l2 = 1.0
            data = lc.get_data()
            dc_params.add_light_curve(lc.band_id, lc.l1, lc.l2, lc.x1, lc.x2, lc.y1, lc.y2, lc.opsf, lc.sigma, lc.ksd,
                                      lc.l3, lc.noise, lc.aextinc, lc.calib, data[0], data[1], data[2],
                                      xunit=lc.xunit, spha1=lc.e1, spha2=lc.e2, spha3=lc.e3, spha4=lc.e4)

        if self.loadobservations_widget.velocity_curves[0] is not None:
            vc = self.loadobservations_widget.velocity_curves[0]
            data = vc.get_data()
            dc_params.add_velocity_curve(1, vc.sigma, vc.ksd, vc.wla, data[0], data[1], data[2])

        if self.loadobservations_widget.velocity_curves[1] is not None:
            vc = self.loadobservations_widget.velocity_curves[1]
            data = vc.get_data()
            dc_params.add_velocity_curve(2, vc.sigma, vc.ksd, vc.wla, data[0], data[1], data[2])

        if self.eclipsetimings_widget.iftime_chk.isChecked():
            data = self.eclipsetimings_widget.get_data()
            dc_params.add_eclipse_times(self.eclipsetimings_widget.sigma_ipt.value(),
                                        self.eclipsetimings_widget.ksd_box.value(),
                                        data[0], data[1], data[2])

        return dc_params

    def apply_constraints(self):
        self.clear_constraints()
        if str(self.mode_combobox.currentText()) == "Mode -1":
            self.pot2_ipt.setDisabled(True)
            self.pot2_ipt.setValue(self.pot1_ipt.value())
            self.dc_widget.pot2_chk.setDisabled(True)
            self.dc_widget.pot2_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 0":
            self.ipb_chk.setChecked(True)
            self.ipb_chk.setDisabled(True)

        if str(self.mode_combobox.currentText()) == "Mode 1":
            # Curve independent params
            self.pot2_ipt.setValue(self.pot1_ipt.value())
            self.t2_ipt.setValue(self.t1_ipt.value())
            self.gr2_ipt.setValue(self.gr1_ipt.value())
            self.alb2_ipt.setValue(self.alb1_ipt.value())

            self.pot2_ipt.setDisabled(True)
            self.dc_widget.pot2_chk.setDisabled(True)
            self.dc_widget.pot2_chk.setChecked(False)

            self.t2_ipt.setDisabled(True)
            self.dc_widget.t2_chk.setDisabled(True)
            self.dc_widget.t2_chk.setChecked(False)

            self.gr2_ipt.setDisabled(True)
            self.dc_widget.g2_chk.setDisabled(True)
            self.dc_widget.g2_chk.setChecked(False)

            self.alb2_ipt.setDisabled(True)
            self.dc_widget.alb2_chk.setDisabled(True)
            self.dc_widget.alb2_chk.setChecked(False)

            # Curve dependent params
            if self.ipb_chk.isChecked() is not True:
                self.dc_widget.l2_chk.setDisabled(True)
                self.dc_widget.l2_chk.setChecked(False)
                if self.ld1_chk.isChecked() and self.ld2_chk.isChecked():
                    for curve in self.loadobservations_widget.light_curves:
                        curve.x2 = curve.x1
                        curve.y2 = curve.y1
            self.dc_widget.x2_chk.setDisabled(True)
            self.dc_widget.x2_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 3":
            self.pot2_ipt.setValue(self.pot1_ipt.value())
            self.pot2_ipt.setDisabled(True)
            self.dc_widget.pot2_chk.setDisabled(True)
            self.dc_widget.pot2_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 2":
            self.ipb_chk.setDisabled(True)
            self.ipb_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 4":
            self.update_potentials()
            self.pot1_ipt.setDisabled(True)
            self.dc_widget.pot1_chk.setDisabled(True)
            self.dc_widget.pot1_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 5":
            self.update_potentials()
            self.pot2_ipt.setDisabled(True)
            self.dc_widget.pot2_chk.setDisabled(True)
            self.dc_widget.pot2_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 6":
            self.update_potentials()
            self.pot1_ipt.setDisabled(True)
            self.pot2_ipt.setDisabled(True)
            self.dc_widget.pot1_chk.setDisabled(True)
            self.dc_widget.pot2_chk.setDisabled(True)
            self.dc_widget.pot1_chk.setChecked(False)
            self.dc_widget.pot2_chk.setChecked(False)

        self.check_ld()
        self.check_ipb()

    def clear_constraints(self):
        self.pot1_ipt.setDisabled(False)
        self.pot2_ipt.setDisabled(False)
        self.dc_widget.pot1_chk.setDisabled(False)
        self.dc_widget.pot2_chk.setDisabled(False)

        self.t2_ipt.setDisabled(False)
        self.dc_widget.t2_chk.setDisabled(False)

        self.gr2_ipt.setDisabled(False)
        self.dc_widget.g2_chk.setDisabled(False)

        self.alb2_ipt.setDisabled(False)
        self.dc_widget.alb2_chk.setDisabled(False)

        self.dc_widget.x1_chk.setDisabled(False)
        self.dc_widget.x2_chk.setDisabled(False)
        self.dc_widget.l2_chk.setDisabled(False)

        self.ipb_chk.setDisabled(False)

    def check_jdphs(self):
        jdphs = str(self.jdphs_combobox.currentText())
        self.pshift_ipt.setDisabled(False)
        self.dc_widget.pshift_chk.setDisabled(False)
        self.dc_widget.jd0_chk.setDisabled(False)
        if jdphs == "Time":
            self.pshift_ipt.setDisabled(True)
            self.pshift_ipt.setValue(0.0)
            self.dc_widget.pshift_chk.setDisabled(True)
            self.dc_widget.pshift_chk.setChecked(False)
        elif jdphs == "Phase":
            self.dc_widget.jd0_chk.setDisabled(True)
            self.dc_widget.jd0_chk.setChecked(False)

    def update_input_pairs(self):
        sender = self.sender()

        if str(self.mode_combobox.currentText()) == "Mode 1":

            input_pairs = {
                self.pot1_ipt: self.pot2_ipt,
                self.t1_ipt: self.t2_ipt,
                self.gr1_ipt: self.gr2_ipt,
                self.alb1_ipt: self.alb2_ipt,
            }

            input_pairs[sender].setValue(sender.value())

        if str(self.mode_combobox.currentText()) == "Mode 3":
            if str(sender.objectName()) == "pot1_ipt":
                self.pot2_ipt.setValue(self.pot1_ipt.value())

    def check_ld(self):
        self.dc_widget.x1_chk.setDisabled(False)
        self.dc_widget.x2_chk.setDisabled(False)

        if self.ld1_chk.isChecked() is not True:
            self.dc_widget.x1_chk.setDisabled(True)
            self.dc_widget.x1_chk.setChecked(False)

        if self.ld2_chk.isChecked() is not True:
            self.dc_widget.x2_chk.setDisabled(True)
            self.dc_widget.x2_chk.setChecked(False)

        if str(self.mode_combobox.currentText()) == "Mode 1":
            if self.ld1_chk.isChecked() and self.ld2_chk.isChecked():
                for curve in self.loadobservations_widget.light_curves:
                    curve.x2 = curve.x1
                    curve.y2 = curve.y1
                self.dc_widget.x2_chk.setDisabled(True)

    def check_ipb(self):
        if self.ipb_chk.isChecked():
            self.dc_widget.l2_chk.setDisabled(False)
        else:
            self.dc_widget.l2_chk.setDisabled(True)
            self.dc_widget.l2_chk.setChecked(False)

    def update_potentials(self):
        inner_potential = float("nan")
        if self.q_ipt.value() != 0.0 and self.e_ipt.value() < 1.0:
            w = self.perr0_ipt.value()
            e = self.e_ipt.value()
            q = self.q_ipt.value()
            phase_shift = self.pshift_ipt.value()
            phase_of_periastron = methods.compute_conjunction_phases(w, e, phase_shift)[4]

            inner_potential = methods.compute_roche_potentials(w, e, q, phase_of_periastron, phase_shift)[0]

        if str(self.mode_combobox.currentText()) == "Mode 4":
            self.pot1_ipt.setValue(inner_potential)
        if str(self.mode_combobox.currentText()) == "Mode 5":
            self.pot2_ipt.setValue(inner_potential)
        if str(self.mode_combobox.currentText()) == "Mode 6":
            self.pot1_ipt.setValue(inner_potential)
            self.pot2_ipt.setValue(inner_potential)
