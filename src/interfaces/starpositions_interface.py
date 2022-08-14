from PyQt5 import QtWidgets, QtGui, QtCore
from gui import starpositions_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src import constants
from src.helpers.wd_utils import wd_io
from src.helpers import methods
import time
from matplotlib import pyplot
import io
import os


class Widget(QtWidgets.QWidget, starpositions_widget.Ui_StarPositionWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.chart = MatplotlibWidget(self.plot_widget, 1, 1, exportable=False)
        self.chart.create_axis(0, 0, labels=("X", "Y"))
        self.chart.axes[0].axis("equal")

        self.start_btn.setIcon(QtGui.QIcon(constants.PLAY_ICON_PATH))
        self.pause_btn.setIcon(QtGui.QIcon(constants.PAUSE_ICON_PATH))
        self.next_btn.setIcon(QtGui.QIcon(constants.NEXT_ICON_PATH))
        self.prev_btn.setIcon(QtGui.QIcon(constants.PREV_ICON_PATH))

        # self.viewport_pixlabel.setScaledContents(True)

        # animation variables
        self.frames_per_second = 25.0  # 25 frames per second
        self.frame_time = 1.0 / self.frames_per_second  # 0.04ms per frame (25fps)
        self.playback = False
        self.Stopwatch = Stopwatch()
        self.Stopwatch.wait = self.frame_time

        self.star_position_data = None
        self.rendered_frames = None
        self.current_frame = None
        self.single_frame = None

        self.connect_signals()

    def connect_signals(self):
        self.roche_groupbox.toggled.connect(self.handle_toggle_roche)
        self.positions_groupbox.toggled.connect(self.handle_toggle_position)
        self.plot_btn.clicked.connect(self.plot)

        # animator
        self.horizontalSlider.sliderMoved.connect(self.slider_move)
        self.horizontalSlider.sliderPressed.connect(self.slider_move)
        self.start_btn.clicked.connect(self.start_playback)
        self.pause_btn.clicked.connect(self.stop_playback)
        self.next_btn.clicked.connect(self.next_frame)
        self.prev_btn.clicked.connect(self.previous_frame)
        self.Stopwatch.tick.connect(self.advance_frame)
        self.render_btn.clicked.connect(self.start_render)
        self.saveframe_btn.clicked.connect(self.save_frame)
        self.saveall_btn.clicked.connect(self.save_all_frames)

    def resizeEvent(self, resize_event):
        self.draw_current_frame()
        super(Widget, self).resizeEvent(resize_event)

    def handle_toggle_roche(self):
        self.roche_groupbox.blockSignals(True)
        self.positions_groupbox.blockSignals(True)

        if self.roche_groupbox.isChecked():
            self.positions_groupbox.setChecked(False)
        elif self.positions_groupbox.isChecked() is False:
            self.roche_groupbox.setChecked(True)

        self.roche_groupbox.blockSignals(False)
        self.positions_groupbox.blockSignals(False)

    def handle_toggle_position(self):
        self.roche_groupbox.blockSignals(True)
        self.positions_groupbox.blockSignals(True)

        if self.positions_groupbox.isChecked():
            self.roche_groupbox.setChecked(False)
        elif self.roche_groupbox.isChecked() is False:
            self.positions_groupbox.setChecked(True)

        self.roche_groupbox.blockSignals(False)
        self.positions_groupbox.blockSignals(False)

    def plot(self):
        if self.positions_groupbox.isChecked():
            self.plot_positions()
        elif self.roche_groupbox.isChecked():
            self.plot_roche()

    def plot_positions(self):
        self.chart.set_axis_title("Star Positions")

        lc_params = self.main_window.get_lc_params()

        lc_params["jdphs"] = 2
        lc_params["phstrt"] = self.phase_spinbox.value()
        lc_params["phstop"] = self.phase_spinbox.value()
        lc_params["phin"] = 0.1

        lc_io = wd_io.LCIO(lc_params,
                           wd_path=self.main_window.lc_path,
                           lc_binary_name=self.main_window.lc_binary)

        results = lc_io.fill_for_star_positions().save().run().read_star_positions()[0]

        self.chart.plot(results[0], results[1], linestyle="", marker="+", markersize=1, color="black")
        self.chart.plot([0], [0], clear=False, linestyle="", marker="+", markersize=10, color=constants.COLOR_RED)

    def plot_roche(self):
        # compute_roche_potentials(w, e, q, phase, phase_shift, plot_elements=None):

        self.chart.clear_all()

        w = self.main_window.perr0_ipt.value()
        e = self.main_window.e_ipt.value()
        q = self.main_window.q_ipt.value()
        phase = self.main_window.p0_ipt.value()
        pshift = self.main_window.pshift_ipt.value()
        pot1 = self.main_window.pot1_ipt.value()
        pot2 = self.main_window.pot2_ipt.value()

        inner_critical, outer_critical = methods.compute_roche_potentials(w, e, q, phase, pshift,
                                                                          plot_elements=[self.chart.axes[0], pot1, pot2])
        self.chart.set_labels("X", "Y")
        self.chart.set_axis_title("Roche Potentials")
        self.chart.redraw()

        self.inner_crit_otpt.setValue(inner_critical)
        self.outer_crit_otpt.setValue(outer_critical)

    def setup_slider(self):
        self.horizontalSlider.setMinimum(0)
        self.horizontalSlider.setMaximum(len(self.rendered_frames) - 1)
        self.horizontalSlider.setTickInterval(1)

    def slider_move(self):
        self.stop_playback()
        if self.rendered_frames is not None:
            self.show_frame(self.rendered_frames[self.horizontalSlider.value()])

    def start_playback(self):
        if self.rendered_frames is not None:
            self.Stopwatch.start()

    def stop_playback(self):
        if self.rendered_frames is not None:
            self.Stopwatch.stop()

    def next_frame(self):
        self.stop_playback()
        self.advance_frame()

    def advance_frame(self):
        if self.rendered_frames is not None:
            index = self.horizontalSlider.value() + 1
            if index == len(self.rendered_frames):
                index = 0

            self.horizontalSlider.setValue(index)
            self.show_frame(self.rendered_frames[index])

    def previous_frame(self):
        self.stop_playback()
        if self.rendered_frames is not None:
            index = self.horizontalSlider.value() - 1
            if index == -1:
                index = len(self.rendered_frames) - 1

            self.horizontalSlider.setValue(index)
            self.show_frame(self.rendered_frames[index])

    def clear_animator(self):
        self.viewport_pixlabel.clear()
        self.rendered_frames = None
        self.star_position_data = None

    def start_render(self):
        self.clear_animator()
        if self.single_chk.isChecked() is not True:
            #if self.render_stars() is not 1:
            if self.render_stars() != 1:
                self.setup_slider()

        else:
            self.render_single()

    def update_message_label(self, msg):
        self.message_label.setText(msg)
        self.message_label.repaint()

    def update_progress_bar(self, value):
        self.progressBar.setValue(value)
        self.progressBar.repaint()

    def render_stars(self):
        increment = None
        fmt = None
        iterations = None

        if self.main_window.jd_radiobtn.isChecked():
            mn = self.main_window.lc_jd_start_ipt.value()
            mx = self.main_window.lc_jd_end_ipt.value()
            increment = self.main_window.lc_jd_incr_ipt.value()
            fmt = "{:7.6f}"

            iterations = int((mx - mn) / increment)

        elif self.main_window.phase_radiobtn.isChecked():
            mn = self.main_window.lc_phs_start_ipt.value()
            mx = self.main_window.lc_phs_stop_ipt.value()
            increment = self.main_window.lc_phs_incr_ipt.value()
            fmt = "{:4.3f}"

            iterations = int((mx - mn) / increment)

        if iterations > 500:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(msg.Warning)
            msg.setText("Expected iteration count is larger than 500. (" + str(iterations) + ")")
            msg.setInformativeText("This might result in a very large lcout file (>100MB), "
                                   "take a long time and might crash LC altogether. "
                                   "Are you sure you want to render the animation?")
            msg.setStandardButtons(msg.Ok | msg.Cancel)
            if msg.exec_() == msg.Cancel:
                return 1

        self.update_message_label("Running LC...")

        lc_params = self.main_window.get_lc_params()

        lc_io = wd_io.LCIO(lc_params,
                           wd_path=self.main_window.lc_path,
                           lc_binary_name=self.main_window.lc_binary)

        results = lc_io.fill_for_star_positions().save().run().read_star_positions()

        self.rendered_frames = []

        self.update_message_label("Rendering plots...")

        progress_increment = 100.0 / float(len(results))
        current_progress = 0.0

        for idx, result in enumerate(results):
            qpixmap = self.render_frame(result[0], result[1], increment * idx, fmt=fmt)
            self.rendered_frames.append(qpixmap)

            current_progress = current_progress + progress_increment
            self.update_progress_bar(current_progress)

        self.show_frame(self.rendered_frames[0])

        self.update_message_label("Done.")
        self.update_progress_bar(100)

        self.single_frame = False

    def render_single(self):
        lc_params = self.main_window.get_lc_params()

        lc_params["jdphs"] = 2
        lc_params["phstrt"] = self.render_phaseSpinbox.value()
        lc_params["phstop"] = self.render_phaseSpinbox.value()
        lc_params["phin"] = 0.1

        lc_io = wd_io.LCIO(lc_params,
                           wd_path=self.main_window.lc_path,
                           lc_binary_name=self.main_window.lc_binary)

        results = lc_io.fill_for_star_positions().save().run().read_star_positions()[0]

        qpixmap = self.render_frame(results[0], results[1], self.render_phaseSpinbox.value())

        self.show_frame(qpixmap)

        self.single_frame = True

    def render_frame(self, x, y, t, fmt="{:4.3f}"):
        pyplot.cla()
        pyplot.axis("equal")
        pyplot.xlabel("X")
        pyplot.ylabel("Y")
        pyplot.xlim(self.min_spinbox.value(), self.max_spinbox.value())
        pyplot.ylim(self.min_spinbox.value(), self.max_spinbox.value())
        pyplot.plot(x, y, 'ko', markersize=0.2, label=fmt.format(t))
        pyplot.legend(loc="upper right")

        pyplot.plot([0], [0], linestyle="", marker="+", markersize=10, color="#ff3a3a")

        image = io.BytesIO()

        dpi_dict = {
            "64dpi": 64,
            "128dpi": 128,
            "196dpi": 196,
            "256dpi": 256
        }

        pyplot.savefig(image, dpi=dpi_dict[str(self.dpi_combobox.currentText())], format="png")
        image.seek(0)

        qbyte = QtCore.QByteArray(image.getvalue())
        qpixmap = QtGui.QPixmap()
        qpixmap.loadFromData(qbyte, "png")

        return qpixmap

    def show_frame(self, qpixmap):
        self.current_frame = qpixmap
        self.draw_current_frame()

    def draw_current_frame(self):
        if self.current_frame is not None:
            w = self.viewport_pixlabel.width()
            h = self.viewport_pixlabel.height()

            # 1 means keep aspect ratio, it is a Qt enum (Qt::KeepAspectRatio)
            self.viewport_pixlabel.setPixmap(QtGui.QPixmap(self.current_frame).scaled(w, h, 1))

    def save_frame(self):

        frame = None

        if self.single_frame:
            frame = self.current_frame

        elif self.rendered_frames is not None:
            frame = self.rendered_frames[self.horizontalSlider.value()]

        if frame is not None:
            dialog = QtWidgets.QFileDialog()
            dialog.setDefaultSuffix("png")
            dialog.setNameFilter("PNG File (*png)")
            dialog.setAcceptMode(1)

            return_code = dialog.exec_()
            file_path = str((dialog.selectedFiles())[0])

            if file_path != "" and return_code != 0:
                frame.save(file_path, "png", 100)
                msg = QtWidgets.QMessageBox()
                msg.setText("Frame is saved into " + file_path)
                msg.exec_()

    def save_all_frames(self):
        if self.rendered_frames is not None:
            dialog = QtWidgets.QFileDialog()
            dialog.setFileMode(2)

            return_code = dialog.exec_()
            file_path = str((dialog.selectedFiles())[0])

            if file_path != "" and return_code != 0:
                for idx, qpixmap in enumerate(self.rendered_frames):
                    qpixmap.save(os.path.join(file_path, "{:0>4d}".format(idx) + ".png"), "png", 100)

                msg = QtWidgets.QMessageBox()
                msg.setText("Frames are saved into " + file_path)
                msg.exec_()


class Stopwatch(QtCore.QThread):
    tick = QtCore.pyqtSignal(name="tick")

    def __init__(self):
        QtCore.QThread.__init__(self)
        self.wait = None
        self.running = True

    def run(self):
        self.running = True
        while self.running:
            time.sleep(self.wait)
            self.tick.emit()

    def stop(self):
        self.running = False
