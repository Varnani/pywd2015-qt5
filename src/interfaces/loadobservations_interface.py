from PyQt5 import QtWidgets, QtGui
from gui import loadobservations_widget
import curveproperties_interface
from matplotlib import pyplot
import os
from src import constants


class Widget(QtWidgets.QWidget, loadobservations_widget.Ui_ObservationWidget):
    def __init__(self, main_window):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = main_window

        self.velocity_curves = [None, None]
        self.light_curves = []

        self.curve_treewidget.header().setSectionResizeMode(3)

        self.connect_signals()

    def connect_signals(self):
        self.add_btn.clicked.connect(self.add_observation)
        self.edit_btn.clicked.connect(self.edit_observation)
        self.remove_btn.clicked.connect(self.remove_observation)
        self.plot_btn.clicked.connect(self.plot_observation)

    def selected_item(self):
        selected_item = self.curve_treewidget.selectedItems()
        if len(selected_item) > 0:
            return selected_item[0]
        else:
            return None

    def add_observation(self):
        menu = QtWidgets.QMenu(self)

        # vc menu
        addvc = QtWidgets.QMenu("Velocity Curve")
        menu.addMenu(addvc)
        pri = addvc.addAction("Primary")
        pri.setObjectName("primary")
        sec = addvc.addAction("Secondary")
        sec.setObjectName("secondary")

        if self.velocity_curves[0] is not None:
            pri.setDisabled(True)

        if self.velocity_curves[1] is not None:
            sec.setDisabled(True)

        # lc action
        addlc = menu.addAction("Light Curve")
        addlc.setObjectName("lightcurve")

        selection = menu.exec_(QtGui.QCursor.pos())

        if selection is not None:
            curve_type = ""
            if selection.objectName() == "primary" or selection.objectName() == "secondary":
                curve_type = curveproperties_interface.TYPE_VELOCITY
            elif selection.objectName() == "lightcurve":
                curve_type = curveproperties_interface.TYPE_LIGHT

            curve = curveproperties_interface.CurveProperty(curve_type, self.main_window)
            curve.show_dialog()

            if curve.unpopulated is not True:
                if selection.objectName() == "primary":
                    self.velocity_curves[0] = curve
                elif selection.objectName() == "secondary":
                    self.velocity_curves[1] = curve
                elif selection.objectName() == "lightcurve":
                    self.light_curves.append(curve)

                self.update_curve_list()

    def edit_observation(self):
        item = self.selected_item()
        if item is not None:
            if item.text(1) == "Velocity Curve (Star 1)":
                self.velocity_curves[0].edit()

            elif item.text(1) == "Velocity Curve (Star 2)":
                self.velocity_curves[1].edit()

            elif item.text(1) == "Light Curve":
                index = self.curve_treewidget.invisibleRootItem().indexOfChild(item)

                if self.velocity_curves[0] is not None:
                    index = index - 1

                if self.velocity_curves[1] is not None:
                    index = index - 1

                self.light_curves[index].edit()

            self.update_curve_list()

    def remove_observation(self):
        item = self.selected_item()
        if item is not None:
            if item.text(1) == "Velocity Curve (Star 1)":
                self.velocity_curves[0] = None

            elif item.text(1) == "Velocity Curve (Star 2)":
                self.velocity_curves[1] = None

            elif item.text(1) == "Light Curve":
                index = self.curve_treewidget.invisibleRootItem().indexOfChild(item)

                if self.velocity_curves[0] is not None:
                    index = index - 1

                if self.velocity_curves[1] is not None:
                    index = index - 1

                self.light_curves.pop(index)

            self.update_curve_list()

    def plot_observation(self):
        item = self.selected_item()
        curve = None
        if item is not None:
            if item.text(1) == "Velocity Curve (Star 1)":
                curve = self.velocity_curves[0]

            elif item.text(1) == "Velocity Curve (Star 2)":
                curve = self.velocity_curves[1]

            elif item.text(1) == "Light Curve":
                index = self.curve_treewidget.invisibleRootItem().indexOfChild(item)

                if self.velocity_curves[0] is not None:
                    index = index - 1

                if self.velocity_curves[1] is not None:
                    index = index - 1

                curve = self.light_curves[index]

            x, y, z = curve.get_data()

            pyplot.plot(x, y, linestyle="", marker="o", markersize=2, color=constants.COLOR_BLUE)
            pyplot.get_current_fig_manager().set_window_title("Matplotlib - " + item.text(0))
            if str(self.main_window.maglite_combobox.currentText()) == "Magnitude" and str(item.text(1)) == "Light Curve":
                pyplot.gca().invert_yaxis()
            pyplot.show()

    def update_curve_list(self):
        self.curve_treewidget.clear()

        if self.velocity_curves[0] is not None:
            itm = QtWidgets.QTreeWidgetItem(self.curve_treewidget)
            itm.setText(0, os.path.basename(self.velocity_curves[0].filepath))
            itm.setText(1, "Velocity Curve (Star 1)")
            itm.setText(2, " - ")
            # itm.setText(2, constants.ID_BANDPASS_DICT[str(self.velocity_curves[0].band_id)])

        if self.velocity_curves[1] is not None:
            itm = QtWidgets.QTreeWidgetItem(self.curve_treewidget)
            itm.setText(0, os.path.basename(self.velocity_curves[1].filepath))
            itm.setText(1, "Velocity Curve (Star 2)")
            itm.setText(2, " - ")
            # itm.setText(2, constants.ID_BANDPASS_DICT[str(self.velocity_curves[1].band_id)])

        for light_curve in self.light_curves:
            itm = QtWidgets.QTreeWidgetItem(self.curve_treewidget)
            itm.setText(0, os.path.basename(light_curve.filepath))
            itm.setText(1, "Light Curve")
            itm.setText(2, constants.ID_BANDPASS_DICT[str(light_curve.band_id)])

        self.main_window.curve_list_changed()

    def write_into_parser(self, parser):
        parser.add_section(constants.CONFIG_SECTION_CURVE_COUNTS)
        parser.set(constants.CONFIG_SECTION_CURVE_COUNTS, "light_curves", str(len(self.light_curves)))

        i = 1
        for light_curve in self.light_curves:
            section = constants.CONFIG_SECTION_LIGHT_CURVE_STUB + str(i)
            parser.add_section(section)
            parser = light_curve.populate_parser(parser, section)
            i = i + 1

        if self.velocity_curves[0] is not None:
            parser.set(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_1", str(True))

            section = constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + "1"
            parser.add_section(section)
            parser = self.velocity_curves[0].populate_parser(parser, section)

        else:
            parser.set(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_1", str(False))

        if self.velocity_curves[1] is not None:
            parser.set(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_2", str(True))

            section = constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + "2"
            parser.add_section(section)
            parser = self.velocity_curves[1].populate_parser(parser, section)

        else:
            parser.set(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_2", str(False))

        return parser

    def read_from_parser(self, parser):

        self.velocity_curves = [None, None]
        self.light_curves = []

        if parser.getboolean(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_1") is True:
            curve = curveproperties_interface.CurveProperty(curveproperties_interface.TYPE_VELOCITY, self.main_window)
            curve.populate_from_parser(parser, constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + "1")
            self.velocity_curves[0] = curve

        if parser.getboolean(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_2") is True:
            curve = curveproperties_interface.CurveProperty(curveproperties_interface.TYPE_VELOCITY, self.main_window)
            curve.populate_from_parser(parser, constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + "2")
            self.velocity_curves[1] = curve

        lc_count = parser.getint(constants.CONFIG_SECTION_CURVE_COUNTS, "light_curves")
        i = 1

        while i <= lc_count:
            curve = curveproperties_interface.CurveProperty(curveproperties_interface.TYPE_LIGHT, self.main_window)
            curve.populate_from_parser(parser, constants.CONFIG_SECTION_LIGHT_CURVE_STUB + str(i))
            self.light_curves.append(curve)
            i = i + 1

        self.update_curve_list()

    def get_all_curves(self):
        curves = []

        if self.velocity_curves[0] is not None:
            curves.append(self.velocity_curves[0])

        if self.velocity_curves[1] is not None:
            curves.append(self.velocity_curves[1])

        for light_curve in self.light_curves:
            curves.append(light_curve)

        return curves
