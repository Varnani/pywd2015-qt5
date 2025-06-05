from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import numpy
import pickle
from src import constants


delimiter = constants.EXPORT_DELIMITER


class MatplotlibWidget:
    def __init__(self, widget, rows, cols, h_ratios=None, w_ratios=None,
                 poppable=True, exportable=True, grid_allowed=True):

        if h_ratios is None:
            h_ratios = numpy.ones(rows)

        if w_ratios is None:
            w_ratios = numpy.ones(cols)

        # properties
        self.rows = rows
        self.cols = cols
        self.widget = widget

        self.axes = []
        self.titles = []
        self.labels = []

        self.figure = Figure()
        self.gridspec = GridSpec(rows, cols, figure=self.figure,
                                 height_ratios=list(h_ratios), width_ratios=list(w_ratios))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, widget)

        # custom ui elements
        self.pop_button = QtWidgets.QPushButton(widget)
        self.pop_button.setText("Pop Plot Window")

        self.export_button = QtWidgets.QPushButton(widget)
        self.export_button.setText("Export Currently Plotted Data")

        self.grid_checkbox = QtWidgets.QCheckBox(widget)
        self.grid_checkbox.setText("Enable Grid")
        self.grid_checkbox.setChecked(False)

        self.popped_windows = []

        if not poppable:
            self.pop_button.hide()

        if not exportable:
            self.export_button.hide()

        if not grid_allowed:
            self.grid_checkbox.hide()

        # set layouts
        self.layout_horizontal = QtWidgets.QHBoxLayout()
        self.layout_horizontal.addWidget(self.pop_button)
        self.layout_horizontal.addWidget(self.export_button)
        self.layout_horizontal.addWidget(self.grid_checkbox)

        self.layout_vertical = QtWidgets.QVBoxLayout(widget)
        self.layout_vertical.addWidget(self.canvas)
        self.layout_vertical.addWidget(self.toolbar)
        self.layout_vertical.addLayout(self.layout_horizontal)

        # connect signals
        self.connect_signals()

    def connect_signals(self):
        self.pop_button.clicked.connect(self.pop)
        self.export_button.clicked.connect(self.export)
        self.grid_checkbox.toggled.connect(self.toggle_grid)

    def clean_popped_windows(self):
        for window in self.popped_windows:
            if not window.isVisible():
                a = self.popped_windows.pop(self.popped_windows.index(window))
                a.setParent(None)
                a.deleteLater()
                del a

    def close_popped_windows(self):
        self.clean_popped_windows()
        for window in self.popped_windows:
            window.close()
        self.popped_windows = []

    def pop(self):
        self.clean_popped_windows()
        window = PopWidget(self)
        self.popped_windows.append(window)
        window.show()

    def export(self):
        output = "# Data exported from a matplotlib figure\n\n"

        for idx, ax in enumerate(self.axes):
            output = output + "# Axis {0}\n".format(str(idx + 1))
            for line in ax.get_lines():
                xdata = line.get_xdata()
                ydata = line.get_ydata()

                try:
                    float(xdata[0])
                    float(ydata[0])

                except:
                    # if we cant cast an array element into a float, most likely we dont need it
                    continue

                for x, y in zip(xdata, ydata):  # "{:0.5f}"
                    output = output + "{0:0.5f}{2}{1:0.10f}\n".format(x, y, delimiter)

            output = output + "\n"

        dialog = QtWidgets.QFileDialog()
        dialog.setDefaultSuffix(".dat")
        dialog.setAcceptMode(1)
        return_code = dialog.exec_()
        if return_code != 0:
            path = str((dialog.selectedFiles())[0])
            with open(path, "w") as destination:
                destination.write(output)

            msg = QtWidgets.QMessageBox()
            msg.setText("File saved:")
            msg.setInformativeText(path)
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.exec_()

    def toggle_grid(self):

        for ax in self.axes:
            ax.grid(visible=self.grid_checkbox.isChecked())

        self.redraw()

    def create_axis(self, row_idx, col_idx, name="", labels=("", ""), *subplot_args, **subplot_kwargs):

        # setup axis
        if isinstance(row_idx, tuple):
            row_idx = slice(row_idx[0], row_idx[1])
        if isinstance(col_idx, tuple):
            col_idx = slice(col_idx[0], col_idx[1])
        self.axes.append(self.figure.add_subplot(self.gridspec[row_idx, col_idx], *subplot_args, **subplot_kwargs))

        # setup title
        self.titles.append(name)
        self.axes[-1].set_title(name)

        # setup labels
        self.labels.append(labels)
        self.axes[-1].set_xlabel(labels[0])
        self.axes[-1].set_ylabel(labels[1])

        self.redraw()

    def set_axis_title(self, title, index=0):
        self.titles[index] = title
        self.axes[index].set_title(title)
        self.redraw()

    def set_labels(self, x_label, y_label, index=0):
        self.labels[index] = (x_label, y_label)
        self.axes[index].set_xlabel(x_label)
        self.axes[index].set_ylabel(y_label)
        self.redraw()

    def plot(self, x, y, index=0, clear=True, update_toolbar=True, *plot_args, **plot_kwargs):

        # clear plot if requested
        if clear:
            self.clear_axis(index)

        # plot the data
        self.axes[index].plot(x, y, *plot_args, **plot_kwargs)

        # reorient axis if requested
        if update_toolbar:
            self.toolbar.update()

        # redraw title
        self.axes[index].set_title(self.titles[index])

        # redraw labels
        self.axes[index].set_xlabel(self.labels[index][0])
        self.axes[index].set_ylabel(self.labels[index][1])
        self.axes[index].ticklabel_format(useOffset=False)

        # finally, draw the canvas
        self.canvas.draw()

    def get_plot_data(self, index=0):
        x = self.axes[index].lines[0].get_xdata()
        y = self.axes[index].lines[0].get_ydata()
        return x, y

    def clear_all(self):
        for axis in self.axes:
            axis.clear()
        self.redraw()

    def clear_axis(self, index):
        self.axes[index].clear()

    def redraw(self):
        self.canvas.draw()

    def __getitem__(self, item):
        return self.axes[item]


class PopWidget(QtWidgets.QWidget):
    def __init__(self, parent):
        super(PopWidget, self).__init__()
        self.parent_widget = parent

        # create a simple container ui
        self.setObjectName("PopPlotDialog")
        self.setWindowTitle("PyWD2015")
        self.resize(900, 500)
        self.setMinimumSize(300, 300)
        h_layout = QtWidgets.QHBoxLayout(self)
        self.plot_widget = QtWidgets.QWidget(self)

        # create a new canvas from the old figure
        # we'll just dump and load the parent figure to create a new figure instance
        # otherwise parent widget's canvas will be broken
        self.canvas = FigureCanvas(figure=pickle.loads(pickle.dumps(self.parent_widget.figure)))
        self.toolbar = NavigationToolbar(self.canvas, self.plot_widget)

        # fill layouts with created widgets
        h_layout.addWidget(self.plot_widget)
        v_layout = QtWidgets.QVBoxLayout(self.plot_widget)
        v_layout.addWidget(self.canvas)
        v_layout.addWidget(self.toolbar)

        QtCore.QMetaObject.connectSlotsByName(self)

    def closeEvent(self, *args, **kwargs):
        for ax in self.canvas.figure.get_axes():
            ax.clear()
        del self.canvas
        del self.toolbar
        self.parent_widget.clean_popped_windows()
