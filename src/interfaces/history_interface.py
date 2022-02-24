from PyQt5 import QtWidgets, QtGui
from gui import history_widget
from src.helpers.matplotlib_embedder import MatplotlibWidget
from src import constants

from matplotlib.ticker import MaxNLocator
from src.helpers import methods
from src.helpers import messenger


delimiter = constants.EXPORT_DELIMITER


class Widget(QtWidgets.QWidget, history_widget.Ui_HistoryWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.ptm_font = QtGui.QFont(parent.monoFont)
        self.ptm_font.setPointSize(9)
        self.history_treewidget.setFont(self.ptm_font)
        self.history_treewidget.header().setSectionResizeMode(3)

        header_item = QtWidgets.QTreeWidgetItem()
        header_item.setText(0, "Iteration #")
        header_item.setText(1, "Parameters")
        self.history_treewidget.setHeaderItem(header_item)

        self.chart = MatplotlibWidget(self.plot_widget, 2, 1, exportable=False)
        self.chart.create_axis(0, 0)
        self.chart.create_axis(1, 0)

        self.solutions = []
        self.stats = []

        self.connect_signals()

    def connect_signals(self):
        self.plot_btn.clicked.connect(self.plot_solutions)
        self.clear_btn.clicked.connect(self.clear)
        self.export_btn.clicked.connect(self.export_data)

    def export_data(self):
        if len(self.solutions) is not None:
            filepath = methods.save_file(self, suffix="txt")
            if filepath is not None:
                output = "# Itr" + delimiter
                for result in self.solutions[0]:
                    output = output + constants.KEEPS_ID_NAME_DICT[result[0]] + delimiter
                output = output + "MRfIV" + delimiter + "MRP" + delimiter + "\n"

                for idx, solution in enumerate(self.solutions):
                    output = output + str(idx + 1) + delimiter
                    for row in solution:
                        output = output + str(row[4]) + delimiter

                    output = output + delimiter + str("{:0.16f}".format(self.stats[idx][0][0])) + delimiter + \
                             str("{:0.16f}".format(self.stats[idx][1][0])) + "\n"

                with open(filepath, "w") as destination:
                    destination.write(output)

                msg = messenger.Messenger("info", "File saved:")
                msg.set_info(filepath)
                msg.show()

    def clear(self):
        self.chart.clear_all()
        self.solutions = []
        self.stats = []
        self.history_treewidget.clear()

    def selected_item(self):
        selecteditem = self.history_treewidget.selectedItems()
        if len(selecteditem) > 0:
            return selecteditem[0]
        else:
            return None

    def receive_solution(self, new_solution, new_stat):
        if len(self.solutions) == 0:
            self.solutions.append(new_solution)
            self.stats.append(new_stat)

        elif len(self.solutions[0]) == len(new_solution):
            same_solution = True
            old_solution = self.solutions[-1]
            for old, new in zip(old_solution, new_solution):
                if old[0] != new[0]:
                    same_solution = False
                    break

            if not same_solution:
                self.clear()
                self.solutions.append(new_solution)
                self.stats.append(new_stat)

            elif same_solution:
                self.solutions.append(new_solution)
                self.stats.append(new_stat)

        else:
            self.clear()
            self.solutions.append(new_solution)
            self.stats.append(new_stat)

        self.update_solutions()

    def update_solutions(self):
        if len(self.solutions) == 1:
            header_item = QtWidgets.QTreeWidgetItem()
            header_item.setText(0, "Itr #")
            for index, row in enumerate(self.solutions[0]):
                par_id = row[0]
                c_id = row[1]

                text = constants.KEEPS_ID_NAME_DICT[par_id]
                if c_id != 0.0:
                    text = text + " (Curve #{0})".format(str(int(c_id)))
                header_item.setText(index + 1, text)
                #header_item.setText(index + 2, "Mean Residual for Input Values")
            self.history_treewidget.setHeaderItem(header_item)

        solution = self.solutions[-1]
        item = QtWidgets.QTreeWidgetItem(self.history_treewidget)
        item.setText(0, str(len(self.solutions)))
        for index, row in enumerate(solution):
            value = row[4]
            if row[0] in (19.0, 20.0):
                value = value * 10000.0
            item.setText(index + 1, "{:0.5f}".format(value))
            #item.setText(index+2, "{:0.16f}".format(self.stats[-1][0][0]))

        if self.auto_chk.isChecked():
            self.plot_solutions()

    def plot_solutions(self):
        selected_item = self.selected_item()
        if selected_item is not None and len(self.solutions) > 0:
            self.chart.clear_all()
            index = self.history_treewidget.currentColumn() - 1
            x = list(range(1, len(self.solutions) + 1))
            y = [solution[index][4] for solution in self.solutions]
            y_w = [solution[index][5] for solution in self.solutions]

            y_s = [stat[0][0] for stat in self.stats]
            y_s2 = [stat[1][0] for stat in self.stats]

            if self.solutions[0][index][0] in (19.0, 20.0):
                _y = [yy * 10000.0 for yy in y]
                _y_w = [yyw * 10000.0 for yyw in y_w]

                y = _y
                y_w = _y_w

            self.chart.axes[0].errorbar(x, y, yerr=y_w, linestyle="-", marker="o", markersize=constants.MARKER_SIZE + 2,
                                        color=constants.COLOR_BLUE, ecolor=constants.COLOR_RED)
            self.chart.axes[0].set_ylabel(constants.KEEPS_ID_NAME_DICT[self.solutions[0][index][0]])
            self.chart.axes[0].ticklabel_format(useOffset=False)
            self.chart.axes[0].get_yaxis().set_major_locator(MaxNLocator(integer=True))
            self.chart.axes[1].plot(x, y_s, marker="o", markersize=constants.MARKER_SIZE + 2,
                                        color=constants.COLOR_BLUE)
            self.chart.axes[1].plot(x, y_s2, marker="o", markersize=constants.MARKER_SIZE + 2,
                                        color=constants.COLOR_RED)
            self.chart.axes[1].set_xlabel("Iteration Number")
            self.chart.axes[1].set_ylabel("Mean Residual\n (Input/Predicted)")
            self.chart.redraw()