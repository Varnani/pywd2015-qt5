from PyQt5 import QtWidgets, QtGui
from gui import single_conjunction
from src import constants
from src.helpers import methods


class Widget(QtWidgets.QWidget, single_conjunction.Ui_ConjunctionWidget):
    def __init__(self, parent):
        super(Widget, self).__init__()
        self.setupUi(self)

        self.setWindowIcon(QtGui.QIcon(constants.MAIN_ICON_PATH))

        self.main_window = parent

        self.conj_btn.clicked.connect(self.compute_conjunctions)

    def compute_conjunctions(self):
        w = self.main_window.perr0_ipt.value()
        e = self.main_window.e_ipt.value()
        phase_shift = self.main_window.pshift_ipt.value()

        conjunctions = methods.compute_conjunction_phases(w, e, phase_shift)

        primary_eclipse, first_quadrature, secondary_eclipse, second_quadrature, periastron, apastron = conjunctions

        self.prieclipse_opt.setText("{:0.8f}".format(primary_eclipse))
        self.firstquad_opt.setText("{:0.8f}".format(first_quadrature))
        self.sececlipse_opt.setText("{:0.8f}".format(secondary_eclipse))
        self.secondquad_opt.setText("{:0.8f}".format(second_quadrature))
        self.periastron_opt.setText("{:0.8f}".format(periastron))
        self.apastron_opt.setText("{:0.8f}".format(apastron))
