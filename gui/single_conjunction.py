# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'single_conjunction.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ConjunctionWidget(object):
    def setupUi(self, ConjunctionWidget):
        ConjunctionWidget.setObjectName("ConjunctionWidget")
        ConjunctionWidget.resize(500, 220)
        ConjunctionWidget.setMinimumSize(QtCore.QSize(500, 220))
        ConjunctionWidget.setMaximumSize(QtCore.QSize(500, 220))
        self.gridLayout = QtWidgets.QGridLayout(ConjunctionWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.groupBox = QtWidgets.QGroupBox(ConjunctionWidget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.primary_label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.primary_label.sizePolicy().hasHeightForWidth())
        self.primary_label.setSizePolicy(sizePolicy)
        self.primary_label.setObjectName("primary_label")
        self.gridLayout_2.addWidget(self.primary_label, 0, 0, 1, 1)
        self.f_quadrature_label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.f_quadrature_label.sizePolicy().hasHeightForWidth())
        self.f_quadrature_label.setSizePolicy(sizePolicy)
        self.f_quadrature_label.setObjectName("f_quadrature_label")
        self.gridLayout_2.addWidget(self.f_quadrature_label, 0, 1, 1, 1)
        self.p_ast_label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.p_ast_label.sizePolicy().hasHeightForWidth())
        self.p_ast_label.setSizePolicy(sizePolicy)
        self.p_ast_label.setObjectName("p_ast_label")
        self.gridLayout_2.addWidget(self.p_ast_label, 0, 2, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.prieclipse_opt = QtWidgets.QLineEdit(self.groupBox)
        self.prieclipse_opt.setReadOnly(True)
        self.prieclipse_opt.setObjectName("prieclipse_opt")
        self.horizontalLayout.addWidget(self.prieclipse_opt)
        self.firstquad_opt = QtWidgets.QLineEdit(self.groupBox)
        self.firstquad_opt.setReadOnly(True)
        self.firstquad_opt.setObjectName("firstquad_opt")
        self.horizontalLayout.addWidget(self.firstquad_opt)
        self.periastron_opt = QtWidgets.QLineEdit(self.groupBox)
        self.periastron_opt.setReadOnly(True)
        self.periastron_opt.setObjectName("periastron_opt")
        self.horizontalLayout.addWidget(self.periastron_opt)
        self.gridLayout_2.addLayout(self.horizontalLayout, 1, 0, 1, 3)
        self.s_eclipse_label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s_eclipse_label.sizePolicy().hasHeightForWidth())
        self.s_eclipse_label.setSizePolicy(sizePolicy)
        self.s_eclipse_label.setObjectName("s_eclipse_label")
        self.gridLayout_2.addWidget(self.s_eclipse_label, 2, 0, 1, 1)
        self.secondary_label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.secondary_label.sizePolicy().hasHeightForWidth())
        self.secondary_label.setSizePolicy(sizePolicy)
        self.secondary_label.setObjectName("secondary_label")
        self.gridLayout_2.addWidget(self.secondary_label, 2, 1, 1, 1)
        self.a_past_label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.a_past_label.sizePolicy().hasHeightForWidth())
        self.a_past_label.setSizePolicy(sizePolicy)
        self.a_past_label.setObjectName("a_past_label")
        self.gridLayout_2.addWidget(self.a_past_label, 2, 2, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.sececlipse_opt = QtWidgets.QLineEdit(self.groupBox)
        self.sececlipse_opt.setReadOnly(True)
        self.sececlipse_opt.setObjectName("sececlipse_opt")
        self.horizontalLayout_2.addWidget(self.sececlipse_opt)
        self.secondquad_opt = QtWidgets.QLineEdit(self.groupBox)
        self.secondquad_opt.setReadOnly(True)
        self.secondquad_opt.setObjectName("secondquad_opt")
        self.horizontalLayout_2.addWidget(self.secondquad_opt)
        self.apastron_opt = QtWidgets.QLineEdit(self.groupBox)
        self.apastron_opt.setReadOnly(True)
        self.apastron_opt.setObjectName("apastron_opt")
        self.horizontalLayout_2.addWidget(self.apastron_opt)
        self.gridLayout_2.addLayout(self.horizontalLayout_2, 3, 0, 1, 3)
        self.conj_btn = QtWidgets.QPushButton(self.groupBox)
        self.conj_btn.setObjectName("conj_btn")
        self.gridLayout_2.addWidget(self.conj_btn, 4, 0, 1, 3)
        self.gridLayout.addWidget(self.groupBox, 0, 0, 1, 1)

        self.retranslateUi(ConjunctionWidget)
        QtCore.QMetaObject.connectSlotsByName(ConjunctionWidget)
        ConjunctionWidget.setTabOrder(self.prieclipse_opt, self.firstquad_opt)
        ConjunctionWidget.setTabOrder(self.firstquad_opt, self.periastron_opt)
        ConjunctionWidget.setTabOrder(self.periastron_opt, self.sececlipse_opt)
        ConjunctionWidget.setTabOrder(self.sececlipse_opt, self.secondquad_opt)
        ConjunctionWidget.setTabOrder(self.secondquad_opt, self.apastron_opt)
        ConjunctionWidget.setTabOrder(self.apastron_opt, self.conj_btn)

    def retranslateUi(self, ConjunctionWidget):
        _translate = QtCore.QCoreApplication.translate
        ConjunctionWidget.setWindowTitle(_translate("ConjunctionWidget", "Computed Conjunctions"))
        self.groupBox.setTitle(_translate("ConjunctionWidget", "Conjunction Phases"))
        self.primary_label.setText(_translate("ConjunctionWidget", "Primary Eclipse"))
        self.f_quadrature_label.setText(_translate("ConjunctionWidget", "First Quadrature"))
        self.p_ast_label.setText(_translate("ConjunctionWidget", "Periastron"))
        self.s_eclipse_label.setText(_translate("ConjunctionWidget", "Secondary Eclipse"))
        self.secondary_label.setText(_translate("ConjunctionWidget", "Second Quadrature"))
        self.a_past_label.setText(_translate("ConjunctionWidget", "Apastron"))
        self.conj_btn.setText(_translate("ConjunctionWidget", "Compute Conjunction Phases"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ConjunctionWidget = QtWidgets.QWidget()
    ui = Ui_ConjunctionWidget()
    ui.setupUi(ConjunctionWidget)
    ConjunctionWidget.show()
    sys.exit(app.exec_())

