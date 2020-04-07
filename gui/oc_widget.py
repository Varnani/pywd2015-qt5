# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'oc_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_OCWidget(object):
    def setupUi(self, OCWidget):
        OCWidget.setObjectName("OCWidget")
        OCWidget.resize(850, 500)
        OCWidget.setMinimumSize(QtCore.QSize(850, 500))
        self.horizontalLayout = QtWidgets.QHBoxLayout(OCWidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.plot_widget = QtWidgets.QWidget(OCWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plot_widget.sizePolicy().hasHeightForWidth())
        self.plot_widget.setSizePolicy(sizePolicy)
        self.plot_widget.setMinimumSize(QtCore.QSize(400, 300))
        self.plot_widget.setObjectName("plot_widget")
        self.gridLayout.addWidget(self.plot_widget, 0, 3, 6, 1)
        self.line_2 = QtWidgets.QFrame(OCWidget)
        self.line_2.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout.addWidget(self.line_2, 0, 2, 6, 1)
        self.groupBox = QtWidgets.QGroupBox(OCWidget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 1, 0, 1, 1)
        self.dt_otpt = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.dt_otpt.setReadOnly(True)
        self.dt_otpt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.dt_otpt.setDecimals(7)
        self.dt_otpt.setMinimum(-9999999.0)
        self.dt_otpt.setMaximum(9999999.0)
        self.dt_otpt.setObjectName("dt_otpt")
        self.gridLayout_2.addWidget(self.dt_otpt, 0, 1, 1, 1)
        self.dp_otpt = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.dp_otpt.setReadOnly(True)
        self.dp_otpt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.dp_otpt.setDecimals(10)
        self.dp_otpt.setMinimum(-9999999.0)
        self.dp_otpt.setMaximum(9999999.0)
        self.dp_otpt.setObjectName("dp_otpt")
        self.gridLayout_2.addWidget(self.dp_otpt, 1, 1, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.calculate_btn = QtWidgets.QPushButton(self.groupBox)
        self.calculate_btn.setObjectName("calculate_btn")
        self.horizontalLayout_2.addWidget(self.calculate_btn)
        self.update_btn = QtWidgets.QPushButton(self.groupBox)
        self.update_btn.setObjectName("update_btn")
        self.horizontalLayout_2.addWidget(self.update_btn)
        self.gridLayout_2.addLayout(self.horizontalLayout_2, 4, 0, 1, 2)
        self.gridLayout.addWidget(self.groupBox, 5, 0, 1, 2)
        self.dpdt_chk = QtWidgets.QCheckBox(OCWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dpdt_chk.sizePolicy().hasHeightForWidth())
        self.dpdt_chk.setSizePolicy(sizePolicy)
        self.dpdt_chk.setObjectName("dpdt_chk")
        self.gridLayout.addWidget(self.dpdt_chk, 2, 1, 1, 1)
        self.line = QtWidgets.QFrame(OCWidget)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 3, 0, 1, 2)
        self.linear_chk = QtWidgets.QCheckBox(OCWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.linear_chk.sizePolicy().hasHeightForWidth())
        self.linear_chk.setSizePolicy(sizePolicy)
        self.linear_chk.setObjectName("linear_chk")
        self.gridLayout.addWidget(self.linear_chk, 2, 0, 1, 1)
        self.data_treewidget = QtWidgets.QTreeWidget(OCWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.data_treewidget.sizePolicy().hasHeightForWidth())
        self.data_treewidget.setSizePolicy(sizePolicy)
        self.data_treewidget.setMinimumSize(QtCore.QSize(350, 0))
        self.data_treewidget.setObjectName("data_treewidget")
        self.gridLayout.addWidget(self.data_treewidget, 4, 0, 1, 2)
        self.line_3 = QtWidgets.QFrame(OCWidget)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.gridLayout.addWidget(self.line_3, 1, 0, 1, 2)
        self.compute_btn = QtWidgets.QPushButton(OCWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.compute_btn.sizePolicy().hasHeightForWidth())
        self.compute_btn.setSizePolicy(sizePolicy)
        self.compute_btn.setObjectName("compute_btn")
        self.gridLayout.addWidget(self.compute_btn, 0, 0, 1, 1)
        self.export_btn = QtWidgets.QPushButton(OCWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.export_btn.sizePolicy().hasHeightForWidth())
        self.export_btn.setSizePolicy(sizePolicy)
        self.export_btn.setObjectName("export_btn")
        self.gridLayout.addWidget(self.export_btn, 0, 1, 1, 1)
        self.horizontalLayout.addLayout(self.gridLayout)

        self.retranslateUi(OCWidget)
        QtCore.QMetaObject.connectSlotsByName(OCWidget)

    def retranslateUi(self, OCWidget):
        _translate = QtCore.QCoreApplication.translate
        OCWidget.setWindowTitle(_translate("OCWidget", "Compute O - C"))
        self.groupBox.setTitle(_translate("OCWidget", "Ephemeris and Period Correction"))
        self.label_2.setText(_translate("OCWidget", "ΔP"))
        self.label.setText(_translate("OCWidget", "ΔT"))
        self.calculate_btn.setText(_translate("OCWidget", "Calculate"))
        self.update_btn.setText(_translate("OCWidget", "Update"))
        self.dpdt_chk.setText(_translate("OCWidget", "Residuals with dP/dt"))
        self.linear_chk.setText(_translate("OCWidget", "Linear Residuals"))
        self.data_treewidget.headerItem().setText(0, _translate("OCWidget", "HJD"))
        self.data_treewidget.headerItem().setText(1, _translate("OCWidget", "Lin. Res."))
        self.data_treewidget.headerItem().setToolTip(1, _translate("OCWidget", "Linear residuals"))
        self.data_treewidget.headerItem().setText(2, _translate("OCWidget", "with dP/dt"))
        self.data_treewidget.headerItem().setToolTip(2, _translate("OCWidget", "Linear residuals with dP/dt"))
        self.compute_btn.setText(_translate("OCWidget", "Compute"))
        self.export_btn.setText(_translate("OCWidget", "Export"))

