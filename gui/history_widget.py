# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'history_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_HistoryWidget(object):
    def setupUi(self, HistoryWidget):
        HistoryWidget.setObjectName("HistoryWidget")
        HistoryWidget.resize(900, 500)
        self.gridLayout_3 = QtWidgets.QGridLayout(HistoryWidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.plot_btn = QtWidgets.QPushButton(HistoryWidget)
        self.plot_btn.setObjectName("plot_btn")
        self.gridLayout_2.addWidget(self.plot_btn, 0, 0, 1, 1)
        self.auto_chk = QtWidgets.QCheckBox(HistoryWidget)
        self.auto_chk.setObjectName("auto_chk")
        self.gridLayout_2.addWidget(self.auto_chk, 0, 1, 1, 1)
        self.clear_btn = QtWidgets.QPushButton(HistoryWidget)
        self.clear_btn.setObjectName("clear_btn")
        self.gridLayout_2.addWidget(self.clear_btn, 1, 0, 1, 1)
        self.export_btn = QtWidgets.QPushButton(HistoryWidget)
        self.export_btn.setObjectName("export_btn")
        self.gridLayout_2.addWidget(self.export_btn, 1, 1, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 2, 0, 1, 1)
        self.history_treewidget = QtWidgets.QTreeWidget(HistoryWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.history_treewidget.sizePolicy().hasHeightForWidth())
        self.history_treewidget.setSizePolicy(sizePolicy)
        self.history_treewidget.setMinimumSize(QtCore.QSize(400, 0))
        self.history_treewidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectColumns)
        self.history_treewidget.setObjectName("history_treewidget")
        self.gridLayout.addWidget(self.history_treewidget, 0, 0, 1, 1)
        self.line_2 = QtWidgets.QFrame(HistoryWidget)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout.addWidget(self.line_2, 1, 0, 1, 1)
        self.plot_widget = QtWidgets.QWidget(HistoryWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plot_widget.sizePolicy().hasHeightForWidth())
        self.plot_widget.setSizePolicy(sizePolicy)
        self.plot_widget.setMinimumSize(QtCore.QSize(300, 300))
        self.plot_widget.setObjectName("plot_widget")
        self.gridLayout.addWidget(self.plot_widget, 0, 2, 3, 1)
        self.line = QtWidgets.QFrame(HistoryWidget)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 0, 1, 3, 1)
        self.gridLayout_3.addLayout(self.gridLayout, 0, 0, 1, 1)

        self.retranslateUi(HistoryWidget)
        QtCore.QMetaObject.connectSlotsByName(HistoryWidget)

    def retranslateUi(self, HistoryWidget):
        _translate = QtCore.QCoreApplication.translate
        HistoryWidget.setWindowTitle(_translate("HistoryWidget", "Solution History"))
        self.plot_btn.setText(_translate("HistoryWidget", "Plot"))
        self.auto_chk.setToolTip(_translate("HistoryWidget", "Automatically plot selected parameter when DC program finishes its calculation"))
        self.auto_chk.setText(_translate("HistoryWidget", "Auto"))
        self.clear_btn.setText(_translate("HistoryWidget", "Clear"))
        self.export_btn.setText(_translate("HistoryWidget", "Export All"))

