# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dimensions_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DimensionWidget(object):
    def setupUi(self, DimensionWidget):
        DimensionWidget.setObjectName("DimensionWidget")
        DimensionWidget.resize(1000, 600)
        self.gridLayout_5 = QtWidgets.QGridLayout(DimensionWidget)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.groupBox_2 = QtWidgets.QGroupBox(DimensionWidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout.setObjectName("verticalLayout")
        self.s1_plot_widget = QtWidgets.QWidget(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s1_plot_widget.sizePolicy().hasHeightForWidth())
        self.s1_plot_widget.setSizePolicy(sizePolicy)
        self.s1_plot_widget.setMinimumSize(QtCore.QSize(300, 200))
        self.s1_plot_widget.setObjectName("s1_plot_widget")
        self.verticalLayout.addWidget(self.s1_plot_widget)
        self.line_2 = QtWidgets.QFrame(self.groupBox_2)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.verticalLayout.addWidget(self.line_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.s1_pole_chk = QtWidgets.QCheckBox(self.groupBox_2)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 48, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 48, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(190, 190, 190))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        self.s1_pole_chk.setPalette(palette)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s1_pole_chk.setFont(font)
        self.s1_pole_chk.setChecked(True)
        self.s1_pole_chk.setObjectName("s1_pole_chk")
        self.horizontalLayout.addWidget(self.s1_pole_chk)
        self.s1_point_chk = QtWidgets.QCheckBox(self.groupBox_2)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s1_point_chk.setFont(font)
        self.s1_point_chk.setObjectName("s1_point_chk")
        self.horizontalLayout.addWidget(self.s1_point_chk)
        self.s1_side_chk = QtWidgets.QCheckBox(self.groupBox_2)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(239, 41, 41))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(239, 41, 41))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(190, 190, 190))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        self.s1_side_chk.setPalette(palette)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s1_side_chk.setFont(font)
        self.s1_side_chk.setObjectName("s1_side_chk")
        self.horizontalLayout.addWidget(self.s1_side_chk)
        self.s1_back_chk = QtWidgets.QCheckBox(self.groupBox_2)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(115, 210, 22))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(115, 210, 22))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(190, 190, 190))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        self.s1_back_chk.setPalette(palette)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s1_back_chk.setFont(font)
        self.s1_back_chk.setObjectName("s1_back_chk")
        self.horizontalLayout.addWidget(self.s1_back_chk)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.gridLayout_5.addWidget(self.groupBox_2, 0, 0, 1, 1)
        self.groupBox = QtWidgets.QGroupBox(DimensionWidget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.s2_pole_chk = QtWidgets.QCheckBox(self.groupBox)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 48, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 48, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(190, 190, 190))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        self.s2_pole_chk.setPalette(palette)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s2_pole_chk.setFont(font)
        self.s2_pole_chk.setChecked(True)
        self.s2_pole_chk.setObjectName("s2_pole_chk")
        self.horizontalLayout_2.addWidget(self.s2_pole_chk)
        self.s2_point_chk = QtWidgets.QCheckBox(self.groupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s2_point_chk.setFont(font)
        self.s2_point_chk.setObjectName("s2_point_chk")
        self.horizontalLayout_2.addWidget(self.s2_point_chk)
        self.s2_side_chk = QtWidgets.QCheckBox(self.groupBox)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(239, 41, 41))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(239, 41, 41))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(190, 190, 190))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        self.s2_side_chk.setPalette(palette)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s2_side_chk.setFont(font)
        self.s2_side_chk.setObjectName("s2_side_chk")
        self.horizontalLayout_2.addWidget(self.s2_side_chk)
        self.s2_back_chk = QtWidgets.QCheckBox(self.groupBox)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(115, 210, 22))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(115, 210, 22))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(190, 190, 190))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        self.s2_back_chk.setPalette(palette)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.s2_back_chk.setFont(font)
        self.s2_back_chk.setObjectName("s2_back_chk")
        self.horizontalLayout_2.addWidget(self.s2_back_chk)
        self.gridLayout_3.addLayout(self.horizontalLayout_2, 2, 0, 1, 1)
        self.s2_plot_widget = QtWidgets.QWidget(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s2_plot_widget.sizePolicy().hasHeightForWidth())
        self.s2_plot_widget.setSizePolicy(sizePolicy)
        self.s2_plot_widget.setMinimumSize(QtCore.QSize(300, 200))
        self.s2_plot_widget.setObjectName("s2_plot_widget")
        self.gridLayout_3.addWidget(self.s2_plot_widget, 0, 0, 1, 1)
        self.line_3 = QtWidgets.QFrame(self.groupBox)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.gridLayout_3.addWidget(self.line_3, 1, 0, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox, 0, 1, 1, 1)
        self.line = QtWidgets.QFrame(DimensionWidget)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout_5.addWidget(self.line, 3, 0, 1, 2)
        self.plot_btn = QtWidgets.QPushButton(DimensionWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plot_btn.sizePolicy().hasHeightForWidth())
        self.plot_btn.setSizePolicy(sizePolicy)
        self.plot_btn.setObjectName("plot_btn")
        self.gridLayout_5.addWidget(self.plot_btn, 4, 0, 1, 2)

        self.retranslateUi(DimensionWidget)
        QtCore.QMetaObject.connectSlotsByName(DimensionWidget)

    def retranslateUi(self, DimensionWidget):
        _translate = QtCore.QCoreApplication.translate
        DimensionWidget.setWindowTitle(_translate("DimensionWidget", "Component Dimensions"))
        self.groupBox_2.setTitle(_translate("DimensionWidget", "Star 1"))
        self.s1_pole_chk.setText(_translate("DimensionWidget", "Pole"))
        self.s1_point_chk.setText(_translate("DimensionWidget", "Point"))
        self.s1_side_chk.setText(_translate("DimensionWidget", "Side"))
        self.s1_back_chk.setText(_translate("DimensionWidget", "Back"))
        self.groupBox.setTitle(_translate("DimensionWidget", "Star 2"))
        self.s2_pole_chk.setText(_translate("DimensionWidget", "Pole"))
        self.s2_point_chk.setText(_translate("DimensionWidget", "Point"))
        self.s2_side_chk.setText(_translate("DimensionWidget", "Side"))
        self.s2_back_chk.setText(_translate("DimensionWidget", "Back"))
        self.plot_btn.setText(_translate("DimensionWidget", "Calculate Component Radii"))
