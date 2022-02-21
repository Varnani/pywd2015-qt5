# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'curveproperties_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_CurvePropertiesDialog(object):
    def setupUi(self, CurvePropertiesDialog):
        CurvePropertiesDialog.setObjectName("CurvePropertiesDialog")
        CurvePropertiesDialog.resize(1000, 550)
        CurvePropertiesDialog.setMinimumSize(QtCore.QSize(1000, 550))
        CurvePropertiesDialog.setMaximumSize(QtCore.QSize(1000, 550))
        self.gridLayout_8 = QtWidgets.QGridLayout(CurvePropertiesDialog)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.whatsthis_btn = QtWidgets.QPushButton(CurvePropertiesDialog)
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        font.setKerning(True)
        self.whatsthis_btn.setFont(font)
        self.whatsthis_btn.setObjectName("whatsthis_btn")
        self.gridLayout_8.addWidget(self.whatsthis_btn, 0, 4, 1, 1)
        self.groupBox_4 = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox_4.setMaximumSize(QtCore.QSize(300, 16777215))
        self.groupBox_4.setObjectName("groupBox_4")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_16 = QtWidgets.QLabel(self.groupBox_4)
        self.label_16.setObjectName("label_16")
        self.gridLayout_4.addWidget(self.label_16, 0, 0, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.groupBox_4)
        self.label_17.setObjectName("label_17")
        self.gridLayout_4.addWidget(self.label_17, 0, 1, 1, 1)
        self.ksd_spinbox = QtWidgets.QSpinBox(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ksd_spinbox.sizePolicy().hasHeightForWidth())
        self.ksd_spinbox.setSizePolicy(sizePolicy)
        self.ksd_spinbox.setButtonSymbols(QtWidgets.QAbstractSpinBox.UpDownArrows)
        self.ksd_spinbox.setMinimum(0)
        self.ksd_spinbox.setMaximum(2)
        self.ksd_spinbox.setProperty("value", 1)
        self.ksd_spinbox.setObjectName("ksd_spinbox")
        self.gridLayout_4.addWidget(self.ksd_spinbox, 1, 0, 1, 1)
        self.noise_combobox = QtWidgets.QComboBox(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.noise_combobox.sizePolicy().hasHeightForWidth())
        self.noise_combobox.setSizePolicy(sizePolicy)
        self.noise_combobox.setObjectName("noise_combobox")
        self.noise_combobox.addItem("")
        self.noise_combobox.addItem("")
        self.noise_combobox.addItem("")
        self.gridLayout_4.addWidget(self.noise_combobox, 1, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.groupBox_4)
        self.label_11.setObjectName("label_11")
        self.gridLayout_4.addWidget(self.label_11, 2, 0, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBox_4)
        self.label_10.setObjectName("label_10")
        self.gridLayout_4.addWidget(self.label_10, 2, 1, 1, 1)
        self.sigma_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_4)
        self.sigma_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.sigma_ipt.setDecimals(6)
        self.sigma_ipt.setMinimum(-99999.0)
        self.sigma_ipt.setMaximum(99999.0)
        self.sigma_ipt.setObjectName("sigma_ipt")
        self.gridLayout_4.addWidget(self.sigma_ipt, 3, 0, 1, 1)
        self.opsf_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_4)
        self.opsf_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.opsf_ipt.setDecimals(6)
        self.opsf_ipt.setMinimum(-99999.0)
        self.opsf_ipt.setMaximum(99999.0)
        self.opsf_ipt.setObjectName("opsf_ipt")
        self.gridLayout_4.addWidget(self.opsf_ipt, 3, 1, 1, 1)
        self.gridLayout_8.addWidget(self.groupBox_4, 3, 2, 1, 1)
        self.groupBox_6 = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox_6.setMaximumSize(QtCore.QSize(300, 16777215))
        self.groupBox_6.setObjectName("groupBox_6")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.groupBox_6)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.label_2 = QtWidgets.QLabel(self.groupBox_6)
        self.label_2.setObjectName("label_2")
        self.gridLayout_6.addWidget(self.label_2, 0, 0, 1, 1)
        self.bandpasscontextlist_btn = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bandpasscontextlist_btn.sizePolicy().hasHeightForWidth())
        self.bandpasscontextlist_btn.setSizePolicy(sizePolicy)
        self.bandpasscontextlist_btn.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.bandpasscontextlist_btn.setObjectName("bandpasscontextlist_btn")
        self.gridLayout_6.addWidget(self.bandpasscontextlist_btn, 0, 3, 2, 1)
        self.band_name_label = QtWidgets.QLabel(self.groupBox_6)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.band_name_label.setFont(font)
        self.band_name_label.setObjectName("band_name_label")
        self.gridLayout_6.addWidget(self.band_name_label, 1, 0, 1, 2, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
        self.band_box = QtWidgets.QSpinBox(self.groupBox_6)
        self.band_box.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.band_box.setMinimum(1)
        self.band_box.setMaximum(95)
        self.band_box.setProperty("value", 7)
        self.band_box.setObjectName("band_box")
        self.gridLayout_6.addWidget(self.band_box, 0, 1, 1, 1)
        self.line = QtWidgets.QFrame(self.groupBox_6)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout_6.addWidget(self.line, 0, 2, 2, 1)
        self.gridLayout_8.addWidget(self.groupBox_6, 2, 2, 1, 1)
        self.groupBox = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox.setMaximumSize(QtCore.QSize(300, 16777215))
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.l1_ipt = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.l1_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.l1_ipt.setDecimals(8)
        self.l1_ipt.setMinimum(-99999.0)
        self.l1_ipt.setMaximum(99999.0)
        self.l1_ipt.setObjectName("l1_ipt")
        self.gridLayout.addWidget(self.l1_ipt, 1, 0, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.groupBox)
        self.label_9.setObjectName("label_9")
        self.gridLayout.addWidget(self.label_9, 0, 2, 1, 1)
        self.l2_ipt = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.l2_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.l2_ipt.setDecimals(8)
        self.l2_ipt.setMinimum(-99999.0)
        self.l2_ipt.setMaximum(99999.0)
        self.l2_ipt.setObjectName("l2_ipt")
        self.gridLayout.addWidget(self.l2_ipt, 1, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 0, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 0, 1, 1, 1)
        self.l3_ipt = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.l3_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.l3_ipt.setDecimals(8)
        self.l3_ipt.setMinimum(-99999.0)
        self.l3_ipt.setMaximum(99999.0)
        self.l3_ipt.setObjectName("l3_ipt")
        self.gridLayout.addWidget(self.l3_ipt, 1, 2, 1, 1)
        self.gridLayout_8.addWidget(self.groupBox, 2, 1, 1, 1)
        self.title_label = QtWidgets.QLabel(CurvePropertiesDialog)
        self.title_label.setObjectName("title_label")
        self.gridLayout_8.addWidget(self.title_label, 0, 1, 1, 2)
        self.line_3 = QtWidgets.QFrame(CurvePropertiesDialog)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.gridLayout_8.addWidget(self.line_3, 1, 1, 1, 4)
        self.line_4 = QtWidgets.QFrame(CurvePropertiesDialog)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.gridLayout_8.addWidget(self.line_4, 0, 3, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox_2.setMaximumSize(QtCore.QSize(300, 16777215))
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_5 = QtWidgets.QLabel(self.groupBox_2)
        self.label_5.setObjectName("label_5")
        self.gridLayout_2.addWidget(self.label_5, 1, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox_2)
        self.label_6.setObjectName("label_6")
        self.gridLayout_2.addWidget(self.label_6, 1, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.groupBox_2)
        self.label_7.setObjectName("label_7")
        self.gridLayout_2.addWidget(self.label_7, 3, 0, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.groupBox_2)
        self.label_8.setObjectName("label_8")
        self.gridLayout_2.addWidget(self.label_8, 3, 1, 1, 1)
        self.x1_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.x1_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.x1_ipt.setDecimals(6)
        self.x1_ipt.setMinimum(-99999.0)
        self.x1_ipt.setMaximum(99999.0)
        self.x1_ipt.setObjectName("x1_ipt")
        self.gridLayout_2.addWidget(self.x1_ipt, 2, 0, 1, 1)
        self.x2_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.x2_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.x2_ipt.setDecimals(6)
        self.x2_ipt.setMinimum(-99999.0)
        self.x2_ipt.setMaximum(99999.0)
        self.x2_ipt.setObjectName("x2_ipt")
        self.gridLayout_2.addWidget(self.x2_ipt, 2, 1, 1, 1)
        self.y1_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.y1_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.y1_ipt.setDecimals(6)
        self.y1_ipt.setMinimum(-99999.0)
        self.y1_ipt.setMaximum(99999.0)
        self.y1_ipt.setObjectName("y1_ipt")
        self.gridLayout_2.addWidget(self.y1_ipt, 4, 0, 1, 1)
        self.y2_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.y2_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.y2_ipt.setDecimals(6)
        self.y2_ipt.setMinimum(-99999.0)
        self.y2_ipt.setMaximum(99999.0)
        self.y2_ipt.setObjectName("y2_ipt")
        self.gridLayout_2.addWidget(self.y2_ipt, 4, 1, 1, 1)
        self.gridLayout_8.addWidget(self.groupBox_2, 3, 1, 1, 1)
        self.line_5 = QtWidgets.QFrame(CurvePropertiesDialog)
        self.line_5.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.gridLayout_8.addWidget(self.line_5, 2, 3, 5, 1)
        self.groupBox_7 = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox_7.setObjectName("groupBox_7")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.groupBox_7)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.time_combobox = QtWidgets.QComboBox(self.groupBox_7)
        self.time_combobox.setObjectName("time_combobox")
        self.horizontalLayout_3.addWidget(self.time_combobox)
        self.obs_combobox = QtWidgets.QComboBox(self.groupBox_7)
        self.obs_combobox.setObjectName("obs_combobox")
        self.horizontalLayout_3.addWidget(self.obs_combobox)
        self.weight_combobox = QtWidgets.QComboBox(self.groupBox_7)
        self.weight_combobox.setObjectName("weight_combobox")
        self.horizontalLayout_3.addWidget(self.weight_combobox)
        self.gridLayout_7.addLayout(self.horizontalLayout_3, 3, 0, 1, 4)
        self.plot_btn = QtWidgets.QPushButton(self.groupBox_7)
        self.plot_btn.setObjectName("plot_btn")
        self.gridLayout_7.addWidget(self.plot_btn, 6, 0, 1, 2)
        self.line_2 = QtWidgets.QFrame(self.groupBox_7)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout_7.addWidget(self.line_2, 7, 0, 1, 4)
        self.repick_btn = QtWidgets.QPushButton(self.groupBox_7)
        self.repick_btn.setObjectName("repick_btn")
        self.gridLayout_7.addWidget(self.repick_btn, 6, 2, 1, 2)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_23 = QtWidgets.QLabel(self.groupBox_7)
        self.label_23.setObjectName("label_23")
        self.horizontalLayout_4.addWidget(self.label_23)
        self.label_24 = QtWidgets.QLabel(self.groupBox_7)
        self.label_24.setObjectName("label_24")
        self.horizontalLayout_4.addWidget(self.label_24)
        self.label_25 = QtWidgets.QLabel(self.groupBox_7)
        self.label_25.setObjectName("label_25")
        self.horizontalLayout_4.addWidget(self.label_25)
        self.gridLayout_7.addLayout(self.horizontalLayout_4, 2, 0, 1, 4)
        self.data_widget = QtWidgets.QTreeWidget(self.groupBox_7)
        self.data_widget.setObjectName("data_widget")
        self.data_widget.headerItem().setText(0, "1")
        self.gridLayout_7.addWidget(self.data_widget, 8, 0, 1, 4)
        self.line_7 = QtWidgets.QFrame(self.groupBox_7)
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.gridLayout_7.addWidget(self.line_7, 1, 0, 1, 4)
        self.filepath_label = QtWidgets.QLabel(self.groupBox_7)
        self.filepath_label.setText("")
        self.filepath_label.setObjectName("filepath_label")
        self.gridLayout_7.addWidget(self.filepath_label, 0, 0, 1, 4)
        self.gridLayout_8.addWidget(self.groupBox_7, 2, 4, 5, 1)
        self.groupBox_5 = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox_5.setMaximumSize(QtCore.QSize(300, 16777215))
        self.groupBox_5.setObjectName("groupBox_5")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_5)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label_19 = QtWidgets.QLabel(self.groupBox_5)
        self.label_19.setObjectName("label_19")
        self.gridLayout_5.addWidget(self.label_19, 0, 0, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.groupBox_5)
        self.label_18.setObjectName("label_18")
        self.gridLayout_5.addWidget(self.label_18, 0, 1, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.groupBox_5)
        self.label_21.setObjectName("label_21")
        self.gridLayout_5.addWidget(self.label_21, 2, 0, 1, 1)
        self.label_20 = QtWidgets.QLabel(self.groupBox_5)
        self.label_20.setObjectName("label_20")
        self.gridLayout_5.addWidget(self.label_20, 2, 1, 1, 1)
        self.wla_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_5)
        self.wla_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.wla_ipt.setDecimals(6)
        self.wla_ipt.setMinimum(-99999.0)
        self.wla_ipt.setMaximum(99999.0)
        self.wla_ipt.setObjectName("wla_ipt")
        self.gridLayout_5.addWidget(self.wla_ipt, 1, 0, 1, 1)
        self.aextinc_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_5)
        self.aextinc_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.aextinc_ipt.setDecimals(6)
        self.aextinc_ipt.setMinimum(-99999.0)
        self.aextinc_ipt.setMaximum(99999.0)
        self.aextinc_ipt.setObjectName("aextinc_ipt")
        self.gridLayout_5.addWidget(self.aextinc_ipt, 1, 1, 1, 1)
        self.xunit_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_5)
        self.xunit_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.xunit_ipt.setDecimals(6)
        self.xunit_ipt.setMinimum(-99999.0)
        self.xunit_ipt.setMaximum(99999.0)
        self.xunit_ipt.setObjectName("xunit_ipt")
        self.gridLayout_5.addWidget(self.xunit_ipt, 3, 0, 1, 1)
        self.calib_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_5)
        self.calib_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.calib_ipt.setDecimals(6)
        self.calib_ipt.setMinimum(-99999.0)
        self.calib_ipt.setMaximum(99999.0)
        self.calib_ipt.setObjectName("calib_ipt")
        self.gridLayout_5.addWidget(self.calib_ipt, 3, 1, 1, 1)
        self.gridLayout_8.addWidget(self.groupBox_5, 4, 2, 1, 1)
        self.groupBox_3 = QtWidgets.QGroupBox(CurvePropertiesDialog)
        self.groupBox_3.setMaximumSize(QtCore.QSize(300, 16777215))
        self.groupBox_3.setObjectName("groupBox_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_3)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_12 = QtWidgets.QLabel(self.groupBox_3)
        self.label_12.setObjectName("label_12")
        self.gridLayout_3.addWidget(self.label_12, 0, 0, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.groupBox_3)
        self.label_14.setObjectName("label_14")
        self.gridLayout_3.addWidget(self.label_14, 0, 1, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.groupBox_3)
        self.label_15.setObjectName("label_15")
        self.gridLayout_3.addWidget(self.label_15, 2, 0, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.groupBox_3)
        self.label_13.setObjectName("label_13")
        self.gridLayout_3.addWidget(self.label_13, 2, 1, 1, 1)
        self.e3_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.e3_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e3_ipt.setDecimals(4)
        self.e3_ipt.setProperty("value", 0.55)
        self.e3_ipt.setObjectName("e3_ipt")
        self.gridLayout_3.addWidget(self.e3_ipt, 1, 1, 1, 1)
        self.e1_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.e1_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e1_ipt.setDecimals(4)
        self.e1_ipt.setObjectName("e1_ipt")
        self.gridLayout_3.addWidget(self.e1_ipt, 1, 0, 1, 1)
        self.e2_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.e2_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e2_ipt.setDecimals(4)
        self.e2_ipt.setProperty("value", 0.05)
        self.e2_ipt.setObjectName("e2_ipt")
        self.gridLayout_3.addWidget(self.e2_ipt, 3, 0, 1, 1)
        self.e4_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.e4_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e4_ipt.setDecimals(4)
        self.e4_ipt.setProperty("value", 0.95)
        self.e4_ipt.setObjectName("e4_ipt")
        self.gridLayout_3.addWidget(self.e4_ipt, 3, 1, 1, 1)
        self.gridLayout_8.addWidget(self.groupBox_3, 4, 1, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.accept_btn = QtWidgets.QPushButton(CurvePropertiesDialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.accept_btn.sizePolicy().hasHeightForWidth())
        self.accept_btn.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.accept_btn.setFont(font)
        self.accept_btn.setObjectName("accept_btn")
        self.horizontalLayout.addWidget(self.accept_btn)
        self.line_6 = QtWidgets.QFrame(CurvePropertiesDialog)
        self.line_6.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.horizontalLayout.addWidget(self.line_6)
        self.discard_btn = QtWidgets.QPushButton(CurvePropertiesDialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.discard_btn.sizePolicy().hasHeightForWidth())
        self.discard_btn.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.discard_btn.setFont(font)
        self.discard_btn.setObjectName("discard_btn")
        self.horizontalLayout.addWidget(self.discard_btn)
        self.gridLayout_8.addLayout(self.horizontalLayout, 5, 1, 2, 2)

        self.retranslateUi(CurvePropertiesDialog)
        self.noise_combobox.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(CurvePropertiesDialog)
        CurvePropertiesDialog.setTabOrder(self.l1_ipt, self.l2_ipt)
        CurvePropertiesDialog.setTabOrder(self.l2_ipt, self.l3_ipt)
        CurvePropertiesDialog.setTabOrder(self.l3_ipt, self.band_box)
        CurvePropertiesDialog.setTabOrder(self.band_box, self.bandpasscontextlist_btn)
        CurvePropertiesDialog.setTabOrder(self.bandpasscontextlist_btn, self.x1_ipt)
        CurvePropertiesDialog.setTabOrder(self.x1_ipt, self.x2_ipt)
        CurvePropertiesDialog.setTabOrder(self.x2_ipt, self.ksd_spinbox)
        CurvePropertiesDialog.setTabOrder(self.ksd_spinbox, self.noise_combobox)
        CurvePropertiesDialog.setTabOrder(self.noise_combobox, self.y1_ipt)
        CurvePropertiesDialog.setTabOrder(self.y1_ipt, self.y2_ipt)
        CurvePropertiesDialog.setTabOrder(self.y2_ipt, self.sigma_ipt)
        CurvePropertiesDialog.setTabOrder(self.sigma_ipt, self.opsf_ipt)
        CurvePropertiesDialog.setTabOrder(self.opsf_ipt, self.e1_ipt)
        CurvePropertiesDialog.setTabOrder(self.e1_ipt, self.e3_ipt)
        CurvePropertiesDialog.setTabOrder(self.e3_ipt, self.wla_ipt)
        CurvePropertiesDialog.setTabOrder(self.wla_ipt, self.aextinc_ipt)
        CurvePropertiesDialog.setTabOrder(self.aextinc_ipt, self.e2_ipt)
        CurvePropertiesDialog.setTabOrder(self.e2_ipt, self.e4_ipt)
        CurvePropertiesDialog.setTabOrder(self.e4_ipt, self.xunit_ipt)
        CurvePropertiesDialog.setTabOrder(self.xunit_ipt, self.calib_ipt)
        CurvePropertiesDialog.setTabOrder(self.calib_ipt, self.time_combobox)
        CurvePropertiesDialog.setTabOrder(self.time_combobox, self.obs_combobox)
        CurvePropertiesDialog.setTabOrder(self.obs_combobox, self.weight_combobox)
        CurvePropertiesDialog.setTabOrder(self.weight_combobox, self.plot_btn)
        CurvePropertiesDialog.setTabOrder(self.plot_btn, self.repick_btn)
        CurvePropertiesDialog.setTabOrder(self.repick_btn, self.data_widget)
        CurvePropertiesDialog.setTabOrder(self.data_widget, self.accept_btn)
        CurvePropertiesDialog.setTabOrder(self.accept_btn, self.discard_btn)
        CurvePropertiesDialog.setTabOrder(self.discard_btn, self.whatsthis_btn)

    def retranslateUi(self, CurvePropertiesDialog):
        _translate = QtCore.QCoreApplication.translate
        CurvePropertiesDialog.setWindowTitle(_translate("CurvePropertiesDialog", "Curve Properties"))
        self.whatsthis_btn.setText(_translate("CurvePropertiesDialog", "?"))
        self.groupBox_4.setTitle(_translate("CurvePropertiesDialog", "Scatter and Attenuation"))
        self.label_16.setToolTip(_translate("CurvePropertiesDialog", "Set standard deviation apply method. [?]"))
        self.label_16.setWhatsThis(_translate("CurvePropertiesDialog", "<html><head/><body><p>An integer array that is 0, 1, or 2 for each input sub-dataset (velocity, light, or eclipse timings). </p><p>The KSDs tell DC whether to apply the input standard deviations (σ’s) to compute curvedependent weights (KSD=0), </p><p>to apply DC’s internally computed σ’s for the weights (KSD=1),</p><p>or to apply σ’s based on one or two restricted phase ranges for the weights (KSD=2).</p><p><span style=\" font-weight:600;\">If unsure, set to 1.</span></p></body></html>"))
        self.label_16.setText(_translate("CurvePropertiesDialog", "KSD"))
        self.label_17.setToolTip(_translate("CurvePropertiesDialog", "Set how observational scatter scales with the light level. [?]"))
        self.label_17.setWhatsThis(_translate("CurvePropertiesDialog", "<html><head/><body><p><span style=\" font-weight:600;\">NOISE</span> should be set to 1 for scatter that scales with the square root of the light level, such as counting statistics, and to 2 for scatter that scales with the light level, such as scintillation noise or fluctuations in sky transparency. If <span style=\" font-weight:600;\">NOISE</span> is set to 0, no level-dependent weighting is applied. </p><p>As a side note, level-dependent weighting does not apply to velocity curves.</p></body></html>"))
        self.label_17.setText(_translate("CurvePropertiesDialog", "NOISE"))
        self.noise_combobox.setItemText(0, _translate("CurvePropertiesDialog", "None"))
        self.noise_combobox.setItemText(1, _translate("CurvePropertiesDialog", "Square Root"))
        self.noise_combobox.setItemText(2, _translate("CurvePropertiesDialog", "Linear"))
        self.label_11.setToolTip(_translate("CurvePropertiesDialog", "Estimated standard deviation of observed light"))
        self.label_11.setText(_translate("CurvePropertiesDialog", "SIGMA"))
        self.label_10.setToolTip(_translate("CurvePropertiesDialog", "Opacity, attenuation of star light by circumstellar matter"))
        self.label_10.setText(_translate("CurvePropertiesDialog", "OPSF"))
        self.groupBox_6.setTitle(_translate("CurvePropertiesDialog", "Band"))
        self.label_2.setToolTip(_translate("CurvePropertiesDialog", "Band identification number. Check [Bandpass List] table for all available bands and their respective #\'s"))
        self.label_2.setText(_translate("CurvePropertiesDialog", "ID"))
        self.bandpasscontextlist_btn.setText(_translate("CurvePropertiesDialog", "List"))
        self.band_name_label.setText(_translate("CurvePropertiesDialog", "Johnson V"))
        self.groupBox.setTitle(_translate("CurvePropertiesDialog", "Luminosity"))
        self.label_9.setToolTip(_translate("CurvePropertiesDialog", "Third light"))
        self.label_9.setText(_translate("CurvePropertiesDialog", "L3"))
        self.label_3.setToolTip(_translate("CurvePropertiesDialog", "Bandpass luminosity for star 1 [HLA]"))
        self.label_3.setText(_translate("CurvePropertiesDialog", "L1"))
        self.label_4.setToolTip(_translate("CurvePropertiesDialog", "Bandpass luminosity for star 2 [CLA]"))
        self.label_4.setText(_translate("CurvePropertiesDialog", "L2"))
        self.title_label.setText(_translate("CurvePropertiesDialog", "TITLE/////////////////"))
        self.groupBox_2.setTitle(_translate("CurvePropertiesDialog", "Limb Darkening"))
        self.label_5.setToolTip(_translate("CurvePropertiesDialog", "Wavelength-specific limb darkening coefficient in linear term, for star 1 [X1A] "))
        self.label_5.setText(_translate("CurvePropertiesDialog", "X1"))
        self.label_6.setToolTip(_translate("CurvePropertiesDialog", "Wavelength-specific limb darkening coefficient in linear term, for star 2 [X2A] "))
        self.label_6.setText(_translate("CurvePropertiesDialog", "X2"))
        self.label_7.setToolTip(_translate("CurvePropertiesDialog", "Bandpass-specific limb darkening coefficient in non-linear term, for star 1 [Y1A] "))
        self.label_7.setText(_translate("CurvePropertiesDialog", "Y1"))
        self.label_8.setToolTip(_translate("CurvePropertiesDialog", "Bandpass-specific limb darkening coefficient in non-linear term, for star 2 [Y2A] "))
        self.label_8.setText(_translate("CurvePropertiesDialog", "Y2"))
        self.groupBox_7.setTitle(_translate("CurvePropertiesDialog", "Data Preview"))
        self.plot_btn.setText(_translate("CurvePropertiesDialog", "Plot"))
        self.repick_btn.setText(_translate("CurvePropertiesDialog", "Pick Again"))
        self.label_23.setText(_translate("CurvePropertiesDialog", "Time / Phase"))
        self.label_24.setText(_translate("CurvePropertiesDialog", "Observation"))
        self.label_25.setText(_translate("CurvePropertiesDialog", "Weight"))
        self.groupBox_5.setTitle(_translate("CurvePropertiesDialog", "Misc."))
        self.label_19.setToolTip(_translate("CurvePropertiesDialog", "<html><head/><body><p>Observational wavelenght in microns. </p><p>Wavelengths are no longer used for light curves, which are based on integrated bandpass radiation since the 2003 revision.</p><p>They are still entered for use as reference wavelengths for line profiles and for opacity computations in circumstellar attenuation.</p></body></html>"))
        self.label_19.setText(_translate("CurvePropertiesDialog", "WLA"))
        self.label_18.setToolTip(_translate("CurvePropertiesDialog", "Interstellar extinction in magnitude in the designated [LINKEXT] photometric band. [?]"))
        self.label_18.setWhatsThis(_translate("CurvePropertiesDialog", "<html><head/><body><p>AEXTINC refers to a definite photometric band (the band designated by control integer LINKEXT), so it has only one value and it is not band-dependent.</p></body></html>"))
        self.label_18.setText(_translate("CurvePropertiesDialog", "AEXTINC"))
        self.label_21.setToolTip(_translate("CurvePropertiesDialog", "The multiplier that restores the flux numbers to directly inverted magnitudes. [?]"))
        self.label_21.setWhatsThis(_translate("CurvePropertiesDialog", "<html><head/><body><p>This quantity is intended only for flux (as opposed to magnitude) input, and allows rescaling of observed fluxes so that the input flux numbers can be in a convenient range. If input fluxes are simply inverted magnitudes (e.g. 10<span style=\" vertical-align:super;\">−0.4V</span> for V magnitudes), then XUNIT should be unity. However, it may be convenient to enter rescaled (i.e. wrongly scaled) fluxes in cases where they would otherwise be of order, say, 10<span style=\" vertical-align:super;\">−10</span> or 10<span style=\" vertical-align:super;\">3</span> . Then XUNIT can be the multiplier that restores the flux numbers to directly inverted magnitudes. If working in (presumably standard) magnitudes with control integer MAGLITE=1, of course, then don’t worry about XUNIT—just set it to 1.0000.</p></body></html>"))
        self.label_21.setText(_translate("CurvePropertiesDialog", "XUNIT"))
        self.label_20.setToolTip(_translate("CurvePropertiesDialog", "Flux calibration constant in erg s −1 cm −3 for a star of magnitude 0.00.  "))
        self.label_20.setText(_translate("CurvePropertiesDialog", "CALIB"))
        self.groupBox_3.setTitle(_translate("CurvePropertiesDialog", "Phase Windows for Std. Dev. Calculation"))
        self.label_12.setToolTip(_translate("CurvePropertiesDialog", "Phase 1"))
        self.label_12.setText(_translate("CurvePropertiesDialog", "E1"))
        self.label_14.setToolTip(_translate("CurvePropertiesDialog", "Phase 3"))
        self.label_14.setText(_translate("CurvePropertiesDialog", "E3"))
        self.label_15.setToolTip(_translate("CurvePropertiesDialog", "Phase 2"))
        self.label_15.setText(_translate("CurvePropertiesDialog", "E2"))
        self.label_13.setToolTip(_translate("CurvePropertiesDialog", "Phase 4"))
        self.label_13.setText(_translate("CurvePropertiesDialog", "E4"))
        self.accept_btn.setText(_translate("CurvePropertiesDialog", "Accept"))
        self.discard_btn.setText(_translate("CurvePropertiesDialog", "Discard"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    CurvePropertiesDialog = QtWidgets.QDialog()
    ui = Ui_CurvePropertiesDialog()
    ui.setupUi(CurvePropertiesDialog)
    CurvePropertiesDialog.show()
    sys.exit(app.exec_())

