# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'configurespots_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_SpotConfigureWidget(object):
    def setupUi(self, SpotConfigureWidget):
        SpotConfigureWidget.setObjectName("SpotConfigureWidget")
        SpotConfigureWidget.resize(900, 650)
        SpotConfigureWidget.setMinimumSize(QtCore.QSize(900, 650))
        self.gridLayout_4 = QtWidgets.QGridLayout(SpotConfigureWidget)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.line = QtWidgets.QFrame(SpotConfigureWidget)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout_4.addWidget(self.line, 2, 3, 2, 1)
        self.spotconfigload_btn = QtWidgets.QPushButton(SpotConfigureWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.spotconfigload_btn.sizePolicy().hasHeightForWidth())
        self.spotconfigload_btn.setSizePolicy(sizePolicy)
        self.spotconfigload_btn.setObjectName("spotconfigload_btn")
        self.gridLayout_4.addWidget(self.spotconfigload_btn, 3, 6, 1, 1)
        self.whatsthis_btn = QtWidgets.QPushButton(SpotConfigureWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.whatsthis_btn.sizePolicy().hasHeightForWidth())
        self.whatsthis_btn.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.whatsthis_btn.setFont(font)
        self.whatsthis_btn.setObjectName("whatsthis_btn")
        self.gridLayout_4.addWidget(self.whatsthis_btn, 2, 8, 2, 1)
        self.kspev_groupbox = QtWidgets.QGroupBox(SpotConfigureWidget)
        self.kspev_groupbox.setCheckable(True)
        self.kspev_groupbox.setObjectName("kspev_groupbox")
        self.gridLayout = QtWidgets.QGridLayout(self.kspev_groupbox)
        self.gridLayout.setObjectName("gridLayout")
        self.label_7 = QtWidgets.QLabel(self.kspev_groupbox)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 0, 0, 1, 1)
        self.nomax_combobox = QtWidgets.QComboBox(self.kspev_groupbox)
        self.nomax_combobox.setObjectName("nomax_combobox")
        self.nomax_combobox.addItem("")
        self.nomax_combobox.addItem("")
        self.gridLayout.addWidget(self.nomax_combobox, 0, 1, 1, 1)
        self.gridLayout_4.addWidget(self.kspev_groupbox, 2, 2, 2, 1)
        self.line_4 = QtWidgets.QFrame(SpotConfigureWidget)
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.gridLayout_4.addWidget(self.line_4, 5, 0, 1, 9)
        self.line_2 = QtWidgets.QFrame(SpotConfigureWidget)
        self.line_2.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout_4.addWidget(self.line_2, 2, 5, 2, 1)
        self.groupBox_3 = QtWidgets.QGroupBox(SpotConfigureWidget)
        self.groupBox_3.setObjectName("groupBox_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_3)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_5 = QtWidgets.QLabel(self.groupBox_3)
        self.label_5.setToolTip("")
        self.label_5.setObjectName("label_5")
        self.gridLayout_3.addWidget(self.label_5, 0, 0, 1, 1)
        self.fspot1_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.fspot1_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.fspot1_ipt.setDecimals(4)
        self.fspot1_ipt.setMaximum(99.0)
        self.fspot1_ipt.setSingleStep(1.0)
        self.fspot1_ipt.setProperty("value", 1.0)
        self.fspot1_ipt.setObjectName("fspot1_ipt")
        self.gridLayout_3.addWidget(self.fspot1_ipt, 0, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox_3)
        self.label_6.setToolTip("")
        self.label_6.setObjectName("label_6")
        self.gridLayout_3.addWidget(self.label_6, 1, 0, 1, 1)
        self.fspot2_ipt = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.fspot2_ipt.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.fspot2_ipt.setDecimals(4)
        self.fspot2_ipt.setMaximum(99.0)
        self.fspot2_ipt.setSingleStep(1.0)
        self.fspot2_ipt.setProperty("value", 1.0)
        self.fspot2_ipt.setObjectName("fspot2_ipt")
        self.gridLayout_3.addWidget(self.fspot2_ipt, 1, 1, 1, 1)
        self.gridLayout_4.addWidget(self.groupBox_3, 2, 0, 2, 1)
        self.star2_groupbox = QtWidgets.QGroupBox(SpotConfigureWidget)
        self.star2_groupbox.setObjectName("star2_groupbox")
        self.gridLayout_4.addWidget(self.star2_groupbox, 7, 0, 1, 9)
        self.spotconfigsave_btn = QtWidgets.QPushButton(SpotConfigureWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.spotconfigsave_btn.sizePolicy().hasHeightForWidth())
        self.spotconfigsave_btn.setSizePolicy(sizePolicy)
        self.spotconfigsave_btn.setObjectName("spotconfigsave_btn")
        self.gridLayout_4.addWidget(self.spotconfigsave_btn, 2, 6, 1, 1)
        self.line_3 = QtWidgets.QFrame(SpotConfigureWidget)
        self.line_3.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.gridLayout_4.addWidget(self.line_3, 2, 7, 2, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(SpotConfigureWidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.ifsmv1_chk = QtWidgets.QCheckBox(self.groupBox_2)
        self.ifsmv1_chk.setObjectName("ifsmv1_chk")
        self.gridLayout_2.addWidget(self.ifsmv1_chk, 0, 0, 1, 1)
        self.ifsmv2_chk = QtWidgets.QCheckBox(self.groupBox_2)
        self.ifsmv2_chk.setObjectName("ifsmv2_chk")
        self.gridLayout_2.addWidget(self.ifsmv2_chk, 1, 0, 1, 1)
        self.gridLayout_4.addWidget(self.groupBox_2, 2, 1, 2, 1)
        self.label_17 = QtWidgets.QLabel(SpotConfigureWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy)
        self.label_17.setObjectName("label_17")
        self.gridLayout_4.addWidget(self.label_17, 0, 0, 1, 5)
        self.line_5 = QtWidgets.QFrame(SpotConfigureWidget)
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.gridLayout_4.addWidget(self.line_5, 1, 0, 1, 9)
        self.star1_groupbox = QtWidgets.QGroupBox(SpotConfigureWidget)
        self.star1_groupbox.setObjectName("star1_groupbox")
        self.gridLayout_4.addWidget(self.star1_groupbox, 6, 0, 1, 9)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.kspot_chk = QtWidgets.QCheckBox(SpotConfigureWidget)
        self.kspot_chk.setObjectName("kspot_chk")
        self.verticalLayout.addWidget(self.kspot_chk)
        self.line_6 = QtWidgets.QFrame(SpotConfigureWidget)
        self.line_6.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.verticalLayout.addWidget(self.line_6)
        self.clear_a_b_button = QtWidgets.QPushButton(SpotConfigureWidget)
        self.clear_a_b_button.setObjectName("clear_a_b_button")
        self.verticalLayout.addWidget(self.clear_a_b_button)
        self.gridLayout_4.addLayout(self.verticalLayout, 2, 4, 2, 1)

        self.retranslateUi(SpotConfigureWidget)
        QtCore.QMetaObject.connectSlotsByName(SpotConfigureWidget)

    def retranslateUi(self, SpotConfigureWidget):
        _translate = QtCore.QCoreApplication.translate
        SpotConfigureWidget.setWindowTitle(_translate("SpotConfigureWidget", "Configure Spots"))
        self.spotconfigload_btn.setText(_translate("SpotConfigureWidget", "Load Spots"))
        self.whatsthis_btn.setText(_translate("SpotConfigureWidget", "?"))
        self.kspev_groupbox.setToolTip(_translate("SpotConfigureWidget", "<html><head/><body><p>Enable spot aging (KSPEV) [?]</p></body></html>"))
        self.kspev_groupbox.setWhatsThis(_translate("SpotConfigureWidget", "<html><head/><body><p>Controls whether spots age (grow and decay) in radius. Currently there is no aging in spot temperature. Uncheck for no aging, check for aging. Solutions for spot aging need good starting estimates for spot parameters and careful monitoring of solution progress.</p></body></html>"))
        self.kspev_groupbox.setTitle(_translate("SpotConfigureWidget", "Spot Aging"))
        self.label_7.setToolTip(_translate("SpotConfigureWidget", "<html><head/><body><p>Spot aging profile (NOMAX) [?]</p></body></html>"))
        self.label_7.setWhatsThis(_translate("SpotConfigureWidget", "<html><head/><body><p>Tells whether the spot growth and decay timewise profile is trapezoidal or triangular. Setting this to triangular eliminates the interval of constant size that otherwise exists at spot maximum.</p></body></html>"))
        self.label_7.setText(_translate("SpotConfigureWidget", "Aging Profile"))
        self.nomax_combobox.setItemText(0, _translate("SpotConfigureWidget", "Triangular"))
        self.nomax_combobox.setItemText(1, _translate("SpotConfigureWidget", "Trapezoidal"))
        self.groupBox_3.setToolTip(_translate("SpotConfigureWidget", "<html><head/><body><p>Spot angular drift rate in longitude for star 1 and 2, rate 1.000 means that drift matches the mean orbital angular rate (IFSMV2)</p></body></html>"))
        self.groupBox_3.setTitle(_translate("SpotConfigureWidget", "Drift Rates"))
        self.label_5.setText(_translate("SpotConfigureWidget", "Star 1"))
        self.label_6.setText(_translate("SpotConfigureWidget", "Star 2"))
        self.star2_groupbox.setTitle(_translate("SpotConfigureWidget", "Star 2"))
        self.spotconfigsave_btn.setText(_translate("SpotConfigureWidget", "Save Spots"))
        self.groupBox_2.setTitle(_translate("SpotConfigureWidget", "Spot Movement"))
        self.ifsmv1_chk.setToolTip(_translate("SpotConfigureWidget", "<html><head/><body><p>Allow spot A to move in longitude (IFSMV1)</p></body></html>"))
        self.ifsmv1_chk.setText(_translate("SpotConfigureWidget", "Allow for Spot A"))
        self.ifsmv2_chk.setToolTip(_translate("SpotConfigureWidget", "<html><head/><body><p>Allow spot B to move in longitude (IFSMV2)</p></body></html>"))
        self.ifsmv2_chk.setText(_translate("SpotConfigureWidget", "Allow for Spot B"))
        self.label_17.setText(_translate("SpotConfigureWidget", "Add spots and configure spot realeted parameters"))
        self.star1_groupbox.setTitle(_translate("SpotConfigureWidget", "Star 1"))
        self.kspot_chk.setToolTip(_translate("SpotConfigureWidget", "<html><head/><body><p>Use &quot;Vector Fractional Area&quot; algorithm (KSPOT) [?]</p></body></html>"))
        self.kspot_chk.setWhatsThis(_translate("SpotConfigureWidget", "<html><head/><body><p>Controls whether the old simple spot algorithm or the much more precise Vector Fractional Area algorithm <span style=\" font-weight:600;\">(Wilson 2012b)</span> is applied.</p></body></html>"))
        self.kspot_chk.setText(_translate("SpotConfigureWidget", "Use \"VFA\""))
        self.clear_a_b_button.setText(_translate("SpotConfigureWidget", "Clear A and B"))
