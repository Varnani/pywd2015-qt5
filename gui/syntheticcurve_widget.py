# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'syntheticcurve_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_SyntheticCurveWidget(object):
    def setupUi(self, SyntheticCurveWidget):
        SyntheticCurveWidget.setObjectName("SyntheticCurveWidget")
        SyntheticCurveWidget.resize(1000, 600)
        self.gridLayout_4 = QtWidgets.QGridLayout(SyntheticCurveWidget)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.tabWidget = QtWidgets.QTabWidget(SyntheticCurveWidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.groupBox_3 = QtWidgets.QGroupBox(self.tab)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.light_plotmodel_chk = QtWidgets.QCheckBox(self.groupBox_3)
        self.light_plotmodel_chk.setChecked(True)
        self.light_plotmodel_chk.setObjectName("light_plotmodel_chk")
        self.verticalLayout_2.addWidget(self.light_plotmodel_chk)
        self.light_plotobs_chk = QtWidgets.QCheckBox(self.groupBox_3)
        self.light_plotobs_chk.setChecked(True)
        self.light_plotobs_chk.setObjectName("light_plotobs_chk")
        self.verticalLayout_2.addWidget(self.light_plotobs_chk)
        self.light_alias_chk = QtWidgets.QCheckBox(self.groupBox_3)
        self.light_alias_chk.setObjectName("light_alias_chk")
        self.verticalLayout_2.addWidget(self.light_alias_chk)
        self.gridLayout_2.addWidget(self.groupBox_3, 0, 2, 1, 2)
        self.groupBox_4 = QtWidgets.QGroupBox(self.tab)
        self.groupBox_4.setObjectName("groupBox_4")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout.setObjectName("gridLayout")
        self.light_treewidget = QtWidgets.QTreeWidget(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.light_treewidget.sizePolicy().hasHeightForWidth())
        self.light_treewidget.setSizePolicy(sizePolicy)
        self.light_treewidget.setMaximumSize(QtCore.QSize(16777215, 175))
        self.light_treewidget.setAlternatingRowColors(True)
        self.light_treewidget.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.light_treewidget.setObjectName("light_treewidget")
        self.light_treewidget.header().setDefaultSectionSize(53)
        self.light_treewidget.header().setMinimumSectionSize(30)
        self.light_treewidget.header().setStretchLastSection(True)
        self.gridLayout.addWidget(self.light_treewidget, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.groupBox_4, 0, 0, 3, 1)
        self.light_plot_btn = QtWidgets.QPushButton(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.light_plot_btn.sizePolicy().hasHeightForWidth())
        self.light_plot_btn.setSizePolicy(sizePolicy)
        self.light_plot_btn.setObjectName("light_plot_btn")
        self.gridLayout_2.addWidget(self.light_plot_btn, 2, 2, 1, 2)
        self.line_3 = QtWidgets.QFrame(self.tab)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.gridLayout_2.addWidget(self.line_3, 3, 0, 1, 4)
        self.light_plot_widget = QtWidgets.QWidget(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.light_plot_widget.sizePolicy().hasHeightForWidth())
        self.light_plot_widget.setSizePolicy(sizePolicy)
        self.light_plot_widget.setMinimumSize(QtCore.QSize(500, 300))
        self.light_plot_widget.setObjectName("light_plot_widget")
        self.gridLayout_2.addWidget(self.light_plot_widget, 4, 0, 1, 1)
        self.line_6 = QtWidgets.QFrame(self.tab)
        self.line_6.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.gridLayout_2.addWidget(self.line_6, 4, 1, 1, 1)
        self.light_treewidget_2 = QtWidgets.QTreeWidget(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.light_treewidget_2.sizePolicy().hasHeightForWidth())
        self.light_treewidget_2.setSizePolicy(sizePolicy)
        self.light_treewidget_2.setMinimumSize(QtCore.QSize(256, 0))
        self.light_treewidget_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.light_treewidget_2.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.light_treewidget_2.setDefaultDropAction(QtCore.Qt.IgnoreAction)
        self.light_treewidget_2.setAlternatingRowColors(True)
        self.light_treewidget_2.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.light_treewidget_2.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        self.light_treewidget_2.setIndentation(5)
        self.light_treewidget_2.setItemsExpandable(True)
        self.light_treewidget_2.setAllColumnsShowFocus(False)
        self.light_treewidget_2.setHeaderHidden(False)
        self.light_treewidget_2.setObjectName("light_treewidget_2")
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        item_0 = QtWidgets.QTreeWidgetItem(self.light_treewidget_2)
        self.light_treewidget_2.header().setCascadingSectionResizes(False)
        self.light_treewidget_2.header().setDefaultSectionSize(110)
        self.light_treewidget_2.header().setHighlightSections(False)
        self.light_treewidget_2.header().setMinimumSectionSize(30)
        self.light_treewidget_2.header().setStretchLastSection(True)
        self.gridLayout_2.addWidget(self.light_treewidget_2, 4, 2, 1, 2)
        self.line_2 = QtWidgets.QFrame(self.tab)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout_2.addWidget(self.line_2, 1, 2, 1, 2)
        self.line = QtWidgets.QFrame(self.tab)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout_2.addWidget(self.line, 0, 1, 3, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.tab_2)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.line_4 = QtWidgets.QFrame(self.tab_2)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.gridLayout_6.addWidget(self.line_4, 0, 2, 5, 1)
        self.vel_plot_widget = QtWidgets.QWidget(self.tab_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.vel_plot_widget.sizePolicy().hasHeightForWidth())
        self.vel_plot_widget.setSizePolicy(sizePolicy)
        self.vel_plot_widget.setObjectName("vel_plot_widget")
        self.gridLayout_6.addWidget(self.vel_plot_widget, 0, 3, 5, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout.setObjectName("verticalLayout")
        self.vel_plotmodel_chk = QtWidgets.QCheckBox(self.groupBox_2)
        self.vel_plotmodel_chk.setChecked(True)
        self.vel_plotmodel_chk.setObjectName("vel_plotmodel_chk")
        self.verticalLayout.addWidget(self.vel_plotmodel_chk)
        self.vel_plotobs_chk = QtWidgets.QCheckBox(self.groupBox_2)
        self.vel_plotobs_chk.setChecked(True)
        self.vel_plotobs_chk.setObjectName("vel_plotobs_chk")
        self.verticalLayout.addWidget(self.vel_plotobs_chk)
        self.vel_alias_chk = QtWidgets.QCheckBox(self.groupBox_2)
        self.vel_alias_chk.setObjectName("vel_alias_chk")
        self.verticalLayout.addWidget(self.vel_alias_chk)
        self.gridLayout_6.addWidget(self.groupBox_2, 1, 0, 1, 2)
        self.line_5 = QtWidgets.QFrame(self.tab_2)
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.gridLayout_6.addWidget(self.line_5, 2, 0, 1, 2)
        self.vel_plot_btn = QtWidgets.QPushButton(self.tab_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.vel_plot_btn.sizePolicy().hasHeightForWidth())
        self.vel_plot_btn.setSizePolicy(sizePolicy)
        self.vel_plot_btn.setObjectName("vel_plot_btn")
        self.gridLayout_6.addWidget(self.vel_plot_btn, 3, 0, 1, 2)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_6.addItem(spacerItem, 4, 0, 1, 2)
        self.groupBox = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.vel_pri_label = QtWidgets.QLabel(self.groupBox)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.vel_pri_label.setFont(font)
        self.vel_pri_label.setObjectName("vel_pri_label")
        self.gridLayout_3.addWidget(self.vel_pri_label, 1, 1, 1, 1)
        self.vel_sec_label = QtWidgets.QLabel(self.groupBox)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.vel_sec_label.setFont(font)
        self.vel_sec_label.setObjectName("vel_sec_label")
        self.gridLayout_3.addWidget(self.vel_sec_label, 3, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(15, 20, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem1, 1, 0, 1, 1)
        self.vel_pri_radiobtn = QtWidgets.QRadioButton(self.groupBox)
        self.vel_pri_radiobtn.setObjectName("vel_pri_radiobtn")
        self.gridLayout_3.addWidget(self.vel_pri_radiobtn, 0, 0, 1, 2)
        self.vel_sec_radiobtn = QtWidgets.QRadioButton(self.groupBox)
        self.vel_sec_radiobtn.setObjectName("vel_sec_radiobtn")
        self.gridLayout_3.addWidget(self.vel_sec_radiobtn, 2, 0, 1, 2)
        self.vel_both_radiobtn = QtWidgets.QRadioButton(self.groupBox)
        self.vel_both_radiobtn.setChecked(True)
        self.vel_both_radiobtn.setObjectName("vel_both_radiobtn")
        self.gridLayout_3.addWidget(self.vel_both_radiobtn, 5, 0, 1, 2)
        self.line_7 = QtWidgets.QFrame(self.groupBox)
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.gridLayout_3.addWidget(self.line_7, 4, 0, 1, 2)
        self.gridLayout_5.addLayout(self.gridLayout_3, 0, 0, 1, 1)
        self.gridLayout_6.addWidget(self.groupBox, 0, 0, 1, 2)
        self.tabWidget.addTab(self.tab_2, "")
        self.gridLayout_4.addWidget(self.tabWidget, 0, 0, 1, 1)

        self.retranslateUi(SyntheticCurveWidget)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(SyntheticCurveWidget)

    def retranslateUi(self, SyntheticCurveWidget):
        _translate = QtCore.QCoreApplication.translate
        SyntheticCurveWidget.setWindowTitle(_translate("SyntheticCurveWidget", "Plot Synthetic Curves"))
        self.groupBox_3.setTitle(_translate("SyntheticCurveWidget", "Options"))
        self.light_plotmodel_chk.setToolTip(_translate("SyntheticCurveWidget", "Compute and plot synthetic curve"))
        self.light_plotmodel_chk.setText(_translate("SyntheticCurveWidget", "Plot Model"))
        self.light_plotobs_chk.setToolTip(_translate("SyntheticCurveWidget", "Include observation in plot"))
        self.light_plotobs_chk.setText(_translate("SyntheticCurveWidget", "Plot Observation"))
        self.light_alias_chk.setToolTip(_translate("SyntheticCurveWidget", "<html><head/><body><p>If working with phase and output is set as phase, observations will be folded or cut to make observations fit into model boundaries.</p><p>If working with phase and output is set to JD, observations will be converted to JD and then folded or cut regarding to model boundaries.</p><p>If working with JD and output is set as phase, observations will be converted to phases and then folded or cut regarding to model boundaries.</p><p>If working with HJD and input is in HJD, nothing will be done.</p></body></html>"))
        self.light_alias_chk.setText(_translate("SyntheticCurveWidget", "Alias Observation with Model"))
        self.groupBox_4.setTitle(_translate("SyntheticCurveWidget", "Curves"))
        self.light_treewidget.headerItem().setText(0, _translate("SyntheticCurveWidget", "File Name"))
        self.light_treewidget.headerItem().setText(1, _translate("SyntheticCurveWidget", "Band"))
        self.light_treewidget.headerItem().setText(2, _translate("SyntheticCurveWidget", "L1"))
        self.light_treewidget.headerItem().setText(3, _translate("SyntheticCurveWidget", "L2"))
        self.light_treewidget.headerItem().setText(4, _translate("SyntheticCurveWidget", "L3"))
        self.light_treewidget.headerItem().setText(5, _translate("SyntheticCurveWidget", "X1"))
        self.light_treewidget.headerItem().setText(6, _translate("SyntheticCurveWidget", "X2"))
        self.light_treewidget.headerItem().setText(7, _translate("SyntheticCurveWidget", "Y1"))
        self.light_treewidget.headerItem().setText(8, _translate("SyntheticCurveWidget", "Y2"))
        self.light_treewidget.headerItem().setText(9, _translate("SyntheticCurveWidget", "Opacity"))
        self.light_treewidget.headerItem().setText(10, _translate("SyntheticCurveWidget", "Extinction"))
        self.light_treewidget.headerItem().setText(11, _translate("SyntheticCurveWidget", "Calibration"))
        self.light_treewidget.headerItem().setText(12, _translate("SyntheticCurveWidget", "Factor"))
        self.light_treewidget.headerItem().setText(13, _translate("SyntheticCurveWidget", "Zero"))
        self.light_plot_btn.setText(_translate("SyntheticCurveWidget", "Plot"))
        self.light_treewidget_2.headerItem().setText(0, _translate("SyntheticCurveWidget", "Parameter"))
        self.light_treewidget_2.headerItem().setText(1, _translate("SyntheticCurveWidget", "Value"))
        __sortingEnabled = self.light_treewidget_2.isSortingEnabled()
        self.light_treewidget_2.setSortingEnabled(False)
        self.light_treewidget_2.topLevelItem(0).setText(0, _translate("SyntheticCurveWidget", "M1 (M☉)"))
        self.light_treewidget_2.topLevelItem(1).setText(0, _translate("SyntheticCurveWidget", "M2 (M☉)"))
        self.light_treewidget_2.topLevelItem(2).setText(0, _translate("SyntheticCurveWidget", "R1 (R☉)"))
        self.light_treewidget_2.topLevelItem(3).setText(0, _translate("SyntheticCurveWidget", "R2 (R☉)"))
        self.light_treewidget_2.topLevelItem(4).setText(0, _translate("SyntheticCurveWidget", "Teff1 (K)"))
        self.light_treewidget_2.topLevelItem(5).setText(0, _translate("SyntheticCurveWidget", "Teff2 (K)"))
        self.light_treewidget_2.topLevelItem(6).setText(0, _translate("SyntheticCurveWidget", "logg1 (cgs)"))
        self.light_treewidget_2.topLevelItem(7).setText(0, _translate("SyntheticCurveWidget", "logg2 (cgs)"))
        self.light_treewidget_2.topLevelItem(8).setText(0, _translate("SyntheticCurveWidget", "Mbol1 (mag)"))
        self.light_treewidget_2.topLevelItem(9).setText(0, _translate("SyntheticCurveWidget", "Mbol2 (mag)"))
        self.light_treewidget_2.topLevelItem(10).setText(0, _translate("SyntheticCurveWidget", "log (L1/L☉)"))
        self.light_treewidget_2.topLevelItem(11).setText(0, _translate("SyntheticCurveWidget", "log (L2/L☉)"))
        self.light_treewidget_2.setSortingEnabled(__sortingEnabled)
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("SyntheticCurveWidget", "Light Curve"))
        self.groupBox_2.setTitle(_translate("SyntheticCurveWidget", "Options"))
        self.vel_plotmodel_chk.setToolTip(_translate("SyntheticCurveWidget", "Compute and plot synthetic curve"))
        self.vel_plotmodel_chk.setText(_translate("SyntheticCurveWidget", "Plot Model"))
        self.vel_plotobs_chk.setToolTip(_translate("SyntheticCurveWidget", "Include observation in plot"))
        self.vel_plotobs_chk.setText(_translate("SyntheticCurveWidget", "Plot Observation"))
        self.vel_alias_chk.setToolTip(_translate("SyntheticCurveWidget", "<html><head/><body><p>If working with phase and output is set as phase, observations will be folded or cut to make observations fit into model boundaries.</p><p>If working with phase and output is set to JD, observations will be converted to JD and then folded or cut regarding to model boundaries.</p><p>If working with JD and output is set as phase, observations will be converted to phases and then folded or cut regarding to model boundaries.</p><p>If working with HJD and input is in HJD, nothing will be done.</p></body></html>"))
        self.vel_alias_chk.setText(_translate("SyntheticCurveWidget", "Alias Observation with Model"))
        self.vel_plot_btn.setToolTip(_translate("SyntheticCurveWidget", "Plot selected item"))
        self.vel_plot_btn.setText(_translate("SyntheticCurveWidget", "Plot"))
        self.groupBox.setTitle(_translate("SyntheticCurveWidget", "Curves"))
        self.vel_pri_label.setText(_translate("SyntheticCurveWidget", "[Synthetic]"))
        self.vel_sec_label.setText(_translate("SyntheticCurveWidget", "[Synthetic]"))
        self.vel_pri_radiobtn.setText(_translate("SyntheticCurveWidget", "Primar&y Curve"))
        self.vel_sec_radiobtn.setText(_translate("SyntheticCurveWidget", "Secondary C&urve"))
        self.vel_both_radiobtn.setText(_translate("SyntheticCurveWidget", "Both"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("SyntheticCurveWidget", "Velocity Curve"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    SyntheticCurveWidget = QtWidgets.QWidget()
    ui = Ui_SyntheticCurveWidget()
    ui.setupUi(SyntheticCurveWidget)
    SyntheticCurveWidget.show()
    sys.exit(app.exec_())

