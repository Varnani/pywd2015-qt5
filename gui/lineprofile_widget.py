# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'lineprofile_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_LineProfileWidget(object):
    def setupUi(self, LineProfileWidget):
        LineProfileWidget.setObjectName("LineProfileWidget")
        LineProfileWidget.resize(1100, 650)
        self.gridLayout_4 = QtWidgets.QGridLayout(LineProfileWidget)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.tabWidget = QtWidgets.QTabWidget(LineProfileWidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.s1_groupbox = QtWidgets.QGroupBox(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s1_groupbox.sizePolicy().hasHeightForWidth())
        self.s1_groupbox.setSizePolicy(sizePolicy)
        self.s1_groupbox.setMinimumSize(QtCore.QSize(500, 0))
        self.s1_groupbox.setMaximumSize(QtCore.QSize(598, 16777215))
        self.s1_groupbox.setCheckable(True)
        self.s1_groupbox.setObjectName("s1_groupbox")
        self.gridLayout = QtWidgets.QGridLayout(self.s1_groupbox)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self.s1_groupbox)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.s1_binwidth_spinbox = QtWidgets.QDoubleSpinBox(self.s1_groupbox)
        self.s1_binwidth_spinbox.setDecimals(6)
        self.s1_binwidth_spinbox.setMinimum(1e-06)
        self.s1_binwidth_spinbox.setMaximum(1.0)
        self.s1_binwidth_spinbox.setSingleStep(0.0001)
        self.s1_binwidth_spinbox.setProperty("value", 1e-05)
        self.s1_binwidth_spinbox.setObjectName("s1_binwidth_spinbox")
        self.gridLayout.addWidget(self.s1_binwidth_spinbox, 0, 1, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.s1_groupbox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 0, 2, 1, 1)
        self.s1_contscale_spinbox = QtWidgets.QDoubleSpinBox(self.s1_groupbox)
        self.s1_contscale_spinbox.setDecimals(5)
        self.s1_contscale_spinbox.setMinimum(0.0)
        self.s1_contscale_spinbox.setMaximum(99.0)
        self.s1_contscale_spinbox.setProperty("value", 1.0)
        self.s1_contscale_spinbox.setObjectName("s1_contscale_spinbox")
        self.gridLayout.addWidget(self.s1_contscale_spinbox, 0, 3, 1, 1)
        self.s1_add_btn = QtWidgets.QPushButton(self.s1_groupbox)
        self.s1_add_btn.setObjectName("s1_add_btn")
        self.gridLayout.addWidget(self.s1_add_btn, 0, 4, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.s1_groupbox)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 1, 0, 1, 1)
        self.s1_subgrid_spinbox = QtWidgets.QSpinBox(self.s1_groupbox)
        self.s1_subgrid_spinbox.setMinimum(1)
        self.s1_subgrid_spinbox.setMaximum(999)
        self.s1_subgrid_spinbox.setProperty("value", 1)
        self.s1_subgrid_spinbox.setObjectName("s1_subgrid_spinbox")
        self.gridLayout.addWidget(self.s1_subgrid_spinbox, 1, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.s1_groupbox)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 1, 2, 1, 1)
        self.s1_contslope_spinbox = QtWidgets.QDoubleSpinBox(self.s1_groupbox)
        self.s1_contslope_spinbox.setDecimals(5)
        self.s1_contslope_spinbox.setObjectName("s1_contslope_spinbox")
        self.gridLayout.addWidget(self.s1_contslope_spinbox, 1, 3, 1, 1)
        self.s1_remove_btn = QtWidgets.QPushButton(self.s1_groupbox)
        self.s1_remove_btn.setObjectName("s1_remove_btn")
        self.gridLayout.addWidget(self.s1_remove_btn, 1, 4, 1, 1)
        self.s1_treewidget = QtWidgets.QTreeWidget(self.s1_groupbox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s1_treewidget.sizePolicy().hasHeightForWidth())
        self.s1_treewidget.setSizePolicy(sizePolicy)
        self.s1_treewidget.setMinimumSize(QtCore.QSize(0, 150))
        self.s1_treewidget.setMaximumSize(QtCore.QSize(16777215, 270))
        self.s1_treewidget.setObjectName("s1_treewidget")
        self.gridLayout.addWidget(self.s1_treewidget, 2, 0, 1, 5)
        self.gridLayout_3.addWidget(self.s1_groupbox, 0, 0, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.plot_btn = QtWidgets.QPushButton(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plot_btn.sizePolicy().hasHeightForWidth())
        self.plot_btn.setSizePolicy(sizePolicy)
        self.plot_btn.setMinimumSize(QtCore.QSize(125, 0))
        self.plot_btn.setObjectName("plot_btn")
        self.horizontalLayout.addWidget(self.plot_btn)
        self.line = QtWidgets.QFrame(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.line.sizePolicy().hasHeightForWidth())
        self.line.setSizePolicy(sizePolicy)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.horizontalLayout.addWidget(self.line)
        self.label_3 = QtWidgets.QLabel(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout.addWidget(self.label_3)
        self.phase_spinbox = QtWidgets.QDoubleSpinBox(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.phase_spinbox.sizePolicy().hasHeightForWidth())
        self.phase_spinbox.setSizePolicy(sizePolicy)
        self.phase_spinbox.setDecimals(4)
        self.phase_spinbox.setSingleStep(0.1)
        self.phase_spinbox.setProperty("value", 0.25)
        self.phase_spinbox.setObjectName("phase_spinbox")
        self.horizontalLayout.addWidget(self.phase_spinbox)
        self.gridLayout_3.addLayout(self.horizontalLayout, 3, 2, 1, 1)
        self.plot_widget = QtWidgets.QWidget(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.plot_widget.sizePolicy().hasHeightForWidth())
        self.plot_widget.setSizePolicy(sizePolicy)
        self.plot_widget.setObjectName("plot_widget")
        self.gridLayout_3.addWidget(self.plot_widget, 0, 2, 2, 1)
        self.s2_groupbox = QtWidgets.QGroupBox(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s2_groupbox.sizePolicy().hasHeightForWidth())
        self.s2_groupbox.setSizePolicy(sizePolicy)
        self.s2_groupbox.setMinimumSize(QtCore.QSize(500, 0))
        self.s2_groupbox.setMaximumSize(QtCore.QSize(598, 16777215))
        self.s2_groupbox.setCheckable(True)
        self.s2_groupbox.setObjectName("s2_groupbox")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.s2_groupbox)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_15 = QtWidgets.QLabel(self.s2_groupbox)
        self.label_15.setObjectName("label_15")
        self.gridLayout_2.addWidget(self.label_15, 0, 0, 1, 1)
        self.s2_binwidth_spinbox = QtWidgets.QDoubleSpinBox(self.s2_groupbox)
        self.s2_binwidth_spinbox.setDecimals(6)
        self.s2_binwidth_spinbox.setMinimum(1e-06)
        self.s2_binwidth_spinbox.setMaximum(1.0)
        self.s2_binwidth_spinbox.setSingleStep(0.0001)
        self.s2_binwidth_spinbox.setProperty("value", 1e-05)
        self.s2_binwidth_spinbox.setObjectName("s2_binwidth_spinbox")
        self.gridLayout_2.addWidget(self.s2_binwidth_spinbox, 0, 1, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.s2_groupbox)
        self.label_12.setObjectName("label_12")
        self.gridLayout_2.addWidget(self.label_12, 0, 2, 1, 1)
        self.s2_contscale_spinbox = QtWidgets.QDoubleSpinBox(self.s2_groupbox)
        self.s2_contscale_spinbox.setDecimals(5)
        self.s2_contscale_spinbox.setMaximum(99.0)
        self.s2_contscale_spinbox.setProperty("value", 1.0)
        self.s2_contscale_spinbox.setObjectName("s2_contscale_spinbox")
        self.gridLayout_2.addWidget(self.s2_contscale_spinbox, 0, 3, 1, 1)
        self.s2_add_btn = QtWidgets.QPushButton(self.s2_groupbox)
        self.s2_add_btn.setObjectName("s2_add_btn")
        self.gridLayout_2.addWidget(self.s2_add_btn, 0, 4, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.s2_groupbox)
        self.label_13.setObjectName("label_13")
        self.gridLayout_2.addWidget(self.label_13, 1, 0, 1, 1)
        self.s2_subgrid_spinbox = QtWidgets.QSpinBox(self.s2_groupbox)
        self.s2_subgrid_spinbox.setMinimum(1)
        self.s2_subgrid_spinbox.setObjectName("s2_subgrid_spinbox")
        self.gridLayout_2.addWidget(self.s2_subgrid_spinbox, 1, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.s2_groupbox)
        self.label_11.setObjectName("label_11")
        self.gridLayout_2.addWidget(self.label_11, 1, 2, 1, 1)
        self.s2_contslope_spinbox = QtWidgets.QDoubleSpinBox(self.s2_groupbox)
        self.s2_contslope_spinbox.setDecimals(5)
        self.s2_contslope_spinbox.setMaximum(1.0)
        self.s2_contslope_spinbox.setObjectName("s2_contslope_spinbox")
        self.gridLayout_2.addWidget(self.s2_contslope_spinbox, 1, 3, 1, 1)
        self.s2_remove_btn = QtWidgets.QPushButton(self.s2_groupbox)
        self.s2_remove_btn.setObjectName("s2_remove_btn")
        self.gridLayout_2.addWidget(self.s2_remove_btn, 1, 4, 1, 1)
        self.s2_treewidget = QtWidgets.QTreeWidget(self.s2_groupbox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.s2_treewidget.sizePolicy().hasHeightForWidth())
        self.s2_treewidget.setSizePolicy(sizePolicy)
        self.s2_treewidget.setMinimumSize(QtCore.QSize(0, 150))
        self.s2_treewidget.setMaximumSize(QtCore.QSize(16777215, 270))
        self.s2_treewidget.setObjectName("s2_treewidget")
        self.gridLayout_2.addWidget(self.s2_treewidget, 2, 0, 1, 5)
        self.gridLayout_3.addWidget(self.s2_groupbox, 1, 0, 3, 1)
        self.line_2 = QtWidgets.QFrame(self.tab)
        self.line_2.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout_3.addWidget(self.line_2, 0, 1, 4, 1)
        self.line_3 = QtWidgets.QFrame(self.tab)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.gridLayout_3.addWidget(self.line_3, 2, 2, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
        self.gridLayout_4.addWidget(self.tabWidget, 0, 0, 1, 1)

        self.retranslateUi(LineProfileWidget)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(LineProfileWidget)

    def retranslateUi(self, LineProfileWidget):
        _translate = QtCore.QCoreApplication.translate
        LineProfileWidget.setWindowTitle(_translate("LineProfileWidget", "Spectral Line Profiles"))
        self.s1_groupbox.setTitle(_translate("LineProfileWidget", "Star 1"))
        self.label.setToolTip(_translate("LineProfileWidget", "The bin width in microns. Too small a bin width gives noisy profiles. Too large a bin\n"
"width gives insufficient spectral resolution."))
        self.label.setText(_translate("LineProfileWidget", "Bin Width"))
        self.label_4.setToolTip(_translate("LineProfileWidget", "The continuum scale (continuum flux at the reference wavelength). The unit is decided by\n"
"the user."))
        self.label_4.setText(_translate("LineProfileWidget", "Cont. Scale"))
        self.s1_add_btn.setText(_translate("LineProfileWidget", "Add"))
        self.label_6.setToolTip(_translate("LineProfileWidget", "<html><head/><body><p>Grid fineness for micro-integration on each surface element. Subgrid = 1 means that there is no micro-integration. Subgrid = n breaks each surface element into n<span style=\" vertical-align:super;\">2</span> pieces, each with its own radial velocity, thus improving integration accuracy.</p></body></html>"))
        self.label_6.setText(_translate("LineProfileWidget", "Subgrid"))
        self.label_5.setToolTip(_translate("LineProfileWidget", "The continuum slope in flux units per micron."))
        self.label_5.setText(_translate("LineProfileWidget", "Cont. Slope"))
        self.s1_remove_btn.setText(_translate("LineProfileWidget", "Remove"))
        self.s1_treewidget.headerItem().setText(0, _translate("LineProfileWidget", "Wavelength"))
        self.s1_treewidget.headerItem().setToolTip(0, _translate("LineProfileWidget", "The line rest wavelength in microns"))
        self.s1_treewidget.headerItem().setText(1, _translate("LineProfileWidget", "Equivalent Width"))
        self.s1_treewidget.headerItem().setToolTip(1, _translate("LineProfileWidget", "<html><head/><body><p>The line equivalent width for a line, in microns—the traditional measure of line strength.</p><p>Absorption and emission lines both have positive equivalent width by program (WD) convention. Whether a line is in absorption or emission is controlled by parameters Rectangle Line Depth.</p></body></html>"))
        self.s1_treewidget.headerItem().setText(2, _translate("LineProfileWidget", "Rect. Line Depth"))
        self.s1_treewidget.headerItem().setToolTip(2, _translate("LineProfileWidget", "<html><head/><body><p>Rectangular line depth for a line. </p><p>Line profiles are formed by binning of Doppler shifted elements that have rectangular profiles, each with a depth and a width. The user supplies the depth and the program (WD) then calculates the rectangular line width needed to reproduce the specified equivalent width. The depth is relative to a unit continuum, so 0.80000 means that 80 percent of the continuum flux is missing within the rectangular profile element, or that the residual flux is 20 percent of the continuum. Negative depths correspond to emission lines, so −0.50000 means 50 percent above the continuum. Depths must be less than 1.0000 (i.e. an absorption line cannot go to zero flux or below), but can be less than −1.0000 (an emission line can go arbitrarily high).</p></body></html>"))
        self.s1_treewidget.headerItem().setText(3, _translate("LineProfileWidget", "KKS"))
        self.s1_treewidget.headerItem().setToolTip(3, _translate("LineProfileWidget", "<html><head/><body><p>This integer specifies a surface region associated with a given spectral line. </p><p>If KKS=0, the line is not specific to a location but applies to the entire star. </p><p>If KKS=1, then the line applies only to the first spot on that star</p><p>if KKS=2 it applies only to the second spot, and so on. </p><p>Naturally the star must have spots for this scheme to work, but the spots need not be hot or cool spots—they can have temperature factors of unity. </p><p>Negative KKS specifies avoidance of regions. Thus KKS=−4 means that the spectral line applies everywhere on the star except within spot 4. If you find this confusing, just set KKS=0 and the line applies in the old simple way—everywhere on the star.</p></body></html>"))
        self.plot_btn.setText(_translate("LineProfileWidget", "Plot"))
        self.label_3.setText(_translate("LineProfileWidget", "Phase"))
        self.s2_groupbox.setTitle(_translate("LineProfileWidget", "Star 2"))
        self.label_15.setToolTip(_translate("LineProfileWidget", "The bin width in microns. Too small a bin width gives noisy profiles. Too large a bin\n"
"width gives insufficient spectral resolution."))
        self.label_15.setText(_translate("LineProfileWidget", "Bin Width"))
        self.label_12.setToolTip(_translate("LineProfileWidget", "The continuum scale (continuum flux at the reference wavelength). The unit is decided by\n"
"the user."))
        self.label_12.setText(_translate("LineProfileWidget", "Cont. Scale"))
        self.s2_add_btn.setText(_translate("LineProfileWidget", "Add"))
        self.label_13.setToolTip(_translate("LineProfileWidget", "<html><head/><body><p>Grid fineness for micro-integration on each surface element. Subgrid = 1 means that there is no micro-integration. Subgrid = n breaks each surface element into n<span style=\" vertical-align:super;\">2</span> pieces, each with its own radial velocity, thus improving integration accuracy.</p></body></html>"))
        self.label_13.setText(_translate("LineProfileWidget", "Subgrid"))
        self.label_11.setToolTip(_translate("LineProfileWidget", "The continuum slope in flux units per micron."))
        self.label_11.setText(_translate("LineProfileWidget", "Cont. Slope"))
        self.s2_remove_btn.setText(_translate("LineProfileWidget", "Remove"))
        self.s2_treewidget.headerItem().setText(0, _translate("LineProfileWidget", "Wavelength"))
        self.s2_treewidget.headerItem().setToolTip(0, _translate("LineProfileWidget", "The line rest wavelength in microns"))
        self.s2_treewidget.headerItem().setText(1, _translate("LineProfileWidget", "Equivalent Width"))
        self.s2_treewidget.headerItem().setToolTip(1, _translate("LineProfileWidget", "<html><head/><body><p>The line equivalent width for a line, in microns—the traditional measure of line strength.</p><p>Absorption and emission lines both have positive equivalent width by program (WD) convention. Whether a line is in absorption or emission is controlled by parameters Rectangle Line Depth.</p></body></html>"))
        self.s2_treewidget.headerItem().setText(2, _translate("LineProfileWidget", "Rect. Line Depth"))
        self.s2_treewidget.headerItem().setToolTip(2, _translate("LineProfileWidget", "<html><head/><body><p>Rectangular line depth for a line. </p><p>Line profiles are formed by binning of Doppler shifted elements that have rectangular profiles, each with a depth and a width. The user supplies the depth and the program (WD) then calculates the rectangular line width needed to reproduce the specified equivalent width. The depth is relative to a unit continuum, so 0.80000 means that 80 percent of the continuum flux is missing within the rectangular profile element, or that the residual flux is 20 percent of the continuum. Negative depths correspond to emission lines, so −0.50000 means 50 percent above the continuum. Depths must be less than 1.0000 (i.e. an absorption line cannot go to zero flux or below), but can be less than −1.0000 (an emission line can go arbitrarily high).</p></body></html>"))
        self.s2_treewidget.headerItem().setText(3, _translate("LineProfileWidget", "KKS"))
        self.s2_treewidget.headerItem().setToolTip(3, _translate("LineProfileWidget", "<html><head/><body><p>This integer specifies a surface region associated with a given spectral line. </p><p>If KKS=0, the line is not specific to a location but applies to the entire star. </p><p>If KKS=1, then the line applies only to the first spot on that star</p><p>if KKS=2 it applies only to the second spot, and so on. </p><p>Naturally the star must have spots for this scheme to work, but the spots need not be hot or cool spots—they can have temperature factors of unity. </p><p>Negative KKS specifies avoidance of regions. Thus KKS=−4 means that the spectral line applies everywhere on the star except within spot 4. If you find this confusing, just set KKS=0 and the line applies in the old simple way—everywhere on the star.</p></body></html>"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("LineProfileWidget", "Single"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("LineProfileWidget", "Animation"))

