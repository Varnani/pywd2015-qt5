# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'outputview_widget.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_OutputView(object):
    def setupUi(self, OutputView):
        OutputView.setObjectName("OutputView")
        OutputView.resize(900, 500)
        OutputView.setMinimumSize(QtCore.QSize(500, 250))
        self.gridLayout = QtWidgets.QGridLayout(OutputView)
        self.gridLayout.setObjectName("gridLayout")
        self.jmp_beginning_btn = QtWidgets.QPushButton(OutputView)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jmp_beginning_btn.sizePolicy().hasHeightForWidth())
        self.jmp_beginning_btn.setSizePolicy(sizePolicy)
        self.jmp_beginning_btn.setObjectName("jmp_beginning_btn")
        self.gridLayout.addWidget(self.jmp_beginning_btn, 0, 0, 1, 1)
        self.jmp_end_btn = QtWidgets.QPushButton(OutputView)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jmp_end_btn.sizePolicy().hasHeightForWidth())
        self.jmp_end_btn.setSizePolicy(sizePolicy)
        self.jmp_end_btn.setObjectName("jmp_end_btn")
        self.gridLayout.addWidget(self.jmp_end_btn, 0, 1, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.jmp_text_btn = QtWidgets.QPushButton(OutputView)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jmp_text_btn.sizePolicy().hasHeightForWidth())
        self.jmp_text_btn.setSizePolicy(sizePolicy)
        self.jmp_text_btn.setObjectName("jmp_text_btn")
        self.horizontalLayout.addWidget(self.jmp_text_btn)
        self.jmp_text_tipt = QtWidgets.QLineEdit(OutputView)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jmp_text_tipt.sizePolicy().hasHeightForWidth())
        self.jmp_text_tipt.setSizePolicy(sizePolicy)
        self.jmp_text_tipt.setObjectName("jmp_text_tipt")
        self.horizontalLayout.addWidget(self.jmp_text_tipt)
        self.gridLayout.addLayout(self.horizontalLayout, 1, 0, 1, 2)
        self.output_textedit = QtWidgets.QPlainTextEdit(OutputView)
        self.output_textedit.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.output_textedit.setReadOnly(True)
        self.output_textedit.setObjectName("output_textedit")
        self.gridLayout.addWidget(self.output_textedit, 2, 0, 1, 2)

        self.retranslateUi(OutputView)
        QtCore.QMetaObject.connectSlotsByName(OutputView)

    def retranslateUi(self, OutputView):
        _translate = QtCore.QCoreApplication.translate
        OutputView.setWindowTitle(_translate("OutputView", "///"))
        self.jmp_beginning_btn.setText(_translate("OutputView", "Jump to Beginning"))
        self.jmp_end_btn.setText(_translate("OutputView", "Jump to End"))
        self.jmp_text_btn.setText(_translate("OutputView", "Jump to Line Containing Text:"))

