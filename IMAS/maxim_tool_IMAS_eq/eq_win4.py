# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'eq_win3.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form_eq(object):
    def buildUI(self, Form_eq):
        Form_eq.setObjectName("Form_eq")
        Form_eq.resize(1080, 720)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form_eq.sizePolicy().hasHeightForWidth())
        Form_eq.setSizePolicy(sizePolicy)
        Form_eq.setMinimumSize(QtCore.QSize(0, 508))
        self.gridLayout = QtWidgets.QGridLayout(Form_eq)
        self.gridLayout.setObjectName("gridLayout")
        self.graplay = QtWidgets.QGridLayout()
        self.graplay.setObjectName("graplay")
        self.gridLayout.addLayout(self.graplay, 0, 0, 1, 1)

        self.retranslateUi(Form_eq)
        QtCore.QMetaObject.connectSlotsByName(Form_eq)

    def retranslateUi(self, Form_eq):
        _translate = QtCore.QCoreApplication.translate
        Form_eq.setWindowTitle(_translate("Form_eq", "EQUILIBRIUM-DINA"))

