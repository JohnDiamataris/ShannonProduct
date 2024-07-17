# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainPaper2extened.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from QtPhenmodelmaker import Phenomenon
from gen_to_phen import phen_to_gen
import matplotlib.pyplot as plt
import itertools as it 
import scipy
from scipy.stats import binom
import sys
import os
import random
import csv
import numpy as np
from scipy.optimize import minimize
from scipy import stats
from PyQt5.QtGui import QPixmap


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(859, 555)
        MainWindow.setStyleSheet("")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.GenoTable = QtWidgets.QTableWidget(self.centralwidget)
        self.GenoTable.setGeometry(QtCore.QRect(20, 178, 161, 291))
        self.GenoTable.setObjectName("GenoTable")
        self.GenoTable.setColumnCount(0)
        self.GenoTable.setRowCount(0)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(20, 150, 151, 21))
        self.label.setObjectName("label")
        self.PhenoButton = QtWidgets.QPushButton(self.centralwidget)
        self.PhenoButton.setGeometry(QtCore.QRect(141, 98, 91, 23))
        self.PhenoButton.setObjectName("PhenoButton")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(26, 79, 151, 21))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(28, 9, 311, 51))
        self.label_3.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_3.setAlignment(QtCore.Qt.AlignCenter)
        self.label_3.setWordWrap(True)
        self.label_3.setObjectName("label_3")
        self.GeneInputTextEdit = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.GeneInputTextEdit.setGeometry(QtCore.QRect(23, 50, 161, 31))
        self.GeneInputTextEdit.setObjectName("GeneInputTextEdit")
        self.GenoTable_2 = QtWidgets.QTableWidget(self.centralwidget)
        self.GenoTable_2.setGeometry(QtCore.QRect(200, 178, 161, 291))
        self.GenoTable_2.setObjectName("GenoTable_2")
        self.GenoTable_2.setColumnCount(0)
        self.GenoTable_2.setRowCount(0)
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(201, 150, 141, 21))
        self.label_4.setObjectName("label_4")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(60, 1, 101, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        font.setKerning(True)
        self.label_6.setFont(font)
        self.label_6.setStyleSheet("color: rgb(55, 58, 255);")
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(240, 1, 101, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        font.setKerning(True)
        self.label_7.setFont(font)
        self.label_7.setStyleSheet("color: rgb(55, 58, 255);")
        self.label_7.setObjectName("label_7")
        self.copyParent = QtWidgets.QPushButton(self.centralwidget)
        self.copyParent.setGeometry(QtCore.QRect(372, 103, 131, 23))
        self.copyParent.setObjectName("copyParent")
        self.SampleSizeSpinBox = QtWidgets.QSpinBox(self.centralwidget)
        self.SampleSizeSpinBox.setGeometry(QtCore.QRect(434, 73, 60, 22))
        self.SampleSizeSpinBox.setSpecialValueText("")
        self.SampleSizeSpinBox.setMaximum(99999)
        self.SampleSizeSpinBox.setSingleStep(1)
        self.SampleSizeSpinBox.setProperty("value", 1)
        self.SampleSizeSpinBox.setObjectName("SampleSizeSpinBox")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(373, 20, 131, 21))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(False)
        font.setUnderline(True)
        font.setWeight(50)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.line_2 = QtWidgets.QFrame(self.centralwidget)
        self.line_2.setGeometry(QtCore.QRect(370, -5, 481, 20))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(30, 470, 311, 31))
        self.label_12.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_12.setAlignment(QtCore.Qt.AlignCenter)
        self.label_12.setWordWrap(True)
        self.label_12.setObjectName("label_12")
        self.label_14 = QtWidgets.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(370, 67, 61, 31))
        self.label_14.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_14.setAlignment(QtCore.Qt.AlignCenter)
        self.label_14.setWordWrap(True)
        self.label_14.setObjectName("label_14")
        self.GenoTable_5 = QtWidgets.QTableWidget(self.centralwidget)
        self.GenoTable_5.setGeometry(QtCore.QRect(705, 30, 130, 161))
        self.GenoTable_5.setObjectName("GenoTable_5")
        self.GenoTable_5.setColumnCount(0)
        self.GenoTable_5.setRowCount(0)
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(705, 12, 130, 21))
        self.label_13.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_13.setAlignment(QtCore.Qt.AlignCenter)
        self.label_13.setWordWrap(True)
        self.label_13.setObjectName("label_13")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(720, 0, 100, 21))
        self.label_15.setObjectName("label_15")
        self.SampleSizeSpinBox_2 = QtWidgets.QSpinBox(self.centralwidget)
        self.SampleSizeSpinBox_2.setGeometry(QtCore.QRect(434, 44, 60, 22))
        self.SampleSizeSpinBox_2.setSpecialValueText("")
        self.SampleSizeSpinBox_2.setMaximum(99999)
        self.SampleSizeSpinBox_2.setSingleStep(1)
        self.SampleSizeSpinBox_2.setProperty("value", 500)
        self.SampleSizeSpinBox_2.setObjectName("SampleSizeSpinBox_2")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(370, 37, 61, 31))
        self.label_16.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_16.setAlignment(QtCore.Qt.AlignCenter)
        self.label_16.setWordWrap(True)
        self.label_16.setObjectName("label_16")
        self.line_3 = QtWidgets.QFrame(self.centralwidget)
        self.line_3.setGeometry(QtCore.QRect(370, 91, 131, 20))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.tableWidget = QtWidgets.QTableWidget(self.centralwidget)
        self.tableWidget.setGeometry(QtCore.QRect(516, 240, 321, 251))
        self.tableWidget.setStyleSheet("color: rgb(0, 0, 255);")
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(520, 218, 81, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.label_10.setFont(font)
        self.label_10.setStyleSheet("color: rgb(0, 0, 255);")
        self.label_10.setObjectName("label_10")
        self.copyParent_2 = QtWidgets.QPushButton(self.centralwidget)
        self.copyParent_2.setGeometry(QtCore.QRect(372, 173, 131, 23))
        self.copyParent_2.setObjectName("copyParent_2")
        self.GeneInputTextEdit_2 = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.GeneInputTextEdit_2.setGeometry(QtCore.QRect(200, 50, 161, 31))
        self.GeneInputTextEdit_2.setObjectName("GeneInputTextEdit_2")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(202, 78, 151, 21))
        self.label_5.setObjectName("label_5")
        self.line_4 = QtWidgets.QFrame(self.centralwidget)
        self.line_4.setGeometry(QtCore.QRect(26, 500, 831, 20))
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.GenoTable_6 = QtWidgets.QTableWidget(self.centralwidget)
        self.GenoTable_6.setGeometry(QtCore.QRect(570, 30, 130, 161))
        self.GenoTable_6.setObjectName("GenoTable_6")
        self.GenoTable_6.setColumnCount(0)
        self.GenoTable_6.setRowCount(0)
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(588, 1, 100, 21))
        self.label_18.setObjectName("label_18")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(570, 12, 130, 21))
        self.label_19.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_19.setAlignment(QtCore.Qt.AlignCenter)
        self.label_19.setWordWrap(True)
        self.label_19.setObjectName("label_19")
        self.Fit = QtWidgets.QPushButton(self.centralwidget)
        self.Fit.setGeometry(QtCore.QRect(372, 126, 131, 23))
        self.Fit.setObjectName("Fit")
        self.copyParent_5 = QtWidgets.QPushButton(self.centralwidget)
        self.copyParent_5.setGeometry(QtCore.QRect(372, 196, 131, 23))
        self.copyParent_5.setObjectName("copyParent_5")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(30, 124, 311, 31))
        self.label_8.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_8.setAlignment(QtCore.Qt.AlignCenter)
        self.label_8.setWordWrap(True)
        self.label_8.setObjectName("label_8")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(5, 0, 361, 511))
        self.groupBox.setStyleSheet("background-color: rgb(255, 255, 220);")
        self.groupBox.setTitle("")
        self.groupBox.setObjectName("groupBox")
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setGeometry(QtCore.QRect(364, 0, 491, 221))
        self.groupBox_2.setStyleSheet("background-color: rgb(212, 215, 255);")
        self.groupBox_2.setTitle("")
        self.groupBox_2.setObjectName("groupBox_2")
        self.label_17 = QtWidgets.QLabel(self.groupBox_2)
        self.label_17.setGeometry(QtCore.QRect(341, 190, 131, 31))
        self.label_17.setStyleSheet("color: rgb(98, 98, 98);")
        self.label_17.setAlignment(QtCore.Qt.AlignCenter)
        self.label_17.setWordWrap(True)
        self.label_17.setObjectName("label_17")
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setGeometry(QtCore.QRect(365, 220, 491, 291))
        self.groupBox_3.setStyleSheet("background-color: rgb(255, 171, 165);")
        self.groupBox_3.setTitle("")
        self.groupBox_3.setObjectName("groupBox_3")
        self.label_11 = QtWidgets.QLabel(self.groupBox_3)
        self.label_11.setGeometry(QtCore.QRect(10, 180, 131, 81))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setScaledContents(False)
        self.label_11.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.label_11.setWordWrap(True)
        self.label_11.setObjectName("label_11")
        self.copyParent_6 = QtWidgets.QPushButton(self.centralwidget)
        self.copyParent_6.setGeometry(QtCore.QRect(372, 149, 131, 23))
        self.copyParent_6.setObjectName("copyParent_6")
        self.groupBox_2.raise_()
        self.groupBox_3.raise_()
        self.groupBox.raise_()
        self.label_3.raise_()
        self.GenoTable.raise_()
        self.label.raise_()
        self.label_2.raise_()
        self.GenoTable_2.raise_()
        self.label_4.raise_()
        self.label_6.raise_()
        self.label_7.raise_()
        self.line_2.raise_()
        self.label_9.raise_()
        self.label_12.raise_()
        self.label_14.raise_()
        self.label_13.raise_()
        self.label_15.raise_()
        self.label_16.raise_()
        self.line_3.raise_()
        self.label_10.raise_()
        self.label_5.raise_()
        self.line_4.raise_()
        self.label_18.raise_()
        self.label_19.raise_()
        self.label_8.raise_()
        self.PhenoButton.raise_()
        self.copyParent_5.raise_()
        self.copyParent.raise_()
        self.copyParent_2.raise_()
        self.Fit.raise_()
        self.SampleSizeSpinBox_2.raise_()
        self.SampleSizeSpinBox.raise_()
        self.GenoTable_5.raise_()
        self.GenoTable_6.raise_()
        self.tableWidget.raise_()
        self.GeneInputTextEdit_2.raise_()
        self.GeneInputTextEdit.raise_()
        self.copyParent_6.raise_()
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 859, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionInput_Sample = QtWidgets.QAction(MainWindow)
        self.actionInput_Sample.setObjectName("actionInput_Sample")
        self.actionMake_Phenotype = QtWidgets.QAction(MainWindow)
        self.actionMake_Phenotype.setObjectName("actionMake_Phenotype")
        self.actionType_in_Sample = QtWidgets.QAction(MainWindow)
        self.actionType_in_Sample.setObjectName("actionType_in_Sample")
       
        self.PhenoButton.clicked.connect(self.makeGenotype)
        self.copyParent.clicked.connect(self.make_offspring)
        self.copyParent_2.clicked.connect(self.checkShannonProb)
        self.copyParent_5.clicked.connect(self.compare_shan_bin)
        self.copyParent_6.clicked.connect(self.fill_shannon_prob)
        self.Fit.clicked.connect(self.Fit_Shannon)
        self.tableWidget.cellClicked.connect(self.showPic)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Shannon Product Tool v4.101"))
        self.label.setText(_translate("MainWindow", "Genotype         |  Individuals"))
        self.PhenoButton.setText(_translate("MainWindow", "MakeGenotypes"))
        self.label_2.setText(_translate("MainWindow", "Type in Alleles for Parent1"))
        self.label_3.setText(_translate("MainWindow", "Alleles should be sperated with commas and loci with semicolons then press make new generation"))
        self.label_4.setText(_translate("MainWindow", "Genotype        |  Individuals"))
        self.label_6.setText(_translate("MainWindow", "P a r e n t   1"))
        self.label_7.setText(_translate("MainWindow", "P a r e n t   2"))
        self.copyParent.setText(_translate("MainWindow", "compute crossover prob"))
        self.label_9.setText(_translate("MainWindow", "Please Enter Steps:"))
        self.label_12.setText(_translate("MainWindow", "If  Parent2 is from the same population as Parent1, the fields of Parent2 may be left empty."))
        self.label_14.setText(_translate("MainWindow", "number of samples "))
        self.label_13.setText(_translate("MainWindow", "input the joint probability"))
        self.label_15.setText(_translate("MainWindow", "Shannon Calculation"))
        self.label_16.setText(_translate("MainWindow", "sample size"))
        self.label_10.setText(_translate("MainWindow", "R e s u l t s"))
        self.copyParent_2.setText(_translate("MainWindow", "Shannon Product"))
        self.label_5.setText(_translate("MainWindow", "Type in Alleles for Parent2"))
        self.label_18.setText(_translate("MainWindow", "crossover probability"))
        self.label_19.setText(_translate("MainWindow", "input map units"))
        self.Fit.setText(_translate("MainWindow", "Fit Shannon Parameters"))
        self.copyParent_5.setText(_translate("MainWindow", "compare distributions"))
        self.label_8.setText(_translate("MainWindow", "Fill in the number of parents with each genotype, then select the sampe size and the number of samples to compute."))
        self.label_17.setText(_translate("MainWindow", "probabilities must sum up to one"))
        self.label_11.setText(_translate("MainWindow", "click on ofsprng to view the plotted charts from the compared distributions  "))
        self.copyParent_6.setText(_translate("MainWindow", "autofill Shanon prob."))
        self.actionInput_Sample.setText(_translate("MainWindow", "Input Sample from file"))
        self.actionMake_Phenotype.setText(_translate("MainWindow", "Make Phenotype"))
        self.actionType_in_Sample.setText(_translate("MainWindow", "Type in Sample"))
        self.Fit.setEnabled(False)
        self.copyParent_6.setEnabled(False)
        self.copyParent_2.setEnabled(False)
        self.copyParent_5.setEnabled(False)

    def makeGenotype(self):

        x.GenoData = self.GeneInputTextEdit.toPlainText().strip()
        x.GenoData2 = self.GeneInputTextEdit_2.toPlainText().strip()
        Genot=phen_to_gen(x.GenoData).genotypes
        Genot2=phen_to_gen(x.GenoData2).genotypes
        if Genot2==[',']:
            Genot2=[]
        if not any(isinstance(i, list) for i in Genot):
            Genlista=[]
            for i in Genot:
                Genlista.append([i])
            Genot=Genlista[:]
        if not any(isinstance(i, list) for i in Genot2):
            Genlista2=[]
            for i in Genot2:
                Genlista2.append([i])
            Genot2=Genlista2[:]
        self.PhenoButton.setText('O.K.')
        self.make_loci()
        self.make_shannon_prob()
        entries=[[i,'0'] for i in Genot]
        entries2=[[i,'0'] for i in Genot2]        
        self.GenoTable.setRowCount(len(entries))
        self.GenoTable.setColumnCount(2)
        self.GenoTable_2.setRowCount(len(entries2))
        self.GenoTable_2.setColumnCount(2)
        for row, i in enumerate(entries):
            i[0]=[k.replace(',','') for k in i[0]]
            i[0]=','.join([str(k) for k in i[0]])
            # if not i[0] in x.Genotypes.keys():
            #     x.Genotypes[i[0]]=None
            for col,j in enumerate(i):
                item1 = QtWidgets.QTableWidgetItem(j)
                self.GenoTable.setItem(row,col,item1)
        for row, i in enumerate(entries2):
            i[0]=[k.replace(',','') for k in i[0]]
            i[0]=','.join([str(k) for k in i[0]])
            # if not i[0] in x.Genotypes2.keys():
            #     x.Genotypes2[i[0]]=None
            for col,j in enumerate(i):
                item2= QtWidgets.QTableWidgetItem(j)
                self.GenoTable_2.setItem(row,col,item2)
        
        self.GenoTable.resizeColumnsToContents()
        self.GenoTable_2.resizeColumnsToContents() 
        self.copyParent.setEnabled(True)
        pass
  
    def make_cMorgan_offspring(self,splitparent):
        split_parent=splitparent
        #choose1= random.uniform(0, 1)
        #choose2= random.uniform(0, 1)
        choose3= [random.uniform(0, 1) for i in range(len(x.CrossProb))]
        chromapairs =[[2,1,0,3],[3,1,2,0],[0,2,1,3],[0,3,2,1]]
        new_split_parent=split_parent[:]
        for n,val in enumerate(x.CrossProb,start=1):
            cp=random.choice(chromapairs)
            if choose3[n-1] <= val[1]:
                for s,v in enumerate(split_parent[n:]):
                    new_split_parent[n+s]=split_parent[n+s][cp[0]]+split_parent[n+s][cp[1]]+split_parent[n+s][cp[2]]+split_parent[n+s][cp[3]]   
            # else:new_split_parent=split_parent[:]
            # if choose2 > val[1][0] and choose2 <= val[1][1]:
            #     cp=random.choice(chromapairs)
            #     split_parent[n]=split_parent[n][cp[0]]+split_parent[n][cp[1]]+split_parent[n][cp[2]]+split_parent[n][cp[3]]   

        return new_split_parent
    
    def get_cross_prob(self):
        listi=[]
        for row in range(self.GenoTable_6.rowCount()):
            it = self.GenoTable_6.item(row, 1)
            if it and it.text():
                listi.append(it.text())
        sumlist=sum([float(n) for n in listi])
        while not sumlist > 0.0:
            self.label_17.setText("Enter map units Please !!!")
            self.label_17.setStyleSheet("color: rgb(255, 0, 0);")
            break 
        else:
            allrowslist=[]
            for row in range(self.GenoTable_6.rowCount()):
                rowlist=[]
                for col in range(self.GenoTable_6.columnCount()):
                    rowlist.append(self.GenoTable_6.item(row,col).text().strip("' "))
                    rowlist[0]=rowlist[0].replace(' and ','')
                x.mapunits.append([rowlist[0],float(rowlist[1])])
                allrowslist.append([rowlist[0],float(rowlist[1])*2])
                # if not rowlist[0] in allrowsdict.keys():
                #     allrowsdict[rowlist[0]]=rowlist[1]
            # distrvalue=[0]
            # for count, value in enumerate(list(allrowsdict.items())):
            #     distrvalue.append(float(value[1])+distrvalue[count])
            #     allrowsdict[value[0]]=[distrvalue[count],float(value[1])+distrvalue[count]]
            x.CrossProb = allrowslist
        return allrowslist

    def fill_shannon_prob(self):
        self.label_17.setText("probabilities must sum up to one")
        self.label_17.setStyleSheet("color: 98,98,98")
        for row in range(self.GenoTable_5.rowCount()):
            itm = self.GenoTable_5.item(row, 0)
            it=itm.text().replace('event ','').replace('-','')
            item=''
            jointProb=1
            for n,a in enumerate(it):
                if a == '0':
                    jointProb*=1-(x.CrossProb[n][1]/2)
                else: jointProb*=(x.CrossProb[n][1]/2)
            item= QtWidgets.QTableWidgetItem(str(jointProb))
            self.GenoTable_5.setItem(row,1,item)
        self.GenoTable_5.resizeColumnsToContents()
        self.copyParent_6.setEnabled(False)
        self.Fit.setEnabled(True)
        self.copyParent_2.setEnabled(True)

    def make_offspring(self):
        cross_prob=self.get_cross_prob()
        self.MakeParentsLists()
        OfSpsteps = self.SampleSizeSpinBox_2.value()
        McSteps = self.SampleSizeSpinBox.value()
        counter = 0
        while counter < McSteps:
            count=0
            OffspringShannon_dict={}
            while count < OfSpsteps:
                Parents=self.chooseParents()
                offsprShan=[]
                for parent in Parents: 
                    split_parent=parent.split(',')
                    for ind, gene in enumerate(split_parent):
                        split_parent[ind]=gene[0]+gene[0]+gene[1]+gene[1]
                    if cross_prob:
                        #split_parentS=self.make_shannon_offsp(split_parent)
                        split_parentS=self.make_cMorgan_offspring(split_parent)
                    else:split_parentS=False
                    haploindex = random.choice([0,1,2,3])
                    if split_parentS:
                        haploS =''
                        for geneS in split_parentS:
                            haploS += geneS[haploindex]
                        offsprShan.append(haploS)

                if len(offsprShan)>0:
                    offspringS=''
                    for k,y in zip(offsprShan[0], offsprShan[1]):
                        offspringS += k+y+','
                    offspringS=offspringS[:-1]
                    if offspringS in OffspringShannon_dict.keys():
                        OffspringShannon_dict[offspringS]+=1
                    else:
                        OffspringShannon_dict[offspringS]=1
                count += 1
            for itm in OffspringShannon_dict.items():
                OffspringShannon_dict[itm[0]]=float(itm[1])#/OfSpsteps
            
            for key in OffspringShannon_dict.items():
                if key[0] in x.AllTogether_dict.keys():
                        x.AllTogether_dict[key[0]].append(key[1])
                else:
                    x.AllTogether_dict[key[0]]=[key[1]]
            del OffspringShannon_dict
            counter += 1
        MeanList=[]#'entries':['gentype','mean','conf interv', 'shann pr']
        
        for key in x.AllTogether_dict.items():
            
            #std_deviation=(np.std((key[1]),ddof=1))#*OfSpsteps
            #variance=std_deviation**2
            #std_error_mean=std_deviation/np.sqrt(len(key[1]))
            logkey1=[np.log(a) for a in key[1]]
            logmean=np.mean(logkey1)
            logvar=np.var(logkey1)
            MV=np.mean(key[1])/self.SampleSizeSpinBox_2.value()
            ## naive method for confidence intervals 
            lowCI=logmean+logvar/2-2.02*np.sqrt(logvar/self.SampleSizeSpinBox.value())
            highCI=logmean+logvar/2+2.02*np.sqrt(logvar/self.SampleSizeSpinBox.value())
            #lowCI=logmean+logvar/2-1.96*np.sqrt(logvar/self.SampleSizeSpinBox.value()+logvar**2/2*(self.SampleSizeSpinBox.value()-1))
            #highCI=logmean+logvar/2+1.96*np.sqrt(logvar/self.SampleSizeSpinBox.value()+logvar**2/2*(self.SampleSizeSpinBox.value()-1))
            
            CI=[np.exp(lowCI)/self.SampleSizeSpinBox_2.value(),np.exp(highCI)/self.SampleSizeSpinBox_2.value()]
            MeanList.append([key[0],'{0:.4f}'.format(MV),'{0:.4f}  ~  {1:.4f}'.format(CI[0],CI[1])])
            #MeanList.append([key[0],'{0:.4f}'.format(MV),'{0:.2f}'.format(std_deviation),'{0:.2f}'.format(variance),'{0:.2f}'.format(std_error_mean)])
        x.meanList=MeanList
        self.copyParent_2.setEnabled(True)
        self.copyParent_6.setEnabled(True)
        self.Fit.setEnabled(True)
        self.make_stats(MeanList)
        self.copyParent.setEnabled(False)

    def make_stats(self,Mlist):
        meanlist=Mlist
        
        self.tableWidget.setRowCount(len(meanlist))
        self.tableWidget.setColumnCount(4)
        self.tableWidget.setHorizontalHeaderLabels(['ofsprg','Mean','confid. intervals','ShanProd'])
        #self.tableWidget.setHorizontalHeaderLabels(['ofsprg','Mean','stdev','var','SEMean','ShanProd'])
        for row, i in enumerate(meanlist):
            for col,j in enumerate(i):
                item = QtWidgets.QTableWidgetItem(str(j))
                self.tableWidget.setItem(row,col,item)
        self.tableWidget.resizeColumnsToContents()
            
    def get_shannon_prob(self):
        listi=[]
        for row in range(self.GenoTable_5.rowCount()):
            it = self.GenoTable_5.item(row, 1)
            if it and it.text():
                listi.append(it.text())
        sumlisti=sum([float(n) for n in listi])
        while sumlisti < 1:
            self.label_17.setText("Probabilities don't sum up to 1 !!!")
            self.label_17.setStyleSheet("color: rgb(255, 0, 0);")
            break 
        else:
            allrowsdict={}
            for row in range(self.GenoTable_5.rowCount()):
                rowlist=[]
                for col in range(self.GenoTable_5.columnCount()):
                    rowlist.append(self.GenoTable_5.item(row,col).text().strip("' "))
                    rowlist[0]=rowlist[0].replace('event ','').replace('-','')
                if not rowlist[0] in allrowsdict.keys():
                    allrowsdict[rowlist[0]]=rowlist[1]
            distrvalue=[0]
            for count, value in enumerate(list(allrowsdict.items())):
                distrvalue.append(float(value[1])+distrvalue[count])
                allrowsdict[value[0]]=[distrvalue[count],float(value[1])+distrvalue[count]]
            x.ShanProb = allrowsdict
        
            return allrowsdict
           
    def make_loci(self):
        loci=[]
        locipair=[]
        for locus in x.GenoData.split(';'):
            loci.append(locus[0].upper())
        #a=[(loci[i],loci[j]) for i in range(len(loci)) for j in range(i+1, len(loci))]
        for i in range(1,len(loci)):
            locipair.append([loci[i-1]+' and '+loci[i], '0.5'])
        self.GenoTable_6.resizeColumnsToContents()
        self.GenoTable_6.setRowCount(len(locipair))
        self.GenoTable_6.setColumnCount(2)
        for row, i in enumerate(locipair):
            for col,j in enumerate(i):
                item = QtWidgets.QTableWidgetItem(j)
                self.GenoTable_6.setItem(row,col,item)
        self.GenoTable_6.resizeColumnsToContents()

    def checkShannonProb(self):
        
        listi=[]
        for row in range(self.GenoTable_6.rowCount()):
            it = self.GenoTable_5.item(row, 1)
            if it and it.text():
                listi.append(it.text())
        if len(listi) == self.GenoTable_6.rowCount():
            self.make_Shannon_product()
        else: 
            self.label_17.setText("Please fill in Shannon probabilities")
            self.label_17.setStyleSheet("color: red")
   
    def make_shannon_prob(self):
        itemlist=[]
        points=len(x.GenoData.split(';'))-1
        prob_points=list(it.product([0,1],repeat=points))
        for i in prob_points:
            event='-'.join(str(s) for s in i)
            itemlist.append(['event ' + event , ''])
        
        self.GenoTable_5.setRowCount(len(itemlist))
        self.GenoTable_5.setColumnCount(2)
        for row, i in enumerate(itemlist):
            for col,j in enumerate(i):
                item = QtWidgets.QTableWidgetItem(j)
                self.GenoTable_5.setItem(row,col,item)
        self.GenoTable_5.resizeColumnsToContents()

    def MakeParentsLists(self):
        SampleSize=0
        for row in range(self.GenoTable.rowCount()):
            SampleSize+=int(self.GenoTable.item(row,1).text())
        for row in range(self.GenoTable.rowCount()):
            rowlist=[]
            for col in range(self.GenoTable.columnCount()):
                rowlist.append(self.GenoTable.item(row,col).text().strip("' "))
            x.Samplefreq[rowlist[0]]=int(rowlist[1])/float(SampleSize)
            x.ShanParents1[rowlist[0]]=int(rowlist[1]) 
            if self.GenoTable_2.rowCount()==0:
                x.Samplefreq2[rowlist[0]]=int(rowlist[1])/SampleSize
                x.ShanParents2[rowlist[0]]=int(rowlist[1])
        if self.GenoTable_2.rowCount()!=0:
            SampleSize2=0
            for row in range(self.GenoTable_2.rowCount()):
                SampleSize2+=int(self.GenoTable_2.item(row,1).text())
            for row in range(self.GenoTable_2.rowCount()):
                rowlist2=[]
                for col in range(self.GenoTable_2.columnCount()):
                    rowlist2.append(self.GenoTable_2.item(row,col).text().strip("' "))
                x.Samplefreq2[rowlist2[0]]=int(rowlist2[1])/SampleSize2
                x.ShanParents2[rowlist2[0]]=int(rowlist2[1])

    def chooseParents(self):
        freq_val_list=list(x.Samplefreq.values())
        freq_key_list=list(x.Samplefreq.keys())
        listi=[]
        for row in range(self.GenoTable_2.rowCount()):
            it = self.GenoTable_2.item(row, 1)
            if it and it.text():
                listi.append(it.text())
        if len(listi)==0:
            Parent1,Parent2=random.choices(freq_key_list,freq_val_list,k=2)
        else:
            freq_val_list2=list(x.Samplefreq2.values())
            freq_key_list2=list(x.Samplefreq2.keys())
            Parent1=random.choices(freq_key_list,freq_val_list,k=1)
            Parent1=Parent1[0]
            Parent2=random.choices(freq_key_list2,freq_val_list2,k=1)
            Parent2=Parent2[0]
        return [Parent1 , Parent2]
   
    def make_Shannon_product(self):
        offspring=[i[0] for i in  x.meanList]
        Genodictlist=[]

        #Parent1=================================================
        GenoData=x.GenoData
        chromo=GenoData.split(';')
        haplo1=[i[0] for i in chromo]
        haplo2=[i[-1] for i in chromo]
        allhaplodict={','.join(haplo1):'',','.join(haplo2):''}
        Genot=phen_to_gen(GenoData).genotypes
        Genodict1={} 
        for b in Genot:
            if tuple(b) not in Genodict1.keys():
                pops=','.join([i.replace(',','') for i in b])
                Genodict1[tuple(b)]=x.Samplefreq[pops]
        GeList=list(Genodict1.keys())
        AllKeylist=[]
        for n,j in enumerate(GeList[:-2]):
            jnp=np.array([c[0].split(',') for c in [[w] for w in j]])
            SumMirrorObj=Genodict1[j]
            keyList=[j]
            if not j in AllKeylist:  
                for t in GeList[n+1:]:
                    npt=np.array([c[0].split(',') for c in [[w] for w in t]])
                    if (jnp[:,0]==npt[:,-1]).all() and (jnp[:,-1]==npt[:,0]).all():
                        SumMirrorObj += Genodict1[t]
                        keyList.append(t)
                AllKeylist.append(keyList) 
            for k in keyList:
                Genodict1[k]=SumMirrorObj/len(keyList)
        Genodictlist.append(Genodict1)
        Genoarray1=[]
        prob_pointsDict={}
        for row in range(self.GenoTable_5.rowCount()):
            rowlist=[]
            for col in range(self.GenoTable_5.columnCount()):
                rowlist.append(self.GenoTable_5.item(row,col).text().strip("' "))
            rl=rowlist[0].replace('event ','').replace('-',',').split(',')
            rowlist[0]=tuple([int(i) for i in rl])
            if not (rowlist[0]) in  prob_pointsDict.keys():
                prob_pointsDict[(rowlist[0])]=float(rowlist[1])
        for g in Genot:
            gt=[i.split(',') for i in g]
            Genoarray1.append(gt)
        Genoarray1 = np.array(Genoarray1)

        #Parent2=========================================
        if x.GenoData2 == '':
            GenoData2=x.GenoData
        else:
            GenoData2=x.GenoData2
        chromo2=GenoData2.split(';')
        haplo12=[i[0] for i in chromo2]
        haplo22=[i[-1] for i in chromo2]
        allhaplo2=[haplo12,haplo22]
        for key in allhaplo2:
            if ','.join(key) not in allhaplodict.keys():
                allhaplodict[','.join(key)]=''
        Genot2=phen_to_gen(GenoData2).genotypes
        Genodict2={} 
        for b in Genot2:
            if tuple(b) not in Genodict2.keys():
                pops=','.join([i.replace(',','') for i in b])
                Genodict2[tuple(b)]=x.Samplefreq2[pops]
        GeList2=list(Genodict2.keys())
        AllKeylist2=[]
        for n,j in enumerate(GeList2[:-2]):
            jnp=np.array([c[0].split(',') for c in [[w] for w in j]])
            SumMirrorObj2=Genodict2[j]
            keyList2=[j]
            if not j in AllKeylist2:
                for t in GeList2[n+1:]:
                    npt=np.array([c[0].split(',') for c in [[w] for w in t]])
                    if (jnp[:,0]==npt[:,-1]).all() and (jnp[:,-1]==npt[:,0]).all():
                        SumMirrorObj2 += Genodict2[t]
                        keyList2.append(t)
                AllKeylist2.append(keyList2)        
                for k in keyList2:
                    Genodict2[k]=SumMirrorObj2/len(keyList2)
        Genodictlist.append(Genodict2)
        Genoarray2=[]
        points=len(chromo2)-1
        prob_points=list(it.product([0,1],repeat=points))
        for g in Genot2:
            gt=[i.split(',') for i in g]
            Genoarray2.append(gt)
        Genoarray2 = np.array(Genoarray2)
        Genoarray=[Genoarray1,Genoarray2]

        #====Shannon Product=============================
        ShannonProductDict={}
        for of in offspring:
            h=of.split(',')
            haplo1=[i[0] for i in h]
            haplo2=[i[-1] for i in h]
            haplo12=[haplo1,haplo2]
            ShannonProduct=0
            for time in [0,1]:
                PeesList=[]
                AllElem=[]
                for z,hp in zip(Genoarray,haplo12):# for parent_i genotypes and for ofspring haplotype_i
                    Pees={}
                    for p in prob_points:
                        if p not in Pees.keys():
                            Pees[p]=[]
                        for k in z: # for each haplotype in Parent_i
                            compList=[k[0][time]]
                            clst=k[1:]
                            for m,n  in enumerate(p):
                                if n == 1:
                                    clst = [i[::-1] for i in clst]
                                    compList.append(clst[m][time])
                                else:   
                                    compList.append(clst[m][time])
                                if compList==hp:
                                    tupk=[]
                                    for q in k:
                                        tupk.append(','.join(q))
                                    if tuple(tupk) not in Pees[p]:
                                        Pees[p].append(tuple(tupk))
                    PeesList.append(Pees)
                
                if not of in x.InferencePees.keys():
                    x.InferencePees[of]=[PeesList]
                else: 
                    x.InferencePees[of].append(PeesList) 

                for Pee,G in zip(PeesList,Genodictlist):
                    Elem=[]
                    for itm in Pee.items():
                        PeeSum=0
                        for v in itm[1]:
                            PeeSum+=G[v]
                        Elem.append(PeeSum*prob_pointsDict[itm[0]])
                    AllElem.append(sum(Elem))
                ShannonProduct += 0.5*np.product(AllElem)
            if of not in ShannonProductDict.keys():
                ShannonProductDict[of] = ShannonProduct 
            else:
                ShannonProductDict[of] += ShannonProduct
        for i in ShannonProductDict.keys():
            if i in offspring:
                shp=ShannonProductDict[i]
                x.meanList[offspring.index(i)].append('{0:.4f}'.format(shp))
            else:
                x.meanList.append([i,'-','-','-','-','{0:.4f}'.format(shp)])
                pass
        self.make_stats(x.meanList)
        results = os.path.dirname(os.path.abspath(__file__))+'\\Results'
        if not os.path.exists(results):
            os.mkdir(results)
        with open (results+'\\results.csv','w') as f01:
            #f01.write('ofsprg,Mean,stdev,var,SEMean,ShanProd')
            w = csv.writer(f01)
            w.writerow( ['genes','mapunits'] )
            w.writerows( x.mapunits )
            w.writerow( ['joint prob','shanon prob'] )
            w.writerows( prob_pointsDict.items())
            w.writerow( ['ofsprg','Mean','stdev','var','SEMean','ShanProd'] )
            w.writerows( x.meanList )
        with open (results+'\\Parent1.csv','w') as f02:
            #f02.write('Genotype ,  # indivuduals')
            w1 = csv.writer(f02)
            w1.writerow(['Genotype' ,  '# indivuduals'])
            w1.writerows(x.ShanParents1.items())
        with open (results+'\\Parent2.csv','w') as f03:
            #f03.write('Genotype ,  # indivuduals')
            w2 = csv.writer(f03)
            w2.writerow(['Genotype' ,  '# indivuduals'])
            w2.writerows(x.ShanParents2.items())
        self.copyParent_5.setEnabled(True)
        self.copyParent_2.setEnabled(False)
        self.Fit.setEnabled(False)
  
    def Fit_Shannon(self): 
        self.label_17.setText("probabilities must sum up to one")
        self.label_17.setStyleSheet("color: 98,98,98")
        offspring=[i[0] for i in  x.meanList]
        Genodictlist=[]
        #Parent1=================================================
        GenoData=x.GenoData
        chromo=GenoData.split(';')
        haplo1=[i[0] for i in chromo]
        haplo2=[i[-1] for i in chromo]
        allhaplodict={','.join(haplo1):'',','.join(haplo2):''}
        Genot=phen_to_gen(GenoData).genotypes
        Genodict1={} 
        for b in Genot:
            if tuple(b) not in Genodict1.keys():
                pops=','.join([i.replace(',','') for i in b])
                Genodict1[tuple(b)]=x.Samplefreq[pops]
        GeList=list(Genodict1.keys())
        AllKeylist=[]
        for n,j in enumerate(GeList[:-2]):
            jnp=np.array([c[0].split(',') for c in [[w] for w in j]])
            SumMirrorObj=Genodict1[j]
            keyList=[j]
            if not j in AllKeylist:  
                for t in GeList[n+1:]:
                    npt=np.array([c[0].split(',') for c in [[w] for w in t]])
                    if (jnp[:,0]==npt[:,-1]).all() and (jnp[:,-1]==npt[:,0]).all():
                        SumMirrorObj += Genodict1[t]
                        keyList.append(t)
                AllKeylist.append(keyList) 
            for k in keyList:
                Genodict1[k]=SumMirrorObj/len(keyList)
        Genodictlist.append(Genodict1)
        Genoarray1=[]
        for g in Genot:
            gt=[i.split(',') for i in g]
            Genoarray1.append(gt)
        Genoarray1 = np.array(Genoarray1)
        #Parent2=========================================
        if x.GenoData2 == '':
            GenoData2=x.GenoData
        else:
            GenoData2=x.GenoData2
        chromo2=GenoData2.split(';')
        haplo12=[i[0] for i in chromo2]
        haplo22=[i[-1] for i in chromo2]
        allhaplo2=[haplo12,haplo22]
        for key in allhaplo2:
            if ','.join(key) not in allhaplodict.keys():
                allhaplodict[','.join(key)]=''
        Genot2=phen_to_gen(GenoData2).genotypes
        Genodict2={} 
        for b in Genot2:
            if tuple(b) not in Genodict2.keys():
                pops=','.join([i.replace(',','') for i in b])
                Genodict2[tuple(b)]=x.Samplefreq2[pops]
        GeList2=list(Genodict2.keys())
        AllKeylist2=[]
        for n,j in enumerate(GeList2[:-2]):
            jnp=np.array([c[0].split(',') for c in [[w] for w in j]])
            SumMirrorObj2=Genodict2[j]
            keyList2=[j]
            if not j in AllKeylist2:
                for t in GeList2[n+1:]:
                    npt=np.array([c[0].split(',') for c in [[w] for w in t]])
                    if (jnp[:,0]==npt[:,-1]).all() and (jnp[:,-1]==npt[:,0]).all():
                        SumMirrorObj2 += Genodict2[t]
                        keyList2.append(t)
                AllKeylist2.append(keyList2)        
                for k in keyList2:
                    Genodict2[k]=SumMirrorObj2/len(keyList2)
        Genodictlist.append(Genodict2)
        Genoarray2=[]
        points=len(chromo2)-1
        prob_points=list(it.product([0,1],repeat=points))
        for g in Genot2:
            gt=[i.split(',') for i in g]
            Genoarray2.append(gt)
        Genoarray2 = np.array(Genoarray2)
        Genoarray=[Genoarray1,Genoarray2]
        #====Shannon Events=============================
        for of in offspring:
            h=of.split(',')
            haplo1=[i[0] for i in h]
            haplo2=[i[-1] for i in h]
            haplo12=[haplo1,haplo2]
            for time in [0,1]:
                PeesList=[]
                for z,hp in zip(Genoarray,haplo12):# for parent_i genotypes and for ofspring haplotype_i
                    Pees={}
                    for p in prob_points:
                        if p not in Pees.keys():
                            Pees[p]=[]
                        for k in z: # for each haplotype in Parent_i
                            compList=[k[0][time]]
                            clst=k[1:]
                            for m,n  in enumerate(p):
                                if n == 1:
                                    clst = [i[::-1] for i in clst]
                                    compList.append(clst[m][time])
                                else:   
                                    compList.append(clst[m][time])
                                if compList==hp:
                                    tupk=[]
                                    for q in k:
                                        tupk.append(','.join(q))
                                    if tuple(tupk) not in Pees[p]:
                                        Pees[p].append(tuple(tupk))
                    PeesList.append(Pees)
                
                if not of in x.InferencePees.keys():
                    x.InferencePees[of]=[PeesList]
                else: 
                    x.InferencePees[of].append(PeesList) 

        for ofspng in x.InferencePees.items():
            for P in ofspng[1]:
                for n,par in enumerate(P):
                    for itm in par.items():
                        parSum=0
                        for elem in itm[1]:
                            parSum+=Genodictlist[n][elem]

                        itm[1].append(parSum)
        noevents=[]
        pureEvents=[]
        bprod=[]
        points=len(x.GenoData.split(';'))-1
        prob_points=list(it.product([0,1],repeat=points))
        for i in prob_points:
            if i.count(1)==1: 
                pureEvents.append(i)
        b=[0,1]
        blist=[b for i in range(len(pureEvents))]
        x0=[0.1]*len(pureEvents)                
        results=[]
        res = minimize(self.Allfunct,x0=x0,method='SLSQP',bounds=blist,tol=1e-6)
        #res = minimize(self.Allfunct,x0=x0,method='Powell',bounds=blist)
        results=res.x
        ################################
        p={}
        count=0
        for i in prob_points:
            event='-'.join(str(s) for s in i)
            if i.count(1)==1:
                p['event {}'.format(event)]=res.x[count] 
                count += 1
        for q in prob_points:
            if q.count(1)>1:
                event='-'.join(str(s) for s in q)
                bprod=[]
                for k,b in enumerate(q):
                    if b==1:
                        for tpl in pureEvents:
                            if tpl[k]==1:
                                bprod.append(tpl)
                bprd=1                
                for e in bprod:
                    evnt='-'.join(str(s) for s in e)
                    bprd *= p['event {}'.format(evnt)]
                p['event {}'.format(event)]=bprd
            
        for f in prob_points:
            event='-'.join(str(s) for s in f)
            if f.count(1)==0: 
                if 1-sum(p.values())<0:
                    return 1
                p['event {}'.format(event)]=1-sum(p.values())
                noevents.append(f)
        ##############################################################      
        itemlist=[]
        for i in sorted(list(p.items())):
            #event='-'.join(str(s) for s in i)
            itemlist.append([i[0] , '{0:.4f}'.format(i[1])])
        self.GenoTable_5.setRowCount(len(itemlist))
        self.GenoTable_5.setColumnCount(2)
        for row, i in enumerate(itemlist):
            for col,j in enumerate(i):
                item = QtWidgets.QTableWidgetItem(j)
                self.GenoTable_5.setItem(row,col,item)
        self.GenoTable_5.resizeColumnsToContents()
        
        self.copyParent_6.setEnabled(True) 
        self.Fit.setEnabled(False)
              
    def Allfunct(self,x0):
        x0=x0
        result=0
        for ofspng in x.InferencePees:
            res=self.funct(x0,ofspng)
            result+=res
        return result

    def funct(self,x0,ofspng):
        ofspg=ofspng
        points=len(x.GenoData.split(';'))-1
        prob_points=list(it.product([0,1],repeat=points))
        noevents=[]
        pureEvents=[]
        bprod=[]
        p={}
        count=0
        for i in prob_points:
            event='-'.join(str(s) for s in i)
            if i.count(1)==1:
                p['p{}'.format(event)]=x0[count] 
                pureEvents.append(i)
                count += 1
        for q in prob_points:
            if q.count(1)>1:
                event='-'.join(str(s) for s in q)
                bprod=[]
                for k,b in enumerate(q):
                    if b==1:
                        for tpl in pureEvents:
                            if tpl[k]==1:
                                bprod.append(tpl)
                bprd=1                
                for e in bprod:
                    evnt='-'.join(str(s) for s in e)
                    bprd *= p['p{}'.format(evnt)]
                p['p{}'.format(event)]=bprd
            
        for f in prob_points:
            event='-'.join(str(s) for s in f)
            if f.count(1)==0: 
                if 1-sum(p.values())<0:
                    return 1
                p['p{}'.format(event)]=1-sum(p.values())
                noevents.append(f)
        Allhalf=[]
        for half in x.InferencePees[ofspg]:
            parAllVal=[]
            for parent in half:
                parVal=[]
                for pee in parent:
                    xpa='-'.join(str(s) for s in pee)
                    parVal.append(('p{}'.format(xpa),parent[pee][-1])) 
                parAllVal.append(parVal)
            Allhalf.append(parAllVal)
        halfinf=[]
        for infP in Allhalf:
            hsumip = 1
            for m in infP:
                sumip=0
                for ip in m:
                    sumip += p[ip[0]]*ip[1]
                hsumip *= sumip
            halfinf.append(hsumip*0.5) 
        result=0 
        for lst in x.meanList:
            if ofspg in lst:
                result= float(lst[1])
        #print(p.items())
        return (sum(halfinf)-result)**2
        
    def compare_shan_bin(self):

        if not os.path.exists(os.path.dirname(__file__)+'\\BarCharts'):
            os.mkdir(os.path.dirname(__file__)+'\\BarCharts')
        filepath2 = os.path.dirname(__file__)+'\\BarCharts'
        shlist=[i[0] for i in x.meanList]
        AllBarData={}
        count=0
        # Filenames={}
        for key in x.AllTogether_dict.items():
            shcompDict={}
            tally=[[t,key[1].count(t)] for t in set(key[1])]
            tallyX=[int(i[0]) for i in tally]
            tallyY=[i[1] for i in tally]
            sp=float(x.meanList[shlist.index(key[0])][1])
            #binomialbin=[]
            for k in range(0,max(tallyX)+10):
                if k not in shcompDict.keys():
                    if k not in tallyX:
                        shcompDict[k]=[0]
                    else: shcompDict[k]= [tallyY[tallyX.index(k)]]

                n=self.SampleSizeSpinBox_2.value() ##sample size
                bi=scipy.stats.binom(n,sp)
                ybi=bi.pmf(k)
                #a=np.math.factorial(n) // np.math.factorial(k) // np.math.factorial(n - k)
                shcompDict[k].append(ybi)
                #shcompDict[k].append(a*sp**k*(1-sp)**(n-k))
                #binomialbin.append(a*sp**k*(1-sp)**(n-k))
                #binomialbin.append(ybi)
            indxlst=[]
            Shannon = []
            Binomial =[]
            for i in shcompDict:
                if shcompDict[i][0] != 0 or shcompDict[i][1] > 1e-7:
                    Shannon.append(shcompDict[i][0])
                    Binomial.append(shcompDict[i][1]*self.SampleSizeSpinBox.value())
                    indxlst.append(list(shcompDict.keys())[i])
            filename2=filepath2 +'\\Bar_{}'.format(key[0])
            # Filenames[key[0]]=filename2

            # save to csv file##########################
            indexlist=[]
            for i in range(0,len(indxlst),5):
                indexlist.extend([indxlst[i],'','','',''])

            header=['Shannon','Binomial']
            shanFname=filepath2+'\\compare_{}.{}.csv'.format(key[0],count)
            newlst=[[i,z] for i,z in zip(Shannon,Binomial)]
            with open(shanFname, 'w', newline='') as f:
                writer = csv.writer(f,delimiter='|')
                writer.writerow(header)
                writer.writerows(newlst)
            f.close()
            # create bar plot ##########################
            
            index = np.arange(len(Shannon))
            bar_width = 0.35
            opacity = 1

            plt.bar(index, Shannon, bar_width,
            alpha=opacity,
            color='b',
            label='Observations')

            plt.bar(index + bar_width, Binomial, bar_width,
            alpha=opacity,
            color='g',
            label='Binomial')
            plt.xlabel('# of observations')
            plt.ylabel('Scores')
            plt.title('{}'.format(key[0]))
            plt.xticks(index + bar_width, list(indexlist))
            plt.xticks(fontsize=8,rotation=30)
            plt.legend()

            plt.tight_layout()
            plt.savefig(filename2+'.{}.jpg'.format(count),dpi=300)
            #plt.show()
            plt.close('all')
            count+=1
            AllBarData[key[0]] = shcompDict
            
    def showPic(self):
        col=self.tableWidget.currentItem().column()
        if col == 0: 
            row = self.tableWidget.currentItem()
            rownum = self.tableWidget.currentItem().row()
            x.fpath = os.path.dirname(__file__)+'\\BarCharts'+'\\'+ 'Bar_'+ row.text()+'.{}'.format(rownum) + '.jpg'
            
            self.Datainput = QtWidgets.QWidget()
            self.ui = Popup()
            self.ui.setupUi(self.Datainput)
            self.Datainput.show()

        

class Popup(object):
    def setupUi(self, Popup):
        Popup.setObjectName("Popup")
        Popup.resize(655, 488)
        self.label = QtWidgets.QLabel(Popup)
        self.label.setGeometry(QtCore.QRect(10, 10, 655, 488))
        self.label.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.label.setObjectName("label")
        self.label.setScaledContents(True)
        self.pixmap = QPixmap(x.fpath)
        self.label.setPixmap(self.pixmap)
        self.closeButton = QtWidgets.QPushButton(Popup)
        self.closeButton.setGeometry(QtCore.QRect(10, 490, 640, 21))
        self.closeButton.setObjectName("pushButton")
        self.closeButton.clicked.connect(Popup.close)
        self.retranslateUi(Popup)
        QtCore.QMetaObject.connectSlotsByName(Popup)

    def retranslateUi(self, Popup):
        _translate = QtCore.QCoreApplication.translate
        Popup.setWindowTitle(_translate("Popup", "compare distributions"))
        self.closeButton.setText(_translate("Popup", "close image"))

class PhenoG():

    def phenomenon(self):
        
        GenTypeAB=phen_to_gen(x.GenoData).genotypes
        if not any(isinstance(i, list) for i in GenTypeAB):
            Genlista=[]
            for i in GenTypeAB:
                Genlista.append([i])
            GenTypeAB=Genlista 
        if x.checkbox:
            start=x.SampleSize-int(x.SampleSize*15/100)
            finish= x.SampleSize+int(x.SampleSize*15/100)
            elem=range(start,finish)
            Sample_Size=random.choice(elem)
        else:  
            Sample_Size = x.SampleSize
            
        if len(x.Weights.values())>0:
            Weights=x.Weights
        else: Weights=None
        a=Phenomenon(GenTypeAB,SampleSize=Sample_Size,Weights=Weights)
        b=a.PhenProb()
        
        while not all(i>0 for i in b[0].values()):
            
            b=a.PhenProb()
        c=a.PopSampleLista
        d=a.PhenTypelist
        x.FinalData=[b,c,d]

class Controller():
    def __init__(self):
        self.SampleSize=0
        self.GenoData=''
        self.GenoData2=''
        self.FinalData=[]
        self.Weights={}
        self.Genotypes={}
        self.Genotypes2={}
        self.checkbox=''
        self.ShanProb={}
        self.mapunits=[]
        self.CrossProb={}
        self.ShanParents1={}
        self.ShanParents2={}
        self.meanList=None
        self.AllTogether_dict={}
        self.InferencePees={}
        self.fpath=''
        self.Samplefreq={}
        self.Samplefreq2={}

    def main(self):
        
        app = QtWidgets.QApplication(sys.argv)
        MainWindow = QtWidgets.QMainWindow()
        ui = Ui_MainWindow()
        ui.setupUi(MainWindow)
        MainWindow.show()
        sys.exit(app.exec_())   

if __name__ == "__main__":
    x=Controller()
    x.main()
    
    
    # app = QtWidgets.QApplication(sys.argv)
    # MainWindow = QtWidgets.QMainWindow()
    # ui = Ui_MainWindow()
    # ui.setupUi(MainWindow)
    # MainWindow.show()
    # sys.exit(app.exec_())

