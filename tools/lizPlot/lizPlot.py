#!/usr/bin/python

import sys
import os
import random
import matplotlib
import h5py as h5
import numpy as np
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QAction, QFileDialog, QTabWidget, QVBoxLayout, QHBoxLayout, QSizePolicy, QMessageBox, QWidget, QSplitter


from HdfTreeModel import *
from PlotWidget import *
from HdfTreeWidget import *

class ApplicationWindow(QMainWindow):

    sigOpen = QtCore.pyqtSignal(list, str)

    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.createActions()
        self.createMenus()


        #> the main layout: left part (for the hdf5 tree view) and right part (for the figures and data)
        self.main_widget = QWidget(self)

        #> if not putting splitter into a layout, the widgets in splitter do not fill the main windows
        #> (may exceed the app windows, so that the figures are partially shown ).
        layout = QVBoxLayout(self.main_widget)
        hSplitter = QSplitter(self.main_widget)
        layout.addWidget(hSplitter)

        #number of space dimension
        self.n_dim = 1

        #> the left part: the hdf5 tree view
        self.filename_data = "data/data0.h5"
        self.filename_grid = "data/grid.h5"

        self.tree_widget = HDFTreeWidget(self.main_widget)
        self.tree_model = HDFTreeModel([])
        if os.path.exists(self.filename_data):
            self.tree_model.openFile(self.filename_data, 'r+')
            h5_file = h5.File(self.filename_data)
            self.n_dim = h5_file.attrs['n_dim']
            print("n_dim: ", self.n_dim)
        if os.path.exists(self.filename_grid):
            self.tree_model.openFile(self.filename_grid, 'r+')
            h5_file = h5.File(self.filename_grid)
            self.n_dim = h5_file.attrs['n_dim']
            print("n_dim: ", self.n_dim)
        self.tree_widget.setModel(self.tree_model)
        hSplitter.addWidget(self.tree_widget)

        self.tree_widget.doubleClicked.connect(self.reload)
        self.sigOpen.connect(self.tree_widget.openFiles)
        
        #> plotWidget
        if self.n_dim == 1:
            self.plot_widget = PlotWidget1D(self.main_widget)
        elif self.n_dim == 2:
            self.plot_widget = PlotWidget2D(self.main_widget)
        elif self.n_dim == 3:
            self.plot_widget = PlotWidget3D(self.main_widget)

        hSplitter.addWidget(self.plot_widget)

        hSplitter.setStretchFactor(0, 2)
        hSplitter.setStretchFactor(1, 3)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.statusBar().showMessage("All hail matplotlib!", 2000)

    def reload(self, index):
    	item = self.tree_model.getItem(index)
    	if (item is not None) and item.isDataset():
            self.plot_widget.reloadData(item.prefix, item.name)

    def createActions(self):
        self.openFileReadOnlyAction = QAction('Open file(s) readonly', self,
                                   # shortcut=QtGui.QKeySequence.Open,
                                   statusTip='Open an HDF5 file for reading',
                                      triggered=self.openFilesReadOnly)
        self.quitAction = QAction('&Quit', self,
                                  statusTip='Quit dataviz',
                                  triggered=self.fileQuit)
        self.aboutAction = QAction('&About', self,
                                  statusTip='Information about the code',
                                  triggered=self.about)


    def createMenus(self):
        self.menuBar().setVisible(True)

        self.fileMenu = self.menuBar().addMenu('&File')
        self.fileMenu.addAction(self.openFileReadOnlyAction)
        self.fileMenu.addAction(self.quitAction)
        self.menuBar().addSeparator()

        self.aboutMenu = self.menuBar().addMenu('&About')
        self.aboutMenu.addAction(self.aboutAction)



    def openFilesReadOnly(self, filePaths=None):
        if filePaths is None or filePaths is False:
            self.fileDialog = QFileDialog(None, 'Open file(s) read-only', './',
                                                 'HDF5 file (*.h5 *.hdf);;All files (*)')
            self.fileDialog.show()
            self.fileDialog.filesSelected.connect(self.openFilesReadOnly)
            return
        # filePaths = QtGui.QFileDialog.getOpenFileNames(self,
        #                                          'Open file(s)', self.lastDir,
        #                                          'HDF5 file (*.h5 *.hdf);;All files (*)')
        filePaths = [str(path) for path in filePaths]   # python2/qt4 compatibility
        if len(filePaths) == 0:
            return
        self.lastDir = QtCore.QFileInfo(filePaths[-1]).dir().absolutePath()
        # TODO handle recent files
        self.sigOpen.emit(filePaths, 'r')


    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QMessageBox.about(self, "About",
"""The code purpose is to draw figures and show data for dataset of hdf5 files
The python code mainly uses the libraries: PyQt5, matplotlib, h5py

Email: wanpengh@gmail.com
"""
)

if __name__ == '__main__':
    app = QApplication(sys.argv)

    aw = ApplicationWindow()
    aw.setWindowTitle("lizPlot")
    aw.show()
    aw.resize(1000, 600)
    #sys.exit(qApp.exec_())
    app.exec_()
