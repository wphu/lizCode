import sys
import random
import matplotlib
import h5py as h5
import numpy as np
from math import isinf, isnan
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.animation as animation

from matplotlib.ticker import MaxNLocator
from matplotlib import cm

from PyQt5 import QtCore
from PyQt5.QtWidgets import QSizePolicy, QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel, QTextEdit, QPushButton, QTabWidget
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem

from collect import collect

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, data4d=None, title = None):
        self.fig = Figure()
        self.fig.subplots_adjust(top=0.85, bottom=0.2, left=0.2)
        self.axes1 = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes1.hold(False)

        self.data4d = data4d

        self.compute_initial_figure(self.data4d, title)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self, data = None, title = None):
        pass

class MyStaticMplCanvas1D(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self, data4d = None, title = None):
        # data2d[itime, ix]
        self.data2d = data4d[:,0,0,:]
        self.title = title
        dx = 1.0
        nx = self.data2d.shape[1]
        self.x = np.linspace(0,100.0,nx)

        self.line1, = self.axes1.plot(self.x, self.data2d[0])
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel(self.title)
        self.axes1.set_title(self.title)
        self.setLimits()


    def update_figure(self, itime):
        line1, = self.axes1.plot(self.x, self.data2d[itime,:])
        self.axes1.set_title(self.title)
        self.setLimits()
        return [line1]

    def setLimits(self):
        if isinf(self.data2d.min()) or isnan(self.data2d.min()) or isinf(self.data2d.max()) or isnan(self.data2d.max()):
            pass
        else:
            self.axes1.axis([self.x.min(), self.x.max(), self.data2d.min(), self.data2d.max()])

    def save_animation(self):
        self.ani = animation.FuncAnimation(self.fig, self.update_figure, self.data2d.shape[0], blit=False, interval=500, repeat=False)
        self.ani.save('movie.gif', writer='imagemagick')

class MyStaticMplCanvas2D(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self, data4d = None, title = None):
        self.fig.clear()
        self.fig.subplots_adjust(top=0.85, bottom=0.2, left=0.2)
        self.axes1 = self.fig.add_subplot(111)

        self.data3d = data4d[0,:,:,:]
        self.title = title
        nx = self.data3d.shape[1]
        ny = self.data3d.shape[2]
        dx=1.0
        dy=1.0
        self.x,self.y=np.mgrid[slice(dx,dx*(nx+1),dx), slice(dy,dy*(ny+0.5),dy)]
        levels=MaxNLocator(nbins=100).tick_values(self.data3d.min(),self.data3d.max())

        if(self.data3d.min() == self.data3d.max()):
        	ticks_val=np.linspace(self.data3d.min(),self.data3d.max()+1.0,5)
        else:
        	ticks_val=np.linspace(self.data3d.min(),self.data3d.max(),5)
        self.cf = self.axes1.contourf(self.x,self.y,self.data3d[10],cmap=cm.get_cmap('jet'),levels=levels)
        self.fig.colorbar(self.cf,ticks=ticks_val)



        self.axes1.axis([self.x.min(),self.x.max(),self.y.min(),self.y.max()])
        self.axes1.set_yticks(np.arange(0,self.y.max(),100))
        self.axes1.set_title(self.title)
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('y(mm)')
        self.axes1.set_aspect('equal', adjustable='box')

    def update_figure(self, itime):
        self.fig.clear()
        self.fig.subplots_adjust(top=0.85, bottom=0.2, left=0.2)
        self.axes1 = self.fig.add_subplot(111)
        levels=MaxNLocator(nbins=100).tick_values(self.data3d.min(),self.data3d.max())

        if(self.data3d.min() == self.data3d.max()):
        	ticks_val=np.linspace(self.data3d.min(),self.data3d.max()+1.0,5)
        else:
        	ticks_val=np.linspace(self.data3d.min(),self.data3d.max(),5)

        self.cf = self.axes1.contourf(self.x,self.y,self.data3d[itime],cmap=cm.get_cmap('jet'),levels=levels)
        self.fig.colorbar(self.cf,ticks=ticks_val)



        self.axes1.axis([self.x.min(),self.x.max(),self.y.min(),self.y.max()])
        self.axes1.set_yticks(np.arange(0,self.y.max(),100))
        self.axes1.set_title(self.title)
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('y(mm)')
        self.axes1.set_aspect('equal', adjustable='box')

    def save_animation(self):
        self.ani = animation.FuncAnimation(self.fig, self.update_figure, self.data3d.shape[0], blit=False, interval=500, repeat=False)
        self.ani.save('movie.gif', writer='imagemagick')

     
class PlotWidget(QWidget):
    """docstring for PlotWidget"""
    def __init__(self, parent = None):
        super(PlotWidget, self).__init__()
        self.parent = parent

        # data4d has four dimensions, the first is time
        try:
            self.prefix = 'grid'
            self.data4d = collect("is_wall", prefix = self.prefix)
            self.dataName = "geometry"
        except:
            self.prefix = 'data'
            self.data4d = collect("/Fields/", "Phi_global_avg", prefix = self.prefix)
            self.dataName = "electric potential"          

        sizePolicy = QSizePolicy();
        sizePolicy.setHorizontalPolicy(QSizePolicy.Expanding);
        sizePolicy.setVerticalPolicy(QSizePolicy.Expanding);
        self.setSizePolicy(sizePolicy);

        self.plotVboxlayout = QVBoxLayout(self)

        self.tab_widget = QTabWidget(self)
        self.plotVboxlayout.addWidget(self.tab_widget)
        
        # matplotlib figure tab
        #self.matplotlib_widget = MatplotlibWidget(self.tab_widget,self.data4d[0,0,:,:])
        if self.data4d.shape[2] == 1:
            self.sc = MyStaticMplCanvas1D(self, self.data4d, title = self.dataName)
        else:
            self.sc = MyStaticMplCanvas2D(self, self.data4d, title = self.dataName)
        self.sc.draw()
        self.tab_widget.addTab(self.sc, "Figures")
        
        # data tab
        self.data_widget = QTableWidget(self)
        self.set_data_widget(self.data4d[0,0,:,:])
        self.tab_widget.addTab(self.data_widget, "data")


        #> The slider
        self.sp_widget = QWidget(self)
        self.sp_layout = QHBoxLayout(self.sp_widget)

        self.label = QLabel("time")
        self.plainTextEdit = QTextEdit("0")
        self.plainTextEdit.setMaximumHeight(20)
        self.plainTextEdit.setMaximumWidth(100)
        self.sp = QSlider(QtCore.Qt.Horizontal)
        self.sp.setMinimum(0)
        self.sp.setMaximum(self.data4d.shape[0]-1)
        self.sp.setTickPosition(QSlider.TicksBelow)
        self.sp.setTickInterval(1)
        self.sp.valueChanged.connect(self.timeChange)
        self.sp_layout.addWidget(self.label)
        self.sp_layout.addWidget(self.plainTextEdit)
        self.sp_layout.addWidget(self.sp)

        self.save_widget = QWidget(self)
        self.save_layout = QHBoxLayout(self.save_widget)
        self.button_saveFig = QPushButton("Save Figure")
        self.button_saveAnimation = QPushButton("Save Animation")
        self.save_layout.addWidget(self.button_saveFig)
        self.save_layout.addWidget(self.button_saveAnimation)
        self.button_saveAnimation.clicked.connect(self.save_animation)
        self.button_saveFig.clicked.connect(self.save_fig)


        #self.plotVboxlayout.addWidget(self.sc)
        self.plotVboxlayout.addWidget(self.sp_widget)
        self.plotVboxlayout.addWidget(self.save_widget)

    def reloadData(self, prefix, dataSetFullName):
        self.prefix = prefix
        self.dataName = dataSetFullName.rsplit('/')[-1]
        self.data4d = collect(dataSetFullName, prefix=self.prefix)
        self.sc.compute_initial_figure(self.data4d, self.dataName)
        self.sc.draw()
        self.set_data_widget(self.data4d[0,0,:,:])
        self.sp.setMaximum(self.data4d.shape[0]-1)
        self.sp.setValue(0)



    def timeChange(self):
        t = self.sp.value()
        self.sc.update_figure(t)
        self.sc.draw()
        self.set_data_widget(self.data4d[t,0,:,:])
        self.plainTextEdit.setText(str(self.sp.value()))

    def set_data_widget(self, data2d):
        self.data_widget.setRowCount(data2d.shape[0])
        self.data_widget.setColumnCount(data2d.shape[1])
        for i in range(0, data2d.shape[0]):
            for j in range(0, data2d.shape[1]):
                newItem = QTableWidgetItem( str( round(data2d[i,j]) ) )
                self.data_widget.setItem(i, j, newItem)

    def save_fig(self):
        filename = "fig" + str(self.sp.value()) + ".pdf"
        self.sc.fig.savefig(filename)

    def save_animation(self):
        self.sc.save_animation()  
