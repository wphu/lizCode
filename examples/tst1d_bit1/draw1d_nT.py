##>>>The code is used to read data from hdf5 file
##>>>and plot on the screen and output figure file using matplotlib-python

import Tkinter as tk
from Tkinter import *

import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import numpy as np
from numpy import arange, sin, pi

import ConfigParser

import h5py as h5
import numpy as np
import matplotlib.pyplot as plt


font={	'family' : 'sans-serif',
	'weight' : 'bold',
	'size' : 8,
	}

mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams['font.weight'] = 'bold'

mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['lines.linewidth'] = 2.0







##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 10


##read data from file
f=h5.File("data_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[3]


dx = 2.0e-2  # unit (mm)
x = np.linspace(0, nx * dx, nx)

xmin = x.min()
xmax = x.max()



##============rho======================================================
sp_temp1=fig.add_subplot(2,1,1)

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "D+1")

val = f["/Fields/Rho_global_D_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "D")



val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C")


val = f["/Fields/Rho_global_C1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C+1")

val = f["/Fields/Rho_global_C2_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C+2")

val = f["/Fields/Rho_global_C3_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C+3")


sp_temp1.legend()
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Number density')



##============ Temperature ======================================================
sp_temp1=fig.add_subplot(2,1,2)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "D+1")

val = f["/Fields/T_global_D_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "D")


val = f["/Fields/T_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C")


val = f["/Fields/T_global_C1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C+1")


val = f["/Fields/T_global_C2_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C+2")


val = f["/Fields/T_global_C3_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "C+3")

sp_temp1.legend()
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Temperature (eV)')



fig.savefig("Profiles.png", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])