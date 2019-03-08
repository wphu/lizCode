from template import *
from collect import collect

import sys

itime = 0
if len(sys.argv) == 2:
	t = int(sys.argv[1])
elif len(sys.argv) > 2:
	print("error: should have one argument, as time step")
else:
	t = itime


legend_fontsize = 14

x_step = 40

##read data from file
path = "Re0.8/data"
val = collect("/Fields/", "Phi_global_avg", path = path)
nx = val.shape[3]

dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)

amplification_factor = 100.0
x = x * amplification_factor
x_less = x[::x_step]

xmin = x.min()
xmax = x.max()


##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

# ===================================== Te ================================
ax0=fig.add_subplot(3,1,1)
##============ case 0 =============
path = "Re0.8/data"
val = collect("/Fields/", "T_global_e_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.80', linestyle = linestyles[0])


##============ case 1 ==========
path = "Re0.85/data"
val = collect("/Fields/", "T_global_e_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.85', linestyle = linestyles[1])


##============ case 2 ============
path = "Re0.9/data"
val = collect("/Fields/", "T_global_e_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.90', linestyle = linestyles[2])


ax0.grid(True)
ax0.legend(loc = 1, framealpha=1, fontsize = legend_fontsize)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 61.0))

major_ticks = np.arange(0, 61, 20)
ax0.set_yticks(major_ticks)

ax0.set_ylabel(r"$T_\mathrm{e} \ \mathrm{(eV)}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)


# ===================================== Ti ================================
ax0=fig.add_subplot(3,1,2)
ax0.yaxis.set_major_formatter(yformatter)
##============ case 0 =============
path = "Re0.8/data"
val = collect("/Fields/", "T_global_D1_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.80$', linestyle = linestyles[0])

##============ case 1 ==========
path = "Re0.85/data"
val = collect("/Fields/", "T_global_D1_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.85$', linestyle = linestyles[1])


##============ case 2 ============
path = "Re0.9/data"
val = collect("/Fields/", "T_global_D1_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.9$', linestyle = linestyles[2])



ax0.grid(True)
#ax0.legend(loc = 1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 230.1))
ax0.set_ylabel(r"$T_{\mathrm{D^+}} \ \mathrm{(eV)}$", fontsize = label_fontsize)

major_ticks = np.arange(0, 230.1, 50)
ax0.set_yticks(major_ticks)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)




# =============================== Potential ================================
ax0=fig.add_subplot(3,1,3)
##============ case 0 ===============
path = "Re0.8/data"
val = collect("/Fields/", "Phi_global_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.80$', linestyle = linestyles[0])


##============ Case1 ============
path = "Re0.85/data"
val = collect("/Fields/", "Phi_global_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.85$', linestyle = linestyles[1])


##============ Case2 ============
path = "Re0.9/data"
val = collect("/Fields/", "Phi_global_avg", path = path)

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.90$', linestyle = linestyles[2])


ax0.grid(True)
#ax0.legend(loc = 1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 80.1))

major_ticks = np.arange(0, 80.1, 20)
#minor_ticks = np.arange(0, 31, 5)
ax0.set_yticks(major_ticks)
#ax0.set_yticks(minor_ticks, minor=True)


#legend1=ax0.legend(loc=(.6,.76),fontsize=16)
#ax0.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
#ax0.set_yticks(np.arange(0,y.max(),100))
ax0.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$\phi\ \mathrm{(V)}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)


fig.savefig("all_TP_D.png", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
