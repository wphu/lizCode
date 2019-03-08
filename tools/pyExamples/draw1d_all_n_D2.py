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


n0 = 1.0e19

x_step = 40

##read data from file
path = "Re0.85/data"
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
fig=plt.figure(figsize=(10,6))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)

# ===================================== Ne ================================
ax0=fig.add_subplot(2,1,1)
##============ case 0 =============
path = "Re0.85/data"
val = collect("/Fields/", "Rho_global_e_avg", path = path)
val = val / n0

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$n_u \mathrm{=1.0\times 10^{19}m^{-3}}$', linestyle = linestyles[0])


##============ case 1 ==========
path = "Re0.85-n1.5/data"
val = collect("/Fields/", "Rho_global_e_avg", path = path)
val = val / n0

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$n_u \mathrm{=1.5\times 10^{19}m^{-3}}$', linestyle = linestyles[1])


##============ case 2 ============
path = "Re0.85-n2.0/data"
val = collect("/Fields/", "Rho_global_e_avg", path = path)
val = val / n0

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$n_u \mathrm{=2.0\times 10^{19}m^{-3}}$', linestyle = linestyles[2])


ax0.grid(True)
ax0.legend(loc = 1, framealpha=1, fontsize = legend_fontsize)
ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((0.0, 30.0))

major_ticks = np.arange(0, 10.01, 2.0)
ax0.set_yticks(major_ticks)

ax0.set_ylabel(r"$n_\mathrm{e} \ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)



# =============================== nD ================================
ax0=fig.add_subplot(2,1,2)
##============ case 0 ===============
path = "Re0.85/data"
val = collect("/Fields/", "Rho_global_D_avg", path = path)
val = val / n0

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$n_u \mathrm{=1.0\times 10^{19}m^{-3}}$', linestyle = linestyles[0])

##============ case 1 ============
path = "Re0.85-n1.5/data"
val = collect("/Fields/", "Rho_global_D_avg", path = path)
val = val / n0

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$n_u \mathrm{=1.5\times 10^{19}m^{-3}}$', linestyle = linestyles[1])

##============ case 2 ============
path = "Re0.85-n2.0/data"
val = collect("/Fields/", "Rho_global_D_avg", path = path)
val = val / n0

val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$n_u \mathrm{=2.0\times 10^{19}m^{-3}}$', linestyle = linestyles[2])


ax0.grid(True)
#ax0.legend(loc = 1)
ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((0.0, 90.0))

major_ticks = np.arange(0, 8.1, 2.0)
ax0.set_yticks(major_ticks)


ax0.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$n_\mathrm{D} \ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)


fig.savefig("all_n_D.png", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
