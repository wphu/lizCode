from template import *
from collect import collect

t = 0
dt = 0.5e-12 * 500000 / 1.0e-6   # unit is us
dt = 1

N0 = 1.0e5
pFlux0 = 1.0e23
hFlux0 = 1.0e7


val1 = collect("/Diagnostic/", "particle_number")

#the first value is initial value, not calculated
nx = (val1.shape)[0] - 1
x = np.linspace(0,nx*dt,nx)

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.6,hspace=0.55)


#============Total particle number=======================================
val = collect("/Diagnostic/", "particle_number")
val = val / N0
val1_1d = val[1:, 0, 0, 0]
val2_1d = val[1:, 0, 0, 1]
print(val.shape)

ax0=fig.add_subplot(3,1,1)


ax0.yaxis.set_major_formatter(yformatter)

line0=ax0.plot(x, val1_1d, label='Electron', linestyle = linestyles[0])
line0=ax0.plot(x, val2_1d, label=r'$\mathrm{D^+}$ ion', linestyle = linestyles[1])

xmin = 0
xmax = x.max()
ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((ymin, ymax))
ax0.set_ylabel(r"$N_\mathrm{M}\ (10^5)$", fontsize = label_fontsize)
#ax0.yaxis.set_label_coords(ylabel_x, 0.5)

#major_ticks = np.arange(0, 4.01, 2.0)
#ax0.set_yticks(major_ticks)


ax0.grid(True)
ax0.legend(loc = 4, framealpha=1)

ax0.annotate(r'$\mathbf{(a)}$', xy=get_axis_limits(ax0), annotation_clip=False)

##============Particle flux======================================================
val = collect("/Diagnostic/", "particle_flux_right")
val = val / pFlux0
val1_1d = val[1:, 0, 0, 0]
val2_1d = val[1:, 0, 0, 1]
print("ion particle flux min and max: ", val2_1d.min(), val2_1d.max())

ax0=fig.add_subplot(3,1,2)

ax0.yaxis.set_major_formatter(yformatter)

line0=ax0.plot(x, val1_1d, label='Electron', linestyle = linestyles[0])
line0=ax0.plot(x, val2_1d, label=r'$D^+ ion$', linestyle = linestyles[1])


ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))
ax0.set_ylabel(r"$\Gamma\ \mathrm{(10^{23}m^{-2}s^{-1})}$", fontsize = label_fontsize)
#ax0.yaxis.set_label_coords(ylabel_x, 0.5)

#major_ticks = np.arange(0, 4.01, 2.0)
#ax0.set_yticks(major_ticks)

ax0.grid(True)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============heat flux======================================================
val = collect("/Diagnostic/", "heat_flux_right")
val = val / hFlux0
val1_1d = val[1:, 0, 0, 0]
val2_1d = val[1:, 0, 0, 1]

ax0=fig.add_subplot(3,1,3)

ax0.yaxis.set_major_formatter(yformatter)

line0=ax0.plot(x, val1_1d, label='Electron', linestyle = linestyles[0])
line0=ax0.plot(x, val2_1d, label=r'$D^+ ion$', linestyle = linestyles[1])


ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((ymin, ymax))
ax0.set_xlabel(r"$t\ \mathrm{(\mu s)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$q\ \mathrm{(10^7Wm^{-2})}$", fontsize = label_fontsize)
#ax0.yaxis.set_label_coords(ylabel_x, 0.5)
ax0.grid(True)
#ax0.legend(loc = 1)

#major_ticks = np.arange(0, 1.51, 0.5)
#ax0.set_yticks(major_ticks)

ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)

fig.savefig("flux.png")
#fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
