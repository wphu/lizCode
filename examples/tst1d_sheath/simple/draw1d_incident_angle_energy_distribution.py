from template import *
from collect import collect

import sys

if len(sys.argv) == 2:
	t = int(sys.argv[1])
elif len(sys.argv) > 2:
	print("error: should have one argument, as time step")
else:
	t = 0

print("time: ", t)

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.6,hspace=0.55)


#============angle distribution=======================================
val = collect("/Diagnostic/", "angle_distribution_right", itime = t)
val1_1d = val[0, 0, :]
val2_1d = val[0, 1, :]
print(val.shape)

nx = (val1_1d.shape)[0]
x = np.linspace(0,nx,nx)


ax0 = fig.add_subplot(2,1,1)

ax0.yaxis.set_major_formatter(yformatter)

line0=ax0.plot(x, val1_1d, label='Electron', linestyle = linestyles[0])
line0=ax0.plot(x, val2_1d, label=r'$\mathrm{D^+}$ ion', linestyle = linestyles[1])

xmin = 0
xmax = x.max()
ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((ymin, ymax))
ax0.set_ylabel(r"$angle\ distribution$", fontsize = label_fontsize)
#ax0.yaxis.set_label_coords(ylabel_x, 0.5)

#major_ticks = np.arange(0, 4.01, 2.0)
#ax0.set_yticks(major_ticks)


ax0.grid(True)
ax0.legend(loc = 4, framealpha=1)

ax0.annotate(r'$\mathbf{(a)}$', xy=get_axis_limits(ax0), annotation_clip=False)

##============Particle flux======================================================
val = collect("/Diagnostic/", "energy_distribution_right", itime = t)
val1_1d = val[0, 0, :]
val2_1d = val[0, 1, :]
print(val.shape)

nx = (val1_1d.shape)[0]
x = np.linspace(0,nx,nx)

ax0=fig.add_subplot(2,1,2)

ax0.yaxis.set_major_formatter(yformatter)

line0=ax0.plot(x, val1_1d, label='Electron', linestyle = linestyles[0])
line0=ax0.plot(x, val2_1d, label=r'$D^+ ion$', linestyle = linestyles[1])


ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

ax0.set_xlim((xmin, xmax))
ax0.set_ylabel(r"$energy\ distribution$", fontsize = label_fontsize)
#ax0.yaxis.set_label_coords(ylabel_x, 0.5)

#major_ticks = np.arange(0, 4.01, 2.0)
#ax0.set_yticks(major_ticks)

ax0.grid(True)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)



pdf_file_name = "incident_angle_energy_distribution" + str(t) + ".png"
fig.savefig(pdf_file_name, dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
