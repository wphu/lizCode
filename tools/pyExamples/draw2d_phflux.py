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



length = collect("length", prefix = "grid")
length = length[0,:]

x = np.zeros(length.shape[0])
x[0] = 0.5 * length[0]
for i in np.arange(1, length.shape[0]):
	x[i] = x[i-1] + 0.5 * (length[i-1] + length[i])
x = x * 1.0e3
xmin = 0.0
xmax = x.max()

particle_flux_normalize = 1.0e23
heat_flux_normalize = 1.0e6


##inite the fig of matplotlib
fig=plt.figure(figsize=(10,6))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

##============ particle flux======================================================
ax0=fig.add_subplot(2,1,1)

val = collect("/Diagnostic/", "particleFlux")
val_2d = val[t, 0, :, :] / particle_flux_normalize

line0=ax0.plot(x, val_2d[0], label = r"$\mathrm{e}$")
line0=ax0.plot(x, val_2d[1], label = r"$\mathrm{D^+}$")

print("e particle flux: ",val_2d[0,0], val_2d[0,20], val_2d[0].max(), val_2d[0].max() / val_2d[0,20])
print("D1 particle flux: ",val_2d[1,0], val_2d[1,20], val_2d[1].max(), val_2d[1].max() / val_2d[1,20])

ax0.grid(True)

ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 4.5))
y_major_ticks = np.arange(0, 4.5, 2.0)
ax0.set_yticks(y_major_ticks)
ax0.set_ylabel(r"$\Gamma\ \mathrm{(10^{23}m^{-2}s^{-1})}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ heat flux ======================================================
ax0=fig.add_subplot(2,1,2)

val = collect("/Diagnostic/", "heatFlux")
val_2d = val[t, 0, :, :] / heat_flux_normalize

line0=ax0.plot(x, val_2d[0], label = r"$\mathrm{e}$")
line0=ax0.plot(x, val_2d[1], label = r"$\mathrm{D^+}$")

print("e heat flux: ",val_2d[0,0], val_2d[0,20], val_2d[0].max(), val_2d[0].max() / val_2d[0,20])
print("D1 heat flux: ",val_2d[1,0], val_2d[1,20], val_2d[1].max(), val_2d[1].max() / val_2d[1,20])

ax0.grid(True)

ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 3.1))
y_major_ticks = np.arange(0, 3.1, 1.0)
ax0.set_yticks(y_major_ticks)
ax0.set_ylabel(r"$q\ \mathrm{(MWm^{-2})}$", fontsize = label_fontsize)

ax0.set_xlabel(r"$\mathrm{x\ (mm)}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)



pdf_file_name = "particle_heat_flux_time" + str(t) + ".png"
fig.savefig(pdf_file_name, dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
