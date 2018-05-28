from template import *
from collect import collect

import sys

itime = 7
if len(sys.argv) == 2:
	t = int(sys.argv[1])
elif len(sys.argv) > 2:
	print("error: should have one argument, as time step")
else:
	t = itime

n0 = 1.0e19
Ex0 = 1.0e6

x_step = 20

amplification_factor = 80.0


data_temp = collect("/Diagnostic/", "radiative_energy_collision")
nx = data_temp.shape[3]


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)
x = x * amplification_factor
x_less = x[::x_step]


xmin = x.min()
xmax = x.max() * 0.02

x0 = int(nx / 2)
x1 = 100




##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

##============rho======================================================
ax0=fig.add_subplot(1,1,1)

val = collect("/Diagnostic/", "radiative_energy_collision")
val = val[:]
val_2d = val[t,0,:,:]
n_collide = val_2d.shape[0]
for i_collide in np.arange(0, n_collide):
	#val_1d = np.transpose(val[t, 0, 0, :])
	#print( "Electron density: ",val_1d[x0], val_1d[x1] )
	val_1d = val_2d[i_collide,:]
	line0=ax0.plot(x, val_1d, label = i_collide)



ax0.grid(True)





ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((0.0, 5.1))
major_ticks = np.arange(0, 5.1, 1.5)
ax0.set_yticks(major_ticks)
ax0.set_ylabel(r"$n\ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)



pdf_file_name = "radiative_energy_collision" + str(t) + ".pdf"
fig.savefig(pdf_file_name, dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
