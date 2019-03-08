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
Ex0 = 1.0e6

x_step = 20
level_num = 30

amplification_factor = 80.0


data_temp = collect("/Fields/", "Phi_global_avg")
nx = data_temp.shape[3]


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)
x = x * amplification_factor
x_less = x[::x_step]


xmin = x.min()
xmax = x.max()

x0 = int(nx / 2)
x1 = 100




##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

##============ potential ======================================================
ax0=fig.add_subplot(2,2,1)

val = collect("/Fields/", "Phi_global_avg")
val_2d = val[t, 0, :, :]

nx = val_2d.shape[0]
ny = val_2d.shape[1]
dx=1.0
dy=1.0
x, y=np.mgrid[slice(0.0,dx*nx,dx), slice(0.0,dy*ny,dy)]


contourf0 = ax0.contourf(x, y, val_2d, level_num, cmap=cm.get_cmap('jet'))

ax0.set_title(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
#ax0.axis('equal')
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlabel('x (mm)')
ax0.set_ylabel('y (mm)')
fig.colorbar(contourf0,)
ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ electric field in x direction ======================================================
ax0=fig.add_subplot(2,2,2)
val = collect("/Fields/", "Ex_global_avg")
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, level_num, cmap=cm.get_cmap('jet'))

ax0.set_title(r"$E \ \mathrm{(10^6V/m)}$", color='#1f77b4', fontsize = label_fontsize)
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlabel('x (mm)')
ax0.set_ylabel('y (mm)')
fig.colorbar(contourf0,)
ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ number density of e ======================================================
ax0=fig.add_subplot(2,2,3)

val = collect("/Fields/", "Rho_global_e_avg")
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, level_num, cmap=cm.get_cmap('jet'))

ax0.set_title(r"$n_e$", color='#1f77b4', fontsize = label_fontsize)
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlabel('x (mm)')
ax0.set_ylabel('y (mm)')
fig.colorbar(contourf0,)
ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ number density of D1 ======================================================
ax0=fig.add_subplot(2,2,4)
val = collect("/Fields/", "Rho_global_D1_avg")
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, level_num, cmap=cm.get_cmap('jet'))

ax0.set_title(r"$n_{D^{+}}$", color='#1f77b4', fontsize = label_fontsize)
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlabel('x (mm)')
ax0.set_ylabel('y (mm)')
fig.colorbar(contourf0,)
ax0.annotate(r"$\mathbf{(d)}$", xy=get_axis_limits(ax0), annotation_clip=False)



pdf_file_name = "Profiles2D_time" + str(t) + ".png"
fig.savefig(pdf_file_name, dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
