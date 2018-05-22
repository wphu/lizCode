from template import *
from collect import collect

import sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

mpl.rcParams['lines.linewidth'] = 1.0


itime = 9
if len(sys.argv) == 2:
	t = int(sys.argv[1])
elif len(sys.argv) > 2:
	print("error: should have one argument, as time step")
else:
	t = itime



angle = 5.0 * math.pi / 180.0
bevel_depth_unit =   0.5 * math.tan(angle)


x_step = 20
level_num = 30

path = "ref/data"
val = collect("/Fields/", "Phi_global_avg", path = path)
val_2d = val[t, 0, :, :]

nx = val_2d.shape[0]
ny = val_2d.shape[1]

dx = 0.35e-2
dy = 0.35e-2

x, y=np.mgrid[slice(0.0,dx*(nx-0.5),dx), slice(0.0,dy*ny,dy)]
#x, y=np.meshgrid(np.arange(0.0,dx*nx,dx), np.arange(0.0,dy*ny,dy))
print("shape of x, y: ", x.shape, y.shape)
print("nx, ny: ", nx, ny)

xmin = 1.0 #(0.5 * nx - 200) * dx
xmax = 2.5 #(0.5 * nx + 200) * dx

ymin = 1.0 #200 * dy
ymax = 3.0 #800 * dy

x_major_ticks = np.arange(1.0, 2.502, 0.5)

zmin = -60.0
zmax = 0.0
levels = np.linspace(zmin, zmax, level_num)


##inite the fig of matplotlib
fig=plt.figure(figsize=(8.5,8))
fig.subplots_adjust(wspace=0.03,hspace=0.35)

fig.suptitle(r"Electric potential $\mathrm{(V)}$", fontsize=18, fontweight='bold')


# axes(position) for colorbar
cax = fig.add_axes([0.85, 0.2, 0.015, 0.6])

#ticks for colorbar
ticks_colorbar = np.linspace(zmin, zmax, 5)


##============ case 0 ======================================================
ax0=fig.add_subplot(2,2,1)

path = "ref/data"
val = collect("/Fields/", "Phi_global_avg", path = path)
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, levels = levels, cmap=cm.get_cmap('jet'))

#ax0.set_title(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
#ax0.axis('equal')
ax0.set_aspect('equal', adjustable='box')
#ax0.set_xlabel('x (mm)')
#ax0.set_ylabel('y (mm)')
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))
ax0.set_xticks(x_major_ticks)

xy_polygon = np.zeros((4, 2))
xy_polygon[0, 0] = 0.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 1.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 1.5
xy_polygon[2, 1] = 2.0105
xy_polygon[3, 0] = 0.0
xy_polygon[3, 1] = 2.0105
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

xy_polygon[0, 0] = 2.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 3.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 3.5
xy_polygon[2, 1] = 2.0105
xy_polygon[3, 0] = 2.0
xy_polygon[3, 1] = 2.0105
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ case 1 ======================================================
ax0=fig.add_subplot(2,2,2)

path = "ref_bevel1/data"
val = collect("/Fields/", "Phi_global_avg", path = path)
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, levels = levels, cmap=cm.get_cmap('jet'))

#ax0.set_title(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
#ax0.axis('equal')
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))
ax0.set_xticks(x_major_ticks)
fig.colorbar(contourf0, cax, ticks = ticks_colorbar)

xy_polygon = np.zeros((4, 2))
xy_polygon[0, 0] = 0.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 1.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 1.5
xy_polygon[2, 1] = 2.0105 + 1.0 * bevel_depth_unit 
xy_polygon[3, 0] = 0.0
xy_polygon[3, 1] = 2.0105 + 0.5 * 1.0 * bevel_depth_unit
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

xy_polygon[0, 0] = 2.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 3.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 3.5
xy_polygon[2, 1] = 2.0105 + 0.5 * 1.0 * bevel_depth_unit
xy_polygon[3, 0] = 2.0
xy_polygon[3, 1] = 2.0105 
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ case 2 ======================================================
ax0=fig.add_subplot(2,2,3)

path = "ref_bevel2/data"
val = collect("/Fields/", "Phi_global_avg", path = path)
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, levels = levels, cmap=cm.get_cmap('jet'))

#ax0.set_title(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
#ax0.axis('equal')
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))
ax0.set_xticks(x_major_ticks)

xy_polygon = np.zeros((4, 2))
xy_polygon[0, 0] = 0.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 1.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 1.5
xy_polygon[2, 1] = 2.0105 + 2.0 * bevel_depth_unit 
xy_polygon[3, 0] = 0.0
xy_polygon[3, 1] = 2.0105 + 0.5 * 2.0 * bevel_depth_unit
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

xy_polygon[0, 0] = 2.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 3.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 3.5
xy_polygon[2, 1] = 2.0105 + 0.5 * 2.0 * bevel_depth_unit
xy_polygon[3, 0] = 2.0
xy_polygon[3, 1] = 2.0105 
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ case 3 ======================================================
ax0=fig.add_subplot(2,2,4)

path = "ref_bevel3/data"
val = collect("/Fields/", "Phi_global_avg", path = path)
val_2d = val[t, 0, :, :]

contourf0 = ax0.contourf(x, y, val_2d, levels = levels, cmap=cm.get_cmap('jet'))

#ax0.set_title(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
#ax0.axis('equal')
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))
ax0.set_xticks(x_major_ticks)

xy_polygon = np.zeros((4, 2))
xy_polygon[0, 0] = 0.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 1.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 1.5
xy_polygon[2, 1] = 2.0105 + 3.0 * bevel_depth_unit 
xy_polygon[3, 0] = 0.0
xy_polygon[3, 1] = 2.0105 + 0.5 * 3.0 * bevel_depth_unit
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

xy_polygon[0, 0] = 2.0
xy_polygon[0, 1] = 0.0
xy_polygon[1, 0] = 3.5
xy_polygon[1, 1] = 0.0
xy_polygon[2, 0] = 3.5
xy_polygon[2, 1] = 2.0105 + 0.5 * 3.0 * bevel_depth_unit
xy_polygon[3, 0] = 2.0
xy_polygon[3, 1] = 2.0105 
polygon = Polygon(xy_polygon, True, color="#aaaaaa")
ax0.add_patch(polygon)

ax0.annotate(r"$\mathbf{(d)}$", xy=get_axis_limits(ax0), annotation_clip=False)


fig.text(0.5, 0.04, 'x (mm)', ha='center', va='center',size=16)
fig.text(0.1, 0.5, 'y (mm)', ha='center', va='center', rotation='vertical',size=16)





pdf_file_name = "all_potential_time" + str(t) + ".png"
fig.savefig(pdf_file_name, dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
