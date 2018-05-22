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




xmin = -1.0
xmax = 1.0

particle_flux_normalize = 1.0e23
heat_flux_normalize = 1.0e6


##inite the fig of matplotlib
fig=plt.figure(figsize=(12,6))
fig.subplots_adjust(left = 0.1, right = 0.8, bottom=0.15,wspace=0.2,hspace=0.45)


# heat_flux_list0: for electron; heat_flux_list1: for D1 ion
# _left: left tile; _right: right tile
x_left_list  = []
x_right_list = []
heat_flux_list0_left  = []
heat_flux_list0_right = []
heat_flux_list1_left  = []
heat_flux_list1_right = []


path_list = ["ref/data", "ref_bevel1/data", "ref_bevel2/data", "ref_bevel3/data", "ref_bevel4/data"]
labels    = [r"$\mathrm{h}\ =\ 0$", r"$\mathrm{h}\ =\ 1$", r"$\mathrm{h}\ =\ 2$", r"$\mathrm{h}\ =\ 3$", r"$\mathrm{h}\ =\ 4$"]

for path in path_list:
	# read segment information
	length = collect("length", prefix = "grid", path = path)
	length = length[0,:]
	length = length * 1.0e3

	n_segments = collect("n_segments", prefix = "grid", path = path)
	n_segments = n_segments[0,:]
	print(n_segments.shape, type(n_segments[0]), n_segments[0])

	nx_start_left 		= 0
	nx_end_left	  		= n_segments[0] + n_segments[1]
	nx_start_right 		= nx_end_left + n_segments[2]
	nx_end_right	  	= nx_start_right + n_segments[3] + n_segments[4]

	x0_left  = 0.0
	x0_right = 0.0

	x = np.zeros(length.shape[0])
	x[0] = 0.5 * length[0]
	for i in np.arange(1, length.shape[0]):
		x[i] = x[i-1] + 0.5 * (length[i-1] + length[i])
		if i == n_segments[0]:
			x0_left = x[i]
		if i == nx_end_right - n_segments[4]:
			x0_right = x[i]
			
	x_left  = x[nx_start_left:nx_end_left]   - x0_left
	x_right = x[nx_start_right:nx_end_right] - x0_right
	x_left_list.append(x_left)
	x_right_list.append(x_right)
	
	
	# read data
	val = collect("/Diagnostic/", "heatFlux", path = path)
	val_2d = val[t, 0, :, :] / heat_flux_normalize

	val_1d = val_2d[0]
	val_1d_left  = val_1d[nx_start_left:nx_end_left]
	val_1d_right = val_1d[nx_start_right:nx_end_right]
	heat_flux_list0_left.append(val_1d_left)
	heat_flux_list0_right.append(val_1d_right)

	val_1d = val_2d[1]
	val_1d_left  = val_1d[nx_start_left:nx_end_left]
	val_1d_right = val_1d[nx_start_right:nx_end_right]
	heat_flux_list1_left.append(val_1d_left)
	heat_flux_list1_right.append(val_1d_right)


##============ heat flux for e at left ======================================================
ax0=fig.add_subplot(2,2,1)

for i in np.arange(0,len(heat_flux_list0_left)):
	line0 = ax0.plot(x_left_list[i], heat_flux_list0_left[i], label = r"$e$")

ax0.grid(True)

#ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 3.01))
y_major_ticks = np.arange(0, 3.1, 1.0)
ax0.set_yticks(y_major_ticks)
ax0.set_ylabel(r"$q_{\mathrm{e}}\ \mathrm{(MWm^{-2})}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(a1)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ heat flux for e at right ======================================================
ax0=fig.add_subplot(2,2,2)

for i in np.arange(0,len(heat_flux_list0_right)):
	line0 = ax0.plot(x_right_list[i], heat_flux_list0_right[i], label = labels[i])

ax0.grid(True)

#ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 3.01))
y_major_ticks = np.arange(0, 3.1, 1.0)
ax0.set_yticks(y_major_ticks)

legend0 = ax0.legend(bbox_to_anchor=(1.06, 0.3, 0.1, 0.8), loc=3, ncol=1, borderaxespad=0., fontsize = 11, fancybox=True)
legend0.get_frame().set_linewidth(0.5)




ax0.annotate(r"$\mathbf{(a2)}$", xy=get_axis_limits(ax0), annotation_clip=False)


##============ heat flux for D1 ion at left ======================================================
ax0=fig.add_subplot(2,2,3)

for i in np.arange(0,len(heat_flux_list1_left)):
	line0 = ax0.plot(x_left_list[i], heat_flux_list1_left[i], label = r"$e$")

ax0.grid(True)

#ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 3.01))
y_major_ticks = np.arange(0, 3.1, 1.0)
ax0.set_yticks(y_major_ticks)
ax0.set_ylabel(r"$q_{\mathrm{D^+}}\ \mathrm{(MWm^{-2})}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(b1)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ heat flux for D1 ion at right ======================================================
ax0=fig.add_subplot(2,2,4)

for i in np.arange(0,len(heat_flux_list1_right)):
	line0 = ax0.plot(x_right_list[i], heat_flux_list1_right[i], label = r"$e$")

ax0.grid(True)

#ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 3.01))
y_major_ticks = np.arange(0, 3.1, 1.0)
ax0.set_yticks(y_major_ticks)

ax0.annotate(r"$\mathbf{(b2)}$", xy=get_axis_limits(ax0), annotation_clip=False)

fig.text(0.45, 0.05, r"$\mathrm{x\ (mm)}$", ha='center', va='center', size = label_fontsize)
#fig.text(0.05, 0.5, r"$q\ \mathrm{(MWm^{-2})}$", ha='center', va='center', rotation='vertical',size=16)

pdf_file_name = "all_heat_flux_time" + str(t) + ".png"
fig.savefig(pdf_file_name, dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
