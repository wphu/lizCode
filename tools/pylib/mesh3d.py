from template import *
import math
from scipy import constants as const


class Grid3D:
    def __init__(self, dx, is_wall, bndr_type, bndr_val):
        self.n_dim = 3
        self.dx = dx
        self.is_wall = is_wall
        self.bndr_type = bndr_type
        self.bndr_val = bndr_val


    def save_grid(self):
        f = h5.File('data/grid.h5', 'w')
        f.attrs['n_dim'] = self.n_dim
        f['is_wall']     = self.is_wall
        f['bndr_type']   = self.bndr_type
        f['bndr_val']    = self.bndr_val
        f.close()

    def save_fig(self):
        level_num = 50

        ##inite the fig of matplotlib
        fig=plt.figure(figsize=(10,8))
        fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

        ##============ south ======================================================
        ax0=fig.add_subplot(2,2,1)

        data = self.bndr_type[:,0,:]
        nx = data.shape[0]
        ny = data.shape[1]
        dx=1.0
        dy=1.0
        x, y=np.mgrid[slice(0.0,dx*nx,dx), slice(0.0,dy*ny,dy)]

        contourf0 = ax0.contourf(x, y, data, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{south}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

        ##============ east ======================================================
        ax0=fig.add_subplot(2,2,2)

        data = self.bndr_type[-1,:,:]
        nx = data.shape[0]
        ny = data.shape[1]
        dx=1.0
        dy=1.0
        x, y=np.mgrid[slice(0.0,dx*nx,dx), slice(0.0,dy*ny,dy)]

        contourf0 = ax0.contourf(x, y, data, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{south}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

        ##============ north ======================================================
        ax0=fig.add_subplot(2,2,3)

        data = self.bndr_type[:,-1,:]
        nx = data.shape[0]
        ny = data.shape[1]
        dx=1.0
        dy=1.0
        x, y=np.mgrid[slice(0.0,dx*nx,dx), slice(0.0,dy*ny,dy)]

        contourf0 = ax0.contourf(x, y, data, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{south}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)

        ##============ west ======================================================
        ax0=fig.add_subplot(2,2,4)

        data = self.bndr_type[0,:,:]
        nx = data.shape[0]
        ny = data.shape[1]
        dx=1.0
        dy=1.0
        x, y=np.mgrid[slice(0.0,dx*nx,dx), slice(0.0,dy*ny,dy)]

        contourf0 = ax0.contourf(x, y, data, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{south}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(d)}$", xy=get_axis_limits(ax0), annotation_clip=False)



        pdf_file_name = "grid" + ".png"
        fig.savefig(pdf_file_name, dpi = 300)
