from template import *
import math
from scipy import constants as const


class Grid2D:
    def __init__(self, dx, is_wall, bndr_type, bndr_val):
        self.dx = dx
        self.is_wall = is_wall
        self.bndr_type = bndr_type
        self.bndr_val = bndr_val


    def save_grid(self):
        f = h5.File('data/grid.h5', 'w')
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


if __name__ == "__main__":
    dx = 1.0e-5 #3.5e-6
    dy = 1.0e-5 #3.5e-6
    dz = 1.0e-5
    lx = 0.5e-3
    ly = 0.5e-3
    lz = 0.5e-3
    nx = int(lx / dx)
    ny = int(ly / dy)
    nz = int(lz / dz)
    nz_source = 20
    nz_base = 3
    n_gap_width = int(0.5e-3 / dx)
    n_tile_half = int(0.25e-3 / dx)

    print("nx, ny, nz is : ", nx, ny, nz)


    nz_wall_boundary = int(0.5e-3 / dz) + nz_base
    nz_wall_max = nz - nz_source - 5

    wall_potential = -60.0

    is_wall     = np.zeros((nx+1, ny+1, nz+1), dtype = 'int')
    bndr_type   = np.zeros((nx+1, ny+1, nz+1), dtype = 'int')
    bndr_val    = np.zeros((nx+1, ny+1, nz+1), dtype = 'double')

    # ============================= is_wall ========================
    is_wall[0,:,:] = 1
    is_wall[nx,:,:] = 1
    is_wall[:,0,:] = 1
    is_wall[:,ny,:] = 1
    is_wall[:,:,0] = 1
    is_wall[:,:,nz] = 1

    # ============================= bndr_type and bndr_val ========================
    # bottom
    bndr_type[:, :, 0] = 1
    bndr_val [:, :, 0] = -10.0

    # up
    bndr_type[:, :, nz] = 1
    bndr_val [:, :, nz] = 0.0

    '''	
    # west
    bndr_type[0, :, :] = 1
    bndr_val [0, :, :] = 0.0

    # east
    bndr_type[nx, :, :] = 1
    bndr_val [nx, :, :] = 0.0

    # south
    bndr_type[:, 0, :] = 1
    bndr_val [:, 0, :] = 0.0

    # north
    bndr_type[:, ny, :] = 1
    bndr_val [:, ny, :] = 0.0
    '''	


    
    # west
    bndr_type[0, :, 1:nz] = 8
    bndr_val [0, :, 1:nz] = 0.0

    # east
    bndr_type[nx, :, 1:nz] = 8
    bndr_val [nx, :, 1:nz] = 0.0

    # south
    bndr_type[:, 0, 1:nz] = 8
    bndr_val [:, 0, 1:nz] = 0.0

    # north
    bndr_type[:, ny, 1:nz] = 8
    bndr_val [:, ny, 1:nz] = 0.0
    

    grid2d = Grid2D(dx, is_wall, bndr_type, bndr_val)
    grid2d.save_grid()
    grid2d.save_fig()