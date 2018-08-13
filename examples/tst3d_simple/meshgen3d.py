from mesh3d import *

if __name__ == "__main__":
    dx = 1.0e-5 #3.5e-6
    dy = 1.0e-5 #3.5e-6
    dz = 1.0e-5
    lx = 0.2e-3
    ly = 0.2e-3
    lz = 0.2e-3
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

    print("gap height is : ", nz_wall_boundary)

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
    bndr_val [:, :, 0] = 0.0

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
    

    grid3d = Grid3D(dx, is_wall, bndr_type, bndr_val)
    grid3d.save_grid()
    grid3d.save_fig()