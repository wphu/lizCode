from  import *

if __name__ == "__main__":
    dx = 1.0e-5 #3.5e-6
    dy = 1.0e-5 #3.5e-6
    dz = 1.0e-5
    lx = 1.0e-3
    ly = 1.0e-3
    lz = 1.0e-3
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

    is_wall[:,:,0:nz_base+1] = 1
    
    # southwestern tile
    for i in np.arange(0, n_tile_half + 1):
        for j in np.arange(0, n_tile_half + 1):
            for k in np.arange(0, nz_wall_boundary + 1):
                is_wall[i,j,k] = 1

    # southeastern tile
    for i in np.arange(n_tile_half + n_gap_width, nx + 1):
        for j in np.arange(0, n_tile_half + 1):
            for k in np.arange(0, nz_wall_boundary+1):
                is_wall[i,j,k] = 1

    # northeastern tile
    for i in np.arange(n_tile_half + n_gap_width, nx + 1):
        for j in np.arange(n_tile_half + n_gap_width, ny + 1):
            for k in np.arange(0, nz_wall_boundary+1):
                is_wall[i,j,k] = 1

    # northwestern tile
    for i in np.arange(0, n_tile_half + 1):
        for j in np.arange(n_tile_half + n_gap_width, nx + 1):
            for k in np.arange(0, nz_wall_boundary+1):
                is_wall[i,j,k] = 1


    # ============================= bndr_type and bndr_val ========================
    # lower boundary of source region
    bndr_type[:, :, nz - nz_source] = 1
    bndr_val [:, :, nz - nz_source] = 0.0

    # source region, not solve
    bndr_type[:, :, nz - nz_source + 1 : nz + 1] = 5
    bndr_val [:, :, nz - nz_source + 1 : nz + 1] = 0.0

    # left and right boundary between source region and wall
    bndr_type[0, :, nz_wall_boundary + 1 : nz - nz_source] = 8
    bndr_type[nx,:, nz_wall_boundary + 1 : nz - nz_source] = 8
    bndr_type[:, 0, nz_wall_boundary + 1 : nz - nz_source] = 8
    bndr_type[:,ny, nz_wall_boundary + 1 : nz - nz_source] = 8

    # base
    bndr_type[0 : nx + 1, 0 : ny + 1, 0 : nz_base] = 5
    bndr_val [0 : nx + 1, 0 : ny + 1, 0 : nz_base] = wall_potential
    bndr_type[0 : nx + 1, 0 : ny + 1, nz_base] = 1
    bndr_val [0 : nx + 1, 0 : ny + 1, nz_base] = wall_potential

    # boundnary line of the wall
    ## west of south surface
    bndr_type[0 : n_tile_half + 1, 0, nz_wall_boundary] = 1
    bndr_val [0 : n_tile_half + 1, 0, nz_wall_boundary] = wall_potential
    bndr_type[n_tile_half, 0, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [n_tile_half, 0, nz_base + 1 : nz_wall_boundary + 1] = wall_potential
    ## east of south surface
    bndr_type[n_tile_half + n_gap_width : nx + 1, 0, nz_wall_boundary] = 1
    bndr_val [n_tile_half + n_gap_width : nx + 1, 0, nz_wall_boundary] = wall_potential
    bndr_type[n_tile_half + n_gap_width, 0, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [n_tile_half + n_gap_width, 0, nz_base + 1 : nz_wall_boundary + 1] = wall_potential
    ## south of east surface
    bndr_type[nx, 0 : n_tile_half + 1, nz_wall_boundary] = 1
    bndr_val [nx, 0 : n_tile_half + 1, nz_wall_boundary] = wall_potential
    bndr_type[nx, n_tile_half, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [nx, n_tile_half, nz_base + 1 : nz_wall_boundary + 1] = wall_potential
    ## north of east surface
    bndr_type[nx, n_tile_half + n_gap_width : ny + 1, nz_wall_boundary] = 1
    bndr_val [nx, n_tile_half + n_gap_width : ny + 1, nz_wall_boundary] = wall_potential
    bndr_type[nx, n_tile_half + n_gap_width, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [nx, n_tile_half + n_gap_width, nz_base + 1 : nz_wall_boundary + 1] = wall_potential
    ## east of north surface
    bndr_type[n_tile_half + n_gap_width : nx + 1, ny, nz_wall_boundary] = 1
    bndr_val [n_tile_half + n_gap_width : nx + 1, ny, nz_wall_boundary] = wall_potential
    bndr_type[n_tile_half, ny, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [n_tile_half, ny, nz_base + 1 : nz_wall_boundary + 1] = wall_potential
    ## west of north surface
    bndr_type[0 : n_tile_half + 1, ny, nz_wall_boundary] = 1
    bndr_val [0 : n_tile_half + 1, ny, nz_wall_boundary] = wall_potential
    bndr_type[n_tile_half + n_gap_width, ny, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [n_tile_half + n_gap_width, ny, nz_base + 1 : nz_wall_boundary + 1] = wall_potential
    ## north of west surface
    bndr_type[0, n_tile_half + n_gap_width : ny, nz_wall_boundary] = 1
    bndr_val [0, n_tile_half + n_gap_width : ny, nz_wall_boundary] = wall_potential
    bndr_type[0 , n_tile_half + n_gap_width, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [0 , n_tile_half + n_gap_width, nz_base + 1 : nz_wall_boundary + 1] = wall_potential  
    ## south of west surface
    bndr_type[0, 0 : n_tile_half + 1, nz_wall_boundary] = 1
    bndr_val [0, 0 : n_tile_half + 1, nz_wall_boundary] = wall_potential
    bndr_type[0 , n_tile_half, nz_base + 1 : nz_wall_boundary + 1] = 1
    bndr_val [0 , n_tile_half, nz_base + 1 : nz_wall_boundary + 1] = wall_potential  

    # boundary surface of the wall
    ## west of south surface
    bndr_type[0 : n_tile_half, 0, 0 : nz_wall_boundary] = 5
    bndr_val [0 : n_tile_half, 0, 0 : nz_wall_boundary] = wall_potential
    ## east of south surface
    bndr_type[n_tile_half + n_gap_width + 1 : nx + 1, 0, 0 : nz_wall_boundary] = 5
    bndr_val [n_tile_half + n_gap_width + 1 : nx + 1, 0, 0 : nz_wall_boundary] = wall_potential
    ## south of east surface
    bndr_type[nx, 0 : n_tile_half, 0 : nz_wall_boundary] = 5
    bndr_val [nx, 0 : n_tile_half, 0 : nz_wall_boundary] = wall_potential
    ## north of east surface
    bndr_type[nx, n_tile_half + n_gap_width + 1 : ny + 1, 0 : nz_wall_boundary] = 5
    bndr_val [nx, n_tile_half + n_gap_width + 1 : ny + 1, 0 : nz_wall_boundary] = wall_potential
    ## east of north surface
    bndr_type[n_tile_half + n_gap_width + 1 : nx + 1, ny, 0 : nz_wall_boundary] = 5
    bndr_val [n_tile_half + n_gap_width + 1 : nx + 1, ny, 0 : nz_wall_boundary] = wall_potential
    ## west of north surface
    bndr_type[0 : n_tile_half, ny, 0 : nz_wall_boundary] = 5
    bndr_val [0 : n_tile_half, ny, 0 : nz_wall_boundary] = wall_potential
    ## north of west surface
    bndr_type[0, n_tile_half + n_gap_width + 1 : ny + 1, 0 : nz_wall_boundary] = 5
    bndr_val [0, n_tile_half + n_gap_width + 1 : ny + 1, 0 : nz_wall_boundary] = wall_potential
    ## south of west surface
    bndr_type[0, 0 : n_tile_half, 0 : nz_wall_boundary] = 5
    bndr_val [0, 0 : n_tile_half, 0 : nz_wall_boundary] = wall_potential

    # the region between tiles of simulation domain surface
    ## south
    bndr_type[n_tile_half + 1 : n_tile_half + n_gap_width,  0, nz_base + 1 : nz_wall_boundary + 1] = 8  
    ## east
    bndr_type[nx, n_tile_half + 1 : n_tile_half + n_gap_width, nz_base + 1 : nz_wall_boundary + 1] = 8
    ## north
    bndr_type[n_tile_half + 1 : n_tile_half + n_gap_width, ny, nz_base + 1 : nz_wall_boundary + 1] = 8    
    ## west
    bndr_type[0 , n_tile_half + 1 : n_tile_half + n_gap_width, nz_base + 1 : nz_wall_boundary + 1] = 8


    # wall surface
    for i in np.arange(1, nx):
        for j in np.arange(1, ny):
            for k in np.arange(1, nz_wall_max):
                if is_wall[i, j, k] == 1 and (is_wall[i-1, j, k] == 0 or is_wall[i+1, j, k] == 0 or is_wall[i, j-1, k] == 0 or is_wall[i, j+1, k] == 0 or is_wall[i, j, k-1] == 0 or is_wall[i, j, k+1] == 0):
                    bndr_type[i, j, k] = 1
                    bndr_val [i, j, k] = wall_potential
                elif is_wall[i, j, k] == 1:
                    bndr_type[i, j, k] = 5
                    bndr_val [i, j, k] = wall_potential
    
    grid3d = Grid3D(dx, is_wall, bndr_type, bndr_val)
    grid3d.save_grid()
    grid3d.save_fig()