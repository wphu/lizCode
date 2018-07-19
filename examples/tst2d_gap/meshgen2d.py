from mesh2d import *

if __name__ == "__main__":
    dx = 3.5e-6
    dy = 3.5e-6
    lx = 3.5e-3
    ly = 5.2e-3
    nx = int(lx / dx)
    ny = int(ly / dy)
    ny_source = 20
    ny_base = 3

    print("nx, ny is : ", nx, ny)


    ny_wall_boundary = int(2.0e-3 / dy) + ny_base
    ny_wall_max = ny - ny_source - 5

    wall_potential = -60.0

    polygon_list = []

    point0 = Point(0.0,     ny_base*dy)
    point1 = Point(1.5e-3,  ny_base*dy)
    point2 = Point(1.5e-3,  ny_base*dy + 2.0e-3)
    point3 = Point(0.0,     ny_base*dy + 2.0e-3)

    line0 = Straight_line(point0, point1)
    line1 = Straight_line(point1, point2)
    line2 = Straight_line(point2, point3)
    line3 = Straight_line(point3, point0)

    polygon0 = Polygon()
    polygon0.add_straight_line(line0)
    polygon0.add_straight_line(line1)
    polygon0.add_straight_line(line2)
    polygon0.add_straight_line(line3)
    polygon_list.append(polygon0)


    point4 = Point(2.0e-3,  ny_base*dy)
    point5 = Point(3.5e-3,  ny_base*dy)
    point6 = Point(3.5e-3,  ny_base*dy + 2.0e-3)
    point7 = Point(2.5e-3,  ny_base*dy + 2.0e-3)
    point8 = Point(2.0e-3,  ny_base*dy + 1.5e-3)

    line4 = Straight_line(point4, point5)
    line5 = Straight_line(point5, point6)
    line6 = Straight_line(point6, point7)
    line7 = Straight_line(point7, point8)
    line8 = Straight_line(point8, point4)

    polygon1 = Polygon()
    polygon1.add_straight_line(line4)
    polygon1.add_straight_line(line5)
    polygon1.add_straight_line(line6)
    polygon1.add_straight_line(line7)
    polygon1.add_straight_line(line8)
    polygon_list.append(polygon1)

    is_wall     = np.zeros((nx+1, ny+1), dtype = 'int')
    bndr_type   = np.zeros((nx+1, ny+1), dtype = 'int')
    bndr_val    = np.zeros((nx+1, ny+1), dtype = 'int')

    # ============================= is_wall ========================
    is_wall[0,:] = 1
    is_wall[nx,:] = 1
    is_wall[:,0] = 1
    is_wall[:,ny] = 1

    is_wall[:, 0:ny_base+1] = 1
    
    for i in np.arange(0, nx+1):
        for j in np.arange(0, ny+1):
            point_temp = Point(i*dx, j*dy)
            is_in_polygon = False
            for polygon_temp in polygon_list:
                if polygon_temp.is_in_polygon(point_temp):
                    is_in_polygon = True
            if is_in_polygon:
                is_wall[i,j] = 1
                
    # ============================= bndr_type and bndr_val ========================
    # lower boundary of source region
    bndr_type[:, ny-ny_source] = 1
    bndr_val [:, ny-ny_source] = 0.0

    # source region, not solve
    bndr_type[:, ny-ny_source+1:ny+1] = 5
    bndr_val [:, ny-ny_source+1:ny+1] = 0.0

    # left and right boundary between source region and wall
    bndr_type[0, ny_wall_boundary+1:ny-ny_source] = 8
    bndr_type[nx,ny_wall_boundary+1:ny-ny_source] = 8

    # corner points of wall at left and right boudnary
    bndr_type[0, ny_wall_boundary] = 1
    bndr_type[nx,ny_wall_boundary] = 1
    bndr_val [0, ny_wall_boundary] = wall_potential
    bndr_val [nx,ny_wall_boundary] = wall_potential

    # left, right and bottom boudnary of the wall
    bndr_type[0, 0:ny_wall_boundary] = 5
    bndr_type[nx,0:ny_wall_boundary] = 5
    bndr_type[0:nx,0] 		       = 5
    bndr_val [0, 0:ny_wall_boundary] = wall_potential
    bndr_val [nx,0:ny_wall_boundary] = wall_potential
    bndr_val [0:nx,0] 		       = wall_potential

    # wall surface
    for i in np.arange(1, nx):
        for j in np.arange(1, ny_wall_max):
            if is_wall[i,j] == 1 and (is_wall[i-1, j] == 0 or is_wall[i+1, j] == 0 or is_wall[i, j-1] == 0 or is_wall[i, j+1] == 0):
                bndr_type[i,j] = 1
                bndr_val [i,j] = wall_potential
            elif is_wall[i,j] == 1:
                bndr_type[i,j] = 5
                bndr_val [i,j] = wall_potential
    
    # ============================= boundary lines ========================
    bndr_line_list = []
    bndr_line0 = Straight_line(point3, point2)
    bndr_line1 = Straight_line(point2, point1)
    bndr_line2 = Straight_line(point1, point4)
    bndr_line3 = Straight_line(point4, point8)
    bndr_line4 = Straight_line(point8, point7)
    bndr_line5 = Straight_line(point7, point6)

    bndr_line_list.append(bndr_line0)
    bndr_line_list.append(bndr_line1)
    bndr_line_list.append(bndr_line2)
    bndr_line_list.append(bndr_line3)
    bndr_line_list.append(bndr_line4)
    bndr_line_list.append(bndr_line5)

    segment_list = []
    n_segments = []
    for bndr_line_temp in bndr_line_list:
        segment_list_temp = bndr_line_temp.generate_segments(dx, dy)
        n_segments.append(len(segment_list_temp))
        for segment_temp in segment_list_temp:
            segment_list.append(segment_temp)

    grid2d = Grid2D(dx, dy, is_wall, bndr_type, bndr_val, n_segments, segment_list)
    grid2d.save_grid()
    grid2d.save_fig()

