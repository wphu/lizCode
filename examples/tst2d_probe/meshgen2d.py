from mesh2d import *


if __name__ == "__main__":
    dx = 0.5e-5 #3.5e-6
    dy = 0.5e-5 #3.5e-6
    lx = 24.0e-3
    ly = 5.2e-3
    nx = int(lx / dx)
    ny = int(ly / dy)
    ny_source = 20
    ny_base = 3

    x1 = 20.0e-3
    x2 = 2.0e-3
    x3 = 10.0e-3
    x4 = 2.0e-3
    y1 = 2.0e-3
    y2 = 1.0e-3


    print("nx, ny is : ", nx, ny)


    ny_wall_boundary = int(y1 / dy) + ny_base
    ny_wall_max = ny - ny_source - 5

    wall_potential = -60.0

    polygon_list = []

    # left tile
    point0 = Point(0.0,     ny_base*dy)
    point1 = Point(x2,      ny_base*dy)
    point2 = Point(x2,      ny_base*dy + y1)
    point3 = Point(0.0,     ny_base*dy + y1)

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

    # right tile
    point4 = Point(x1 + x2, ny_base*dy)
    point5 = Point(lx,      ny_base*dy)
    point6 = Point(lx,      ny_base*dy + y1)
    point7 = Point(x1 + x2, ny_base*dy + y1)

    line4 = Straight_line(point4, point5)
    line5 = Straight_line(point5, point6)
    line6 = Straight_line(point6, point7)
    line7 = Straight_line(point7, point4)

    polygon1 = Polygon()
    polygon1.add_straight_line(line4)
    polygon1.add_straight_line(line5)
    polygon1.add_straight_line(line6)
    polygon1.add_straight_line(line7)
    polygon_list.append(polygon1)


    # the probe in the gap
    point8 = Point(x2 + x3,         ny_base*dy)
    point9 = Point(x2 + x3 + x4,    ny_base*dy)
    point10 = Point(x2 + x3 + x4,   ny_base*dy + y2)
    point11 = Point(x2 + x3,        ny_base*dy + y2)

    line8 = Straight_line(point8, point9)
    line9 = Straight_line(point9, point10)
    line10 = Straight_line(point10, point11)
    line11 = Straight_line(point11, point8)

    polygon1 = Polygon()
    polygon1.add_straight_line(line8)
    polygon1.add_straight_line(line9)
    polygon1.add_straight_line(line10)
    polygon1.add_straight_line(line11)
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
    bndr_line2 = Straight_line(point1, point8)
    bndr_line3 = Straight_line(point8, point11)
    bndr_line4 = Straight_line(point11, point10)
    bndr_line5 = Straight_line(point10, point9)
    bndr_line6 = Straight_line(point9, point4)
    bndr_line7 = Straight_line(point4, point7)
    bndr_line8 = Straight_line(point7, point6)

    bndr_line_list.append(bndr_line0)
    bndr_line_list.append(bndr_line1)
    bndr_line_list.append(bndr_line2)
    bndr_line_list.append(bndr_line3)
    bndr_line_list.append(bndr_line4)
    bndr_line_list.append(bndr_line5)
    bndr_line_list.append(bndr_line6)
    bndr_line_list.append(bndr_line7)
    bndr_line_list.append(bndr_line8)

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
