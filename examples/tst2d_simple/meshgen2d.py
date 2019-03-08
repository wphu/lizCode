from mesh2d import *


if __name__ == "__main__":
    dx = 1.0e-5 #3.5e-6
    dy = 1.0e-5 #3.5e-6
    nx = 100
    ny = 200
    lx = dx * nx
    ly = dy * ny


    print("nx, ny is : ", nx, ny)
    print("lx, ly is : ", lx, ly)

    polygon_list = []

    point0 = Point(0.0, 0.0)
    point1 = Point(lx,  0.0)
    point2 = Point(lx,  ly)
    point3 = Point(0.0, ly)
    
    is_wall     = np.zeros((nx+1, ny+1), dtype = 'int')
    bndr_type   = np.zeros((nx+1, ny+1), dtype = 'int')
    bndr_val    = np.zeros((nx+1, ny+1), dtype = 'int')

    # ============================= is_wall ========================
    is_wall[0,:] = 1
    is_wall[nx,:] = 1
    is_wall[:,0] = 1
    is_wall[:,ny] = 1

           
    # ============================= bndr_type and bndr_val ========================
    # south and north boundary
    bndr_type[0 : nx+1, 0] = 1
    bndr_val [0 : nx+1, 0] = 0.0
    bndr_type[0 : nx+1, ny] = 1
    bndr_val [0 : nx+1, ny] = 0.0

    # west and east boundary
    bndr_type[0, 1 : ny] = 8
    bndr_type[nx, 1 : ny] = 8
    
    # ============================= boundary lines ========================
    bndr_line_list = []
    bndr_line0 = Straight_line(point0, point1)
    bndr_line1 = Straight_line(point3, point2)


    bndr_line_list.append(bndr_line0)
    bndr_line_list.append(bndr_line1)


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
