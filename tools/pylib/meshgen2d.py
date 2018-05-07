from template import *

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z - z

class Segment:
    def __init__(self):
        self.start_point = None
        self.end_point = None
        self.length = None
        self.grid_point = None
        self.normal = None

    def cal_length():
        x_dif = self.end_point.x - self.start_point.x
        y_dif = self.end_point.y - self.start_point.y
        self.length = math.sqrt(math.pow(x_dif, 2) + math.pow(y_dif, 2))

class Straight_line:
    def __init__(self, start_point, end_point):
        self.start_point = start_point
        self.end_point = end_point
    def generate_segments(self, dx, dy):
        self.segment_list = []
        # calculate a and b for slope-intercept form of a line: y = ax + b
        a = 
        b = 

        i_min = int(self.start_point.x / dx)
        i_max = int(self.end_point.x / dx)
        j_min = int(self.start_point.y / dy)
        j_max = int(self.end_point.y / dy)
        if i_min > i_max:
            i_temp = i_min
            i_min = i_max
            i_max = i_temp
        if j_min > j_max:
            j_temp = j_min
            j_min = j_max
            j_max = j_temp
        
        for  i in np.arange(i_min, imax + 1):
            y0 = a * i * dx + b
            y1 = a * (i+1) * dx + b
            j0 = int(y_lower / dy)
            j1 = int(y_upper / dy)
            if j0 == j1:
                segment_temp = Segment()
                segment_temp.start_point = Point(i * dx, y0)
                segment_temp.end_point   = Point((i + 1) * dx, y1)
                segment_temp.grid_point  = Point(i, j0)
                segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                segment_temp.cal_length()
                self.segment_list.append(segment_temp)
            elif j0 < j1:
                for iSegment in np.arange(0, j1 - j0 + 1):
                    if iSegment == 0:
                        segment_temp = Segment()
                        segment_temp.start_point = Point(i * dx, y0)
                        segment_temp.end_point   = Point(((j0 + iSegment + 1) * dy - b) / a, (j0 + iSegment + 1) * dy)
                        segment_temp.grid_point  = Point(i, j0 + iSegment)
                        segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                        segment_temp.cal_length()
                        self.segment_list.append(segment_temp)
                    elif iSegment == j1 - j0:
                        segment_temp = Segment()
                        segment_temp.start_point = Point(i * dx, y0)
                        segment_temp.end_point   = Point(((j0 + iSegment + 1) * dy - b) / a, (j0 + iSegment + 1) * dy)
                        segment_temp.grid_point  = Point(i, j0 + iSegment)
                        segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                        segment_temp.cal_length()
                        self.segment_list.append(segment_temp)
                    else:
                        segment_temp = Segment()
                        segment_temp.start_point = Point(i * dx, y0)
                        segment_temp.end_point   = Point(((j0 + iSegment + 1) * dy - b) / a, (j0 + iSegment + 1) * dy)
                        segment_temp.grid_point  = Point(i, j0 + iSegment)
                        segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                        segment_temp.cal_length()
                        self.segment_list.append(segment_temp)                        
            elif j0 > j1:
                for iSegment in np.arange(0, j1 - j0 + 1):
                    if iSegment == 0:
                        segment_temp = Segment()
                        segment_temp.start_point = Point(i * dx, y0)
                        segment_temp.end_point   = Point(((j0 + iSegment + 1) * dy - b) / a, (j0 + iSegment + 1) * dy)
                        segment_temp.grid_point  = Point(i, j0 + iSegment)
                        segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                        segment_temp.cal_length()
                        self.segment_list.append(segment_temp)
                    elif iSegment == j1 - j0:
                        segment_temp = Segment()
                        segment_temp.start_point = Point(i * dx, y0)
                        segment_temp.end_point   = Point(((j0 + iSegment + 1) * dy - b) / a, (j0 + iSegment + 1) * dy)
                        segment_temp.grid_point  = Point(i, j0 + iSegment)
                        segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                        segment_temp.cal_length()
                        self.segment_list.append(segment_temp)
                    else:
                        segment_temp = Segment()
                        segment_temp.start_point = Point(i * dx, y0)
                        segment_temp.end_point   = Point(((j0 + iSegment + 1) * dy - b) / a, (j0 + iSegment + 1) * dy)
                        segment_temp.grid_point  = Point(i, j0 + iSegment)
                        segment_temp.normal      = Vector(0.0, 0.0, 0.0)
                        segment_temp.cal_length()
                        self.segment_list.append(segment_temp) 


# determine if two straight lines cross
# ref: https://www.cnblogs.com/wuwangchuxin0924/p/6218494.html
def  is_cross(line0, line1):
    a = line0.start_point
    b = line0.end_point
    c = line1.start_point
    d = line1.end_point
    if not(min(a.x,b.x)<=max(c.x,d.x) and min(c.y,d.y)<=max(a.y,b.y) and min(c.x,d.x)<=max(a.x,b.x) and min(a.y,b.y)<=max(c.y,d.y)):
        return False
    u = (c.x-a.x)*(b.y-a.y)-(b.x-a.x)*(c.y-a.y)
    v = (d.x-a.x)*(b.y-a.y)-(b.x-a.x)*(d.y-a.y)
    w = (a.x-c.x)*(d.y-c.y)-(d.x-c.x)*(a.y-c.y)
    z = (b.x-c.x)*(d.y-c.y)-(d.x-c.x)*(b.y-c.y)
    return (u*v<=0.0 and w*z<=0.0)


class Polygon:
    def __init__(self):
        self.line_list = []
        self.radial_line_end_point = Point(0.0, 1.0)

    def add_straight_line(self, straight_line):
        self.line_list.append(straight_line)

    def add_arc_line(self, arc_line):
        pass

    # determine is a point is in the polygon by radial line method
    # ref: https://www.cnblogs.com/On1Key/p/5673886.html
    def is_in_polygon(self, point):
        self.radial_line = Straight_line(point, self.radial_line_end_point)
        self.cross_number = 0
        for line_temp in self.line_list:
            if is_cross(self.radial_line, line_temp):
                self.cross_number += 1
        if self.cross_number%2 == 0:
            return False
        else:
            return True







def save_grid(segment_list):
    pass



if __name__ == "__main__":
    dx = 1.0e-5
    dy = 1.0e-5
    nx = 500
    ny = 500
    lx = dx * nx
    ly = dy * ny

    ny_source = 20
    ny_base = 3

    ny_wall_boundary = 150
    ny_wall_max = 200

    wall_potential = 60.0

    polygon_list = []

    point0 = Point(0.0,     ny_base*dy)
    point1 = Point(200*dx,  ny_base*dy)
    point2 = Point(0.0,     (ny_base+200)*dy)
    point3 = Point(200*dx,  (ny_base+200)*dy)

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


    point4 = Point(300*dx,  ny_base*dy)
    point5 = Point(500*dx,  ny_base*dy)
    point6 = Point(500*dx,  (ny_base+200)*dy)
    point7 = Point(350*dx,  (ny_base+200)*dy)
    point8 = Point(300*dx,  (ny_base+150)*dy)

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

    is_wall[:, 0:ny_base] = 1
    
    for i in np.arange(0, nx+1):
        for j in np.arange(0, ny+1):
            point_temp = Point(i*dx, j*dy)
            is_in_polygon = False
            for polygon_temp in polygon_list:
                if polygon_temp.is_in_polygon(point_temp):
                    is_in_polygon = True
            if is_in_polygon:
                is_wall[i,j] = 1
                
    # ============================= bndr_type ========================
    # lower boundary of source region
    bndr_type[:, ny-ny_source:ny] = 1
    bndr_val [:, ny-ny_source:ny] = 0.0

    # left and right boundary between source region and wall
    bndr_type[0, ny_wall_boundary+1:ny-ny_source-1] = 8
    bndr_type[nx,ny_wall_boundary+1:ny-ny_source-1] = 8

    # corner points of wall at left and right boudnary
    bndr_type[0, ny_wall_boundary] = 1
    bndr_type[nx,ny_wall_boundary] = 1
    bndr_val [0, ny_wall_boundary] = wall_potential
    bndr_val [nx,ny_wall_boundary] = wall_potential

    # wall surface
    for i in np.arange(1, nx):
        for j in np.arange(1, ny_wall_max):
            if is_wall[i,j] == 1 and (is_wall[i-1, j] == 0 or is_wall[i+1, j] == 0 or is_wall[i, j-1] == 0 or is_wall[i, j+1] == 0):
                bndr_type[i,j] = 1
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
    for bndr_line_temp in bndr_line_list:
        segments_temp = bndr_line_temp.generate_segments()
        for segment_temp in segments_temp:
            segment_list.append(segment_temp)

    save_grid(segment_list)
