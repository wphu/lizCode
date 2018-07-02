from template import *
import math
from scipy import constants as const

level_num = 50


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Segment:
    def __init__(self):
        self.start_point = None
        self.end_point = None
        self.length = None
        self.grid_point0 = None
        self.grid_point1 = None
        self.normal = None

    def cal_length(self):
        x_dif = self.end_point.x - self.start_point.x
        y_dif = self.end_point.y - self.start_point.y
        self.length = math.sqrt(math.pow(x_dif, 2) + math.pow(y_dif, 2))

class Straight_line:
    def __init__(self, start_point, end_point):
        self.start_point = start_point
        self.end_point = end_point
    def generate_segments(self, dx, dy):
        self.segment_list = []
        
        # the line is in y direction
        if self.start_point.x == self.end_point.x:
            normal_x = 1.0
            normal_y = 0.0
            normal_z = 0.0

            i = int(self.start_point.x / dx)
            j_start = int(self.start_point.y / dy)
            j_end = int(self.end_point.y / dy)
            if j_start <= j_end:
                for j in np.arange(j_start, j_end + 1):
                    segment_temp = Segment()
                    segment_temp.start_point = Point(i * dx, j * dy)
                    segment_temp.end_point   = Point(i * dx, (j + 1) * dy)
                    segment_temp.grid_point0  = Point(i, j)
                    segment_temp.grid_point1  = Point(i, j)
                    segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
                    segment_temp.cal_length()
                    self.segment_list.append(segment_temp)
            elif j_start > j_end:
                for j in np.arange(j_start, j_end - 1, -1):
                    segment_temp = Segment()
                    segment_temp.start_point = Point(i * dx, j * dy)
                    segment_temp.end_point   = Point(i * dx, (j - 1) * dy)
                    segment_temp.grid_point0  = Point(i, j)
                    segment_temp.grid_point1  = Point(i, j)
                    segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
                    segment_temp.cal_length()
                    self.segment_list.append(segment_temp)

        # calculate a and b for slope-intercept form of a line: y = ax + b
        else:
            a = (self.end_point.y - self.start_point.y) / (self.end_point.x - self.start_point.x)
            b = -self.start_point.x * (self.end_point.y - self.start_point.y) / (self.end_point.x - self.start_point.x) + self.start_point.y
            normal_x = -a
            normal_y = 1.0
            normal_z = 0.0
            normal_length = math.sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z)
            normal_x /= normal_length
            normal_y /= normal_length
            normal_z /= normal_length

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
            if abs(a) <= 1.0:
                for  i in np.arange(i_min, i_max + 1):
                    y0 = a * i * dx + b
                    y1 = a * (i+1) * dx + b
                    j0 = int(y0 / dy)
                    j1 = int(y1 / dy)

                    segment_temp = Segment()
                    segment_temp.start_point = Point(i * dx, y0)
                    segment_temp.end_point   = Point((i + 1) * dx, y1)
                    segment_temp.grid_point0 = Point(i, j0)
                    segment_temp.grid_point1 = Point(i, j1)
                    segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
                    segment_temp.cal_length()
                    self.segment_list.append(segment_temp) 
            else:
                for  j in np.arange(j_min, j_max + 1):
                    x0 = (j * dy - b) / a
                    x1 = ((j + 1) * dy - b) / a
                    i0 = int(x0 / dx)
                    i1 = int(x1 / dx)

                    segment_temp = Segment()
                    segment_temp.start_point = Point(x0, j * dy)
                    segment_temp.end_point   = Point(x1, (j + 1) * dy)
                    segment_temp.grid_point0 = Point(i0, j)
                    segment_temp.grid_point1 = Point(i1, j)
                    segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
                    segment_temp.cal_length()
                    self.segment_list.append(segment_temp)
        return self.segment_list

class Arc_line:
    def __init__(self, center_point, radius, start_angle, end_angle):
        self.center_point = center_point
        self.radius       = radius
        self.start_angle  = start_angle
        self.end_angle    = end_angle

    def generate_segments_for_circle(self, dx, dy):
        self.segment_list_for_circle = []
        self.angle0_list = []
        self.angle1_list = []

        # generate segment_list_for_circle

        # 0 ~ 0.25pi
        y0 = self.center_point.y
        y1 = self.center_point.y + self.radius * math.sin(0.25 * const.pi)
        j0 = int(y0 / dy)
        j1 = int(y1 / dy)
        
        for j in np.arange(j0, j1 + 1):
            y0_temp = j * dy
            y1_temp = (j + 1) * dy
            alpha0 = math.asin((y0_temp - self.center_point.y) / self.radius)
            alpha1 = math.asin((y1_temp - self.center_point.y) / self.radius)
            x0_temp = self.center_point.x + self.radius * math.cos(alpha0)
            x1_temp = self.center_point.x + self.radius * math.cos(alpha1)
            i0_temp = int(x0_temp / dx)
            i1_temp = int(x1_temp / dx)
            normal_x = 1.0
            normal_y = 0.0
            normal_z = 0.0

            segment_temp = Segment()
            segment_temp.start_point = Point(x0_temp, y0_temp)
            segment_temp.end_point   = Point(x1_temp, y1_temp)
            segment_temp.grid_point0  = Point(min(i0_temp, i1_temp), j)
            segment_temp.grid_point1  = Point(max(i0_temp, i1_temp), j)
            segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
            segment_temp.cal_length()
            self.segment_list_for_circle.append(segment_temp)
            self.angle0_list.append(alpha0)
            self.angle1_list.append(alpha1)

        # 0.25pi ~ 0.75pi
        x0 = self.center_point.x + self.radius * math.cos(0.25 * const.pi)
        x1 = self.center_point.x + self.radius * math.cos(0.75 * const.pi)
        i0 = int(x0 / dx)
        i1 = int(x1 / dx)
        
        for i in np.arange(i0, i1 - 1, -1):
            x0_temp = i * dx
            x1_temp = (i - 1) * dx
            alpha0 = math.acos((x0_temp - self.center_point.x) / self.radius)
            alpha1 = math.acos((x1_temp - self.center_point.x) / self.radius)
            y0_temp = self.center_point.y + self.radius * math.sin(alpha0)
            y1_temp = self.center_point.y + self.radius * math.sin(alpha1)
            j0_temp = int(y0_temp / dy)
            j1_temp = int(y1_temp / dy)
            normal_x = 1.0
            normal_y = 0.0
            normal_z = 0.0
            
            segment_temp = Segment()
            segment_temp.start_point = Point(x0_temp, y0_temp)
            segment_temp.end_point   = Point(x1_temp, y1_temp)
            segment_temp.grid_point0  = Point(i, min(j0_temp, j1_temp))
            segment_temp.grid_point1  = Point(i, max((j0_temp, j1_temp)))
            segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
            segment_temp.cal_length()
            self.segment_list_for_circle.append(segment_temp)
            self.angle0_list.append(alpha0)
            self.angle1_list.append(alpha1)

        # 0.75pi ~ 1.25pi
        y0 = self.center_point.y
        y1 = self.center_point.y + self.radius * math.sin(1.25 * const.pi)
        j0 = int(y0 / dy)
        j1 = int(y1 / dy)
        
        for j in np.arange(j0, j1 - 1, -1):
            y0_temp = j * dy
            y1_temp = (j - 1) * dy
            alpha0 = const.pi - math.asin((y0_temp - self.center_point.y) / self.radius)
            alpha1 = const.pi - math.asin((y1_temp - self.center_point.y) / self.radius)
            x0_temp = self.center_point.x + self.radius * math.cos(alpha0)
            x1_temp = self.center_point.x + self.radius * math.cos(alpha1)
            i0_temp = int(x0_temp / dx)
            i1_temp = int(x1_temp / dx)
            normal_x = 1.0
            normal_y = 0.0
            normal_z = 0.0

            segment_temp = Segment()
            segment_temp.start_point = Point(x0_temp, y0_temp)
            segment_temp.end_point   = Point(x1_temp, y1_temp)
            segment_temp.grid_point0  = Point(i0_temp, j)
            segment_temp.grid_point1  = Point(i1_temp, j)
            segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
            segment_temp.cal_length()
            self.segment_list_for_circle.append(segment_temp)
            self.angle0_list.append(alpha0)
            self.angle1_list.append(alpha1)

        # 1.25pi ~ 1.75pi
        x0 = self.center_point.x + self.radius * math.cos(1.25 * const.pi)
        x1 = self.center_point.x + self.radius * math.cos(1.75 * const.pi)
        i0 = int(x0 / dx)
        i1 = int(x1 / dx)
        
        for i in np.arange(i0, i1 + 1):
            x0_temp = i * dx
            x1_temp = (i + 1) * dx
            alpha0 = 2.0 * const.pi - math.acos((x0_temp - self.center_point.x) / self.radius)
            alpha1 = 2.0 * const.pi - math.acos((x1_temp - self.center_point.x) / self.radius)
            y0_temp = self.center_point.y + self.radius * math.sin(alpha0)
            y1_temp = self.center_point.y + self.radius * math.sin(alpha1)
            j0_temp = int(y0_temp / dy)
            j1_temp = int(y1_temp / dy)
            normal_x = 1.0
            normal_y = 0.0
            normal_z = 0.0
            
            segment_temp = Segment()
            segment_temp.start_point = Point(x0_temp, y0_temp)
            segment_temp.end_point   = Point(x1_temp, y1_temp)
            segment_temp.grid_point0  = Point(i, min(j0_temp, j1_temp))
            segment_temp.grid_point1  = Point(i, max(j0_temp, j1_temp))
            segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
            segment_temp.cal_length()
            self.segment_list_for_circle.append(segment_temp)
            self.angle0_list.append(alpha0)
            self.angle1_list.append(alpha1)

        # 1.75 ~ 2.0pi
        y0 = self.center_point.y + self.radius * math.sin(1.75 * const.pi)
        y1 = self.center_point.y + self.radius * math.sin(2.0  * const.pi)
        j0 = int(y0 / dy)
        j1 = int(y1 / dy)
        
        for j in np.arange(j0, j1 + 1):
            y0_temp = j * dy
            y1_temp = (j + 1) * dy
            alpha0 = 2.0 * const.pi + math.asin((y0_temp - self.center_point.y) / self.radius)
            alpha1 = 2.0 * const.pi + math.asin((y1_temp - self.center_point.y) / self.radius)
            x0_temp = self.center_point.x + self.radius * math.cos(alpha0)
            x1_temp = self.center_point.x + self.radius * math.cos(alpha1)
            i0_temp = int(x0_temp / dx)
            i1_temp = int(x1_temp / dx)
            normal_x = 1.0
            normal_y = 0.0
            normal_z = 0.0

            segment_temp = Segment()
            segment_temp.start_point = Point(x0_temp, y0_temp)
            segment_temp.end_point   = Point(x1_temp, y1_temp)
            segment_temp.grid_point0  = Point(min(i0_temp, i1_temp), j)
            segment_temp.grid_point1  = Point(max(i0_temp, i1_temp), j)
            segment_temp.normal      = Vector(normal_x, normal_y, normal_z)
            segment_temp.cal_length()
            self.segment_list_for_circle.append(segment_temp)
            self.angle0_list.append(alpha0)
            self.angle1_list.append(alpha1)

    def generate_segments(self, dx, dy):
        self.segment_list = []
        i_begin = 0
        i_end   = 0
        for i in np.arange(0,len(self.segment_list_for_circle)):
            if self.start_angle >= self.angle0_list[i] and self.start_angle <= self.angle1_list[i]:
                i_begin = i
            if self.end_angle >= self.angle0_list[i] and self.end_angle <= self.angle1_list[i]:
                i_end   = i
        for i in np.arange(i_begin, i_end+1):
            self.segment_list.append(self.segment_list_for_circle[i])

        return self.segment_list


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


class Grid2D:
    def __init__(self, dx, dy, is_wall, bndr_type, bndr_val, n_segments, segment_list):
        self.dx = dx
        self.dy = dy
        self.is_wall = is_wall[np.newaxis,:,:]
        self.bndr_type = bndr_type[np.newaxis,:,:]
        self.bndr_val = bndr_val[np.newaxis,:,:]
        self.n_segments = n_segments
        self.segment_list = segment_list

    def save_grid(self):
        self.start_point = np.zeros((len(self.segment_list), 2), dtype = 'double')
        self.end_point   = np.zeros((len(self.segment_list), 2), dtype = 'double')
        self.length      = np.zeros(len(self.segment_list), dtype = 'double')
        self.normal      = np.zeros((len(self.segment_list), 3), dtype = 'double')
        self.grid_point0 = np.zeros((len(self.segment_list), 2), dtype = 'int')
        self.grid_point1 = np.zeros((len(self.segment_list), 2), dtype = 'int')
        for i in np.arange(0, len(self.segment_list)):
            self.start_point[i,0] = self.segment_list[i].start_point.x
            self.start_point[i,1] = self.segment_list[i].start_point.y
            self.end_point[i,0]   = self.segment_list[i].end_point.x
            self.end_point[i,1]   = self.segment_list[i].end_point.y
            self.length[i]        = self.segment_list[i].length
            self.normal[i,0]      = self.segment_list[i].normal.x
            self.normal[i,1]      = self.segment_list[i].normal.y
            self.normal[i,2]      = self.segment_list[i].normal.z
            self.grid_point0[i,0] = self.segment_list[i].grid_point0.x
            self.grid_point0[i,1] = self.segment_list[i].grid_point0.y
            self.grid_point1[i,0] = self.segment_list[i].grid_point1.x
            self.grid_point1[i,1] = self.segment_list[i].grid_point1.y
        
        f = h5.File('data/grid.h5', 'w')
        f['is_wall']     = self.is_wall
        f['bndr_type']   = self.bndr_type
        f['bndr_val']    = self.bndr_val
        f['n_segments']  = self.n_segments
        f['start_point'] = self.start_point
        f['end_point']   = self.end_point
        f['length']      = self.length
        f['normal']      = self.normal
        f['grid_point0'] = self.grid_point0
        f['grid_point1'] = self.grid_point1
        f.close()

    def save_fig(self):
        ##inite the fig of matplotlib
        fig=plt.figure(figsize=(10,8))
        fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

        ##============ is_wall ======================================================
        nx = is_wall.shape[0]
        ny = is_wall.shape[1]
        dx=1.0
        dy=1.0
        x, y=np.mgrid[slice(0.0,dx*nx,dx), slice(0.0,dy*ny,dy)]

        ax0=fig.add_subplot(2,2,1)
        contourf0 = ax0.contourf(x, y, is_wall, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{is\_wall}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

        ##============ bndr_lines ======================================================
        ax0=fig.add_subplot(2,2,2)
        for i in np.arange(0, len(self.segment_list)):
            segment_x = [self.segment_list[i].start_point.x, self.segment_list[i].end_point.x]
            segment_y = [self.segment_list[i].start_point.y, self.segment_list[i].end_point.y]
            ax0.plot(segment_x, segment_y)

        ax0.set_title(r"$\mathrm{bndr\_lines}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        #ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

        ##============ bndr_type ======================================================
        ax0=fig.add_subplot(2,2,3)
        contourf0 = ax0.contourf(x, y, bndr_type, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{bndr\_type}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)

        ##============ bndr_val ======================================================
        ax0=fig.add_subplot(2,2,4)
        contourf0 = ax0.contourf(x, y, bndr_val, level_num, cmap=cm.get_cmap('jet'))

        ax0.set_title(r"$\mathrm{bndr\_val}$", color='#1f77b4', fontsize = label_fontsize)
        #ax0.axis('equal')
        ax0.set_aspect('equal', adjustable='box')
        ax0.set_xlabel('x ')
        ax0.set_ylabel('y ')
        ax0.annotate(r"$\mathbf{(d)}$", xy=get_axis_limits(ax0), annotation_clip=False)



        pdf_file_name = "grid" + ".png"
        fig.savefig(pdf_file_name, dpi = 300)


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
