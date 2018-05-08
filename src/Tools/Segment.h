#ifndef SEGMENT_H
#define SEGMENT_H

// define line segment for boundary lines
struct segment
{
    double start_point[2];
    double end_point[2];
    double length;
    // grid_point0: left-lower grid point of the rectangle region where the segment locates
    int grid_point0[2];
    // grid_point0: right-upper grid point of the rectangle region where the segment locates
    int grid_point1[2];
    // the normal vector of the segment
    double normal[3];
};


#endif