#ifndef SEGMENT_H
#define SEGMENT_H

// define line segment for boundary lines
struct segment
{
    double start_point[2];
    double end_point[2];
    double length;
    int grid_point[2];
    double normal[3];
};


#endif