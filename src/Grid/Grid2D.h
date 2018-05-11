#ifndef GRID2D_H
#define GRID2D_H

#include <cmath>

#include <vector>
#include <string>

#include "Grid.h"
#include "Segment.h"

using namespace std;

class SmileiMPI;

// class Grid2D used to defined a 2d vector
class Grid2D : public Grid
{

public:
    // Constructor for Grid2D: no input argument
    Grid2D();

    // Constructor for Grid2D: with the vector dimension as input argument
    Grid2D(
        PicParams &params,
        SmileiMPI *smpi,
        string grid_type,
        string gap_kind,
        int ny_source_temp,
        int ny_gapHeight_temp,
        int nx_gapWeight_temp,
        int ny_bevel_depth_temp,
        double potential_wall_temp);

    // Destructor for Grid2D
    ~Grid2D();

    // Method used to allocate a Grid2D
    void allocateDims();
    void geometry();
    void geometry_gap();
    void geometry_iter_gap();
    void computeNcp();
    void compute();

    int **iswall_2D;
  	int **iswall_global_2D;
    int **bndr_global_2D;
    double **bndrVal_global_2D;

    // define boundary lines, lines[iLine][iSegment]
    vector< vector<segment> > lines;
    int n_line;
    int n_segment_total;
    vector<int> n_segments;

    // The number of the current point in the discrete Poisson Eqution left coefficient matrix
    int **numcp_global_2D;

    // Tomakak divertor gap geometry Parameters
    std::string gapKind;
    int ny_source;
    int ny_gapHeight;
    int nx_gapWeight;
    double potential_wall;

    // ITER divetor gap 
    double bevel_depth;
    int ny_bevel_depth;


private:
    // todo{Comment what are these stuffs (MG for JD)}
    // double *data_2D;
    double dx, dy;
};

#endif
