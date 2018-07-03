#ifndef GRID3D_H
#define GRID3D_H

#include <cmath>

#include <vector>
#include <string>

#include "Grid.h"
#include "Array3D.h"

using namespace std;

class SmileiMPI;

//! class Grid3D used to defined a 3D vector
class Grid3D : public Grid
{

public:
    //! Constructor for Grid3D: no input argument
    Grid3D();

    //! Constructor for Grid3D: with the vector dimension as input argument
    Grid3D(
        PicParams &params,
        SmileiMPI *smpi,
        string grid_type,
        string gap_kind,
        int ny_source_temp,
        int ny_gapHeight_temp,
        int nx_gapWeight_temp,
        int ny_bevel_depth_temp,
        double potential_wall_temp);

    // Destructor for Grid3D
    ~Grid3D();

    // Method used to allocate a Grid3D
    void allocateDims();
    void geometry();
    void geometry_gap();
    void computeNcp();
    Array3D<int> iswall_3D;
  	Array3D<int> iswall_global_3D;
    Array3D<int> bndr_global_3D;
    Array3D<double> bndrVal_global_3D;
    // The number of the current point in the discrete Poisson Eqution left coefficient matrix
    Array3D<int> numcp_global_3D;

    double potential_wall;


private:
    // \todo{Comment what are these stuffs (MG for JD)}
    //double *data_3D;
};

#endif
