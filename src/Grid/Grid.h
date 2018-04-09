#ifndef GRID_H
#define GRID_H

#include <cmath>

#include <vector>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"
#include "PicParams.h"


//! Class Grid: generic class allowing to define complex geometry and boundary, also used
//! to decide if the particle hit the wall
class Grid
{

public:



    //! Constructor for Grid: with no input argument
    Grid(){};

    Grid(PicParams &params){};

    //! Destructor for Grid
    virtual ~Grid() {
        ;
    } ;

    //! Virtual method used to allocate Grid
    virtual void allocateDims(){};
    virtual void geometry(){};
    virtual void computeNcp(){};
    virtual void writeGrid(){};
    virtual void readGrid(){};

    //! vector containing the dimensions of the Grid
    //! \todo private/friend/modify SmileiMPI* (JD)
    std::vector<int> dims_;
    std::vector<int> globalDims_;

    //! returns the dimension of the Grid
    inline std::vector<int> dims () {return dims_;}
    //! All arrays may be viewed as a 1D array
    //! Linearized diags


    //! pointer to the linearized array
    int* iswall_;

    int* iswall_global_;
    int* bndr_global_;
    double* bndrVal_global_;

    Field* normal_x;      // x component of the surface normal
    Field* normal_y;      // y component of the surface normal 
    Field* normal_z;      // z component of the surface normal 

    Field* normal_x_global;
    Field* normal_y_global;
    Field* normal_z_global;

    //! The number of the current point in the discrete Poisson Eqution left coefficient matrix
    int* numcp_global_;
    //>>>total number of numcp_global_ points
    int ncp;
    std::vector<int> dims_source;
    int nx,ny,nz;
    std::string gridType;


protected:

private:

};

#endif
