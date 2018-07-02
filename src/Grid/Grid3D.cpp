#include "Grid3D.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <fstream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Grid3D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Grid3D::Grid3D() : Grid()
{


}

// with the dimensions as input argument
Grid3D::Grid3D(
    PicParams &params,
    string grid_type,
    string gap_kind,
    int ny_source_temp,
    int ny_gapHeight_temp,
    int nx_gapWeight_temp,
    double potential_wall_temp):
    Grid(params)
{
    gridType = grid_type;
    potential_wall = potential_wall_temp;

    // number of nodes of the grid in the x-direction
    dims_.resize(3);
    globalDims_.resize(3);

    dims_[0] = params.n_space[0]+1+2*params.oversize[0];
    dims_[1] = params.n_space[1]+1+2*params.oversize[1];
    dims_[1] = params.n_space[2]+1+2*params.oversize[2];

    globalDims_[0]=params.n_space_global[0]+1;
    globalDims_[1]=params.n_space_global[1]+1;
    globalDims_[2]=params.n_space_global[2]+1;

    nx=globalDims_[0];
    ny=globalDims_[1];
    nz=globalDims_[2];

    iswall_3D.allocate_dims(dims_);
    iswall_global_3D.allocate_dims(globalDims_);
    bndr_global_3D.allocate_dims(globalDims_);
    bndrVal_global_3D.allocate_dims(globalDims_);

    if(gridType == "rectangle")
    {
            geometry();
    }
    else if(gridType == "gap")
    {
        geometry_gap();
    }
    computeNcp();
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Grid3D
// ---------------------------------------------------------------------------------------------------------------------
Grid3D::~Grid3D()
{

}


void Grid3D::allocateDims( )
{
}


//>>>no gap geometry, with source in x direction
void Grid3D::geometry( )
{
}


//>>>classical gap geometry, with source in x direction
void Grid3D::geometry_gap( )
{
}



void Grid3D::computeNcp()
{
    ncp=0;
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            if( bndr_global_3D[i][j]==0 || bndr_global_3D[i][j]==1
             || bndr_global_3D[i][j]==2 || bndr_global_3D[i][j]==8)
            {
                ncp++;
                numcp_global_3D[i][j]=ncp-1;
            }
        }
    }

}
