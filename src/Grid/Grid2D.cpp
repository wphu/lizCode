#include "Grid2D.h"
#include "Field2D.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <fstream>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creators for Grid2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Grid2D::Grid2D() : Grid()
{


}

// with the dimensions as input argument
Grid2D::Grid2D(
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
    gapKind = gap_kind;
    ny_source = ny_source_temp;
    ny_gapHeight = ny_gapHeight_temp;
    nx_gapWeight = nx_gapWeight_temp;
    potential_wall = potential_wall_temp;

    dx = params.cell_length[0];
    dy = params.cell_length[1];

    // number of nodes of the grid in the x-direction
    dims_.resize(2);
    globalDims_.resize(2);

    dims_[0] = params.n_space[0]+1+2*params.oversize[0];
    dims_[1] = params.n_space[1]+1+2*params.oversize[1];

    globalDims_[0]=params.n_space_global[0]+1;
    globalDims_[1]=params.n_space_global[1]+1;

    nx=globalDims_[0];
    ny=globalDims_[1];

    normal_x = new Field2D(dims_, "normal_x");
    normal_y = new Field2D(dims_, "normal_y");
    normal_z = new Field2D(dims_, "normal_z");
    normal_x_global = new Field2D(globalDims_, "normal_x_global");
    normal_y_global = new Field2D(globalDims_, "normal_y_global");
    normal_z_global = new Field2D(globalDims_, "normal_z_global");

    allocateDims();
    if(gridType == "rectangle")
    {
            geometry();
    }
    else if(gridType == "gap")
    {
        geometry_gap();
    }
    else if(gridType == "iter_gap")
    {
        geometry_iter_gap();
    }
    computeNcp();

    n_line = lines.size();
    n_segments.resize(n_line);
    n_segment_total = 0;
    for(int iLine = 0; iLine < n_line; iLine++)
    {
        n_segments[iLine] = lines[iLine].size();
        n_segment_total += n_segments[iLine];
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Grid2D
// ---------------------------------------------------------------------------------------------------------------------
Grid2D::~Grid2D()
{

}


void Grid2D::allocateDims( )
{
    iswall_             = new int[dims_[0]*dims_[1]];
    iswall_global_      = new int[globalDims_[0]*globalDims_[1]];
    bndr_global_        = new int[globalDims_[0]*globalDims_[1]];
    bndrVal_global_     = new double[globalDims_[0]*globalDims_[1]];
    numcp_global_       = new int[globalDims_[0]*globalDims_[1]];

    iswall_2D           =new int*[dims_[0]];
    for (int i=0; i<dims_[0]; i++)  {
        iswall_2D[i] = iswall_ + i*dims_[1];
        for (int j=0; j<dims_[1]; j++) iswall_2D[i][j] = 0;
    }

    iswall_global_2D    =new int*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        iswall_global_2D[i] = iswall_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) iswall_global_2D[i][j] = 0;
    }

    bndr_global_2D      =new int*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        bndr_global_2D[i] = bndr_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) bndr_global_2D[i][j] = 0;
    }

    bndrVal_global_2D   =new double*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        bndrVal_global_2D[i] = bndrVal_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) bndrVal_global_2D[i][j] = 0.0;
    }

    numcp_global_2D    =new int*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        numcp_global_2D[i] = numcp_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) numcp_global_2D[i][j] = 0;
    }

}


// rectangle region
void Grid2D::geometry( )
{
    dims_source.resize(2);
    dims_source[0]=0;
    dims_source[1]=0;

    // iswall_global_2D is for particle moving, if four points of one grid is 1, then the grid is wall
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        iswall_global_2D[i][j]=0;
      }

    for(int i=0; i<nx; i++){
      iswall_global_2D[i][0]=1;
      iswall_global_2D[i][ny-1]=1;
    }

    for(int j=0; j<ny; j++){
      iswall_global_2D[0][j]=1;
      iswall_global_2D[nx-1][j]=1;
    }

    //>>>struct boundary condition
    // bndr* is for electric potential solving
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=0;
      }

    // 5 is the particle source region, usually the region has no need to solve the electric potential
    // The electric potential in the region can be set to some constant
    for(int i=0; i<dims_source[0]; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=5;
      }

    for(int i=0; i<nx; i++)
      for(int j=0; j<dims_source[1]; j++){
        bndr_global_2D[i][j]=5;
      }

    /*
    // 8 is the periodic boundary condition
    for(int i=dims_source[0]; i<nx; i++){
      bndr_global_2D[i][0]=8;
      bndr_global_2D[i][ny-1]=8;
    }
    */

    for(int i=dims_source[0]; i<nx; i++){
      bndr_global_2D[i][0] = 1;
      bndrVal_global_2D[i][0] = 0.0;

      bndr_global_2D[i][ny-1] = 1;
      bndrVal_global_2D[i][ny-1] = 0.0;
    }


    // 1 is the Dirchlet boundary condition
    for(int j=0; j<ny; j++){
      bndr_global_2D[dims_source[0]][j]=1;
      bndrVal_global_2D[dims_source[0]][j]=0.0;

      bndr_global_2D[nx-1][j]=1;
      bndrVal_global_2D[nx-1][j]=0.0;
    }


}


// classical gap geometry, with source in y direction, the source region is a rectangle below the upper boundary
void Grid2D::geometry_gap( )
{
    int bottomWall_thickness = 3;
    int nx_left_tile = 0.5*nx - 0.5*nx_gapWeight;
    int nx_right_tile = nx - nx_gapWeight;
    int ny_gap_bottom = bottomWall_thickness;
    int ny_gap_top = bottomWall_thickness + ny_gapHeight;
    // electric potential at the wall surface
    double val1 = potential_wall;
    // electric potential at source region
    double val1_source = 0.0;


    // iswall_global_2D is for particle moving, if four points of one cell are all 1, then the cell is wall
    // firstly, set the whole region as Calculation Region:  iswall_global_2D[i][j]=0
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            iswall_global_2D[i][j]=0;
        }
    }

    // set upper boundary as Wall
    for(int i=0; i<nx; i++)
    {
        iswall_global_2D[i][ny-1]=1;
    }

    // set lower boundary as Wall, the lower boundary has thickness, two or three cell should be OK, defined by bottomWall_thickness
    for(int i=0; i<nx; i++)
    {
        for(int j = 0; j < bottomWall_thickness + 1; j++)
        {
            iswall_global_2D[i][j] = 1;
        }
        
    }

    // set left and right boundary as Wall
    for(int j=0; j<ny; j++)
    {
        iswall_global_2D[0][j]=1;
        iswall_global_2D[nx-1][j]=1;
    }

    // set the left tile as Wall
    for(int i = 0; i < 0.5*nx - 0.5*nx_gapWeight; i++)
    {
        for(int j = 0; j < bottomWall_thickness + ny_gapHeight; j++)
        {
            iswall_global_2D[i][j] = 1;
        }
    }

    // set the right tile as Wall
    for(int i =  0.5*nx + 0.5*nx_gapWeight; i < nx; i++)
    {
        for(int j = 0; j < bottomWall_thickness + ny_gapHeight; j++)
        {
            iswall_global_2D[i][j] = 1;
        }
    }


    //============== construct boundary condition ===================================
    // bndr* is for electric potential solving 
    // 5 is the particle source region or wall region, usually the region has no need to solve the electric potential
    //      The electric potential in the region can be set to some constant   
    // 8 is the periodic boundary condition
    // 1 is the Dirichlet boundary condition
    // 0 is the Calculation Region  
    
    // firstly, set the whole region as Calculation Resion
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            bndr_global_2D[i][j] = 0;
        }
    }

    // set the source region, a rectangle below the upper boundary
    for(int i=0; i<nx; i++)
    {
        for(int j = ny - ny_source; j < ny; j++)
        {
            bndr_global_2D[i][j] = 5;
            bndrVal_global_2D[i][j] = val1_source;
        }
    }

    // set the left and right boudary
    for(int j = 0; j < ny - ny_source; j++){
      bndr_global_2D[0][j]=8;
      bndr_global_2D[nx-1][j]=8;
    }

    // set the source region surface boudary
    for(int i = 0; i < nx; i++){
      bndr_global_2D    [i][ny - ny_source -1] = 1;
      bndrVal_global_2D [i][ny - ny_source -1] = val1_source;
    }

    // set the wall surface boundary
    for(int j = 0; j < ny_gapHeight; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            // set the boundary corner points of two tiles
            if( (i == 0 || i == nx - 1) && j == ny_gapHeight -1 ) 
            {
                bndr_global_2D[i][j] = 1;
                bndrVal_global_2D[i][j] = val1;
            }
            // set the left and right boundary of tile surface
            else if( (i == 0 || i == nx - 1) && j < ny_gapHeight -1 ) 
            {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
            // set the lower boundary, which is never used for field solving, convenient for particle moving
            else if( j == 0 ) 
            {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
            // set tile surface as the Direchlet boundary condition
            else if( iswall_global_2D[i][j] == 1 && ( iswall_global_2D[i-1][j] == 0 || iswall_global_2D[i][j-1] == 0
             || iswall_global_2D[i+1][j] == 0 || iswall_global_2D[i][j+1] == 0 ) ) 
            {
                bndr_global_2D[i][j] = 1;
                bndrVal_global_2D[i][j] = val1;
            }
            // set the remaining part as Wall Interior Domain
            else if( iswall_global_2D[i][j] == 1 ) 
            {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
        }
    }


    //============== construct the unit vector of tile surface normal ===================================
    // set the top surface of the left tile
    for(int i = 0; i < nx_left_tile; i++)
    {
        (*normal_x_global)(i, ny_gap_top) = 0.0;
        (*normal_y_global)(i, ny_gap_top) = 1.0;
        (*normal_z_global)(i, ny_gap_top) = 0.0;
    }
    // set the right surface of the left tile
    for(int j = ny_gap_bottom + 1; j <= ny_gap_top; j++)
    {
        (*normal_x_global)(nx_left_tile, j) = 1.0;
        (*normal_y_global)(nx_left_tile, j) = 0.0;
        (*normal_z_global)(nx_left_tile, j) = 0.0;
    }
    // set the bottom surface of the gap
    for(int i = nx_left_tile; i < nx_right_tile; i++)
    {
        (*normal_x_global)(i, ny_gap_bottom) = 0.0;
        (*normal_y_global)(i, ny_gap_bottom) = 1.0;
        (*normal_z_global)(i, ny_gap_bottom) = 0.0;
    }
    // set the left surface of the right tile
    for(int j = ny_gap_bottom; j < ny_gap_top; j++)
    {
        (*normal_x_global)(nx_right_tile, j) = -1.0;
        (*normal_y_global)(nx_right_tile, j) = 0.0;
        (*normal_z_global)(nx_right_tile, j) = 0.0;
    }
    // set the top surface of the right tile
    for(int i = nx_right_tile; i < nx; i++)
    {
        (*normal_x_global)(i, ny_gap_top) = 0.0;
        (*normal_y_global)(i, ny_gap_top) = 1.0;
        (*normal_z_global)(i, ny_gap_top) = 0.0;
    }

}

// iter divetor gap geometry, with bevel top surface in toroidal direction
void Grid2D::geometry_iter_gap( )
{
    bevel_depth = 20.0 * dx;

    int bottomWall_thickness = 3;
    int nx_left_tile = 0.5*nx - 0.5*nx_gapWeight;
    int nx_right_tile = nx - nx_gapWeight;
    int ny_gap_bottom = bottomWall_thickness;
    int ny_gap_top0 = bottomWall_thickness + ny_gapHeight;
    int ny_gap_top1 = bottomWall_thickness + ny_gapHeight + 0.5 * bevel_depth / dy;
    int ny_gap_top2 = bottomWall_thickness + ny_gapHeight + bevel_depth / dy;
    // electric potential at the wall surface
    double val1 = potential_wall;
    // electric potential at source region
    double val1_source = 0.0;

    double x1, x2, y1, y2, a, b;
    double normal_x, normal_y, normal_z, normal_length;
    double y_lower, y_upper;
    int j_lower, j_upper;
    int nSegment;

    // iswall_global_2D is for particle moving, if four points of one cell are all 1, then the cell is wall
    // firstly, set the whole region as Calculation Region:  iswall_global_2D[i][j]=0
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            iswall_global_2D[i][j]=0;
        }
    }

    // set upper boundary as Wall
    for(int i=0; i<nx; i++)
    {
        iswall_global_2D[i][ny-1]=1;
    }

    // set lower boundary as Wall, the lower boundary has thickness, two or three cell should be OK, defined by bottomWall_thickness
    for(int i=0; i<nx; i++)
    {
        for(int j = 0; j < bottomWall_thickness + 1; j++)
        {
            iswall_global_2D[i][j] = 1;
        }
        
    }

    // set left and right boundary as Wall
    for(int j=0; j<ny; j++)
    {
        iswall_global_2D[0][j]=1;
        iswall_global_2D[nx-1][j]=1;
    }

    // set the left tile as Wall
    x1 = 0.0;
    y1 = ny_gap_top1 * dy;
    x2 = nx_left_tile * dx;
    y2 = ny_gap_top2 * dy;
    a = (y2-y1)/(x2-x1);
    b = -x1*(y2-y1)/(x2-x1) + y1;
    for(int i = 0; i <= nx_left_tile; i++)
    {
        for(int j = 0; j <= ny_gap_top2; j++)
        {
            if(j * dy <= a * i * dx + b)
            {
                iswall_global_2D[i][j] = 1;
            }
        }
    }

    // set the right tile as Wall
    x1 = nx_right_tile * dx;
    y1 = ny_gap_top0 * dy;
    x2 = (nx - 1) * dx;
    y2 = ny_gap_top1 * dy;
    a = (y2-y1)/(x2-x1);
    b = -x1*(y2-y1)/(x2-x1) + y1;
    for(int i =  nx_right_tile; i < nx; i++)
    {
        for(int j = 0; j <= ny_gap_top1; j++)
        {
            if(j * dy <= a * i * dx + b)
            {
                iswall_global_2D[i][j] = 1;
            }
        }
    }


    //============== construct boundary condition ===================================
    // bndr* is for electric potential solving 
    // 5 is the particle source region or wall region, usually the region has no need to solve the electric potential
    //      The electric potential in the region can be set to some constant   
    // 8 is the periodic boundary condition
    // 1 is the Dirichlet boundary condition
    // 0 is the Calculation Region  
    
    // firstly, set the whole region as Calculation Region
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            bndr_global_2D[i][j] = 0;
        }
    }

    // set the source region, a rectangle below the upper boundary
    for(int i=0; i<nx; i++)
    {
        for(int j = ny - ny_source; j < ny; j++)
        {
            bndr_global_2D[i][j] = 5;
            bndrVal_global_2D[i][j] = val1_source;
        }
    }

    // set the left and right boudary
    for(int j = 0; j < ny - ny_source; j++){
      bndr_global_2D[0][j]=8;
      bndr_global_2D[nx-1][j]=8;
    }

    // set the source region surface boudary
    for(int i = 0; i < nx; i++){
      bndr_global_2D    [i][ny - ny_source -1] = 1;
      bndrVal_global_2D [i][ny - ny_source -1] = val1_source;
    }

    // set the wall surface boundary
    for(int j = 0; j <= ny_gap_top2; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            if( i == 0 || i == nx - 1) 
            {
                // set the boundary corner points of two tiles
                if(j == ny_gap_top1)
                {
                    bndr_global_2D[i][j] = 1;
                    bndrVal_global_2D[i][j] = val1;
                }
                // set the left and right boundary of tile surface
                else if(j < ny_gap_top1)
                {
                    bndr_global_2D[i][j] = 5;
                    bndrVal_global_2D[i][j] = val1;
                }
            }
            // set the lower boundary, which is never used for field solving, convenient for particle moving
            else if( j == 0 ) 
            {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
            // set tile surface as the Direchlet boundary condition
            else if( iswall_global_2D[i][j] == 1 && ( iswall_global_2D[i-1][j] == 0 || iswall_global_2D[i][j-1] == 0
             || iswall_global_2D[i+1][j] == 0 || iswall_global_2D[i][j+1] == 0 ) ) 
            {
                bndr_global_2D[i][j] = 1;
                bndrVal_global_2D[i][j] = val1;
            }
            // set the remaining part as Wall Interior Domain
            else if( iswall_global_2D[i][j] == 1 ) 
            {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
        }
    }


    //============== construct boundary lines ===================================
    lines.resize(5);

    // set the top surface of the left tile
    x1 = 0.0;
    y1 = ny_gap_top1 * dy;
    x2 = nx_left_tile * dx;
    y2 = ny_gap_top2 * dy;
    a = (y2-y1)/(x2-x1);
    b = -x1*(y2-y1)/(x2-x1) + y1;
    normal_x = -bevel_depth * bevel_depth / nx_left_tile;
    normal_y = bevel_depth;
    normal_z = 0.0;
    normal_length = sqrt( normal_x * normal_x + normal_y * normal_y );
    normal_x /= normal_length;
    normal_y /= normal_length;
    for(int i = 0; i < nx_left_tile; i++)
    {
        y_lower = a * i * dx + b;
        y_upper = a * (i+1) * dx + b;
        j_lower = y_lower / dy;
        j_upper = y_upper / dy;
        if(j_upper == j_lower)
        {
            segment seg;
            seg.start_point[0] = i * dx;
            seg.start_point[1] = y_lower;
            seg.end_point[0] = (i+1) * dx;
            seg.end_point[1] = y_upper;
            seg.grid_point[0] = i;
            seg.grid_point[1] = j_lower;
            seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
            seg.normal[0] = normal_x;
            seg.normal[1] = normal_y;
            seg.normal[2] = normal_z;
            if(seg.length == 0.0)
            {
                continue;
            }
            lines[0].push_back(seg);
        }
        else
        {
            for(int iSegment = 0; iSegment < j_upper - j_lower + 1; iSegment++)
            {
                if(iSegment == 0)
                {
                    segment seg;
                    seg.start_point[0] = i * dx;
                    seg.start_point[1] = y_lower;
                    seg.end_point[0] = ((j_lower + iSegment + 1) * dy - b) / a; 
                    seg.end_point[1] = (j_lower + iSegment + 1) * dy;
                    seg.grid_point[0] = i;
                    seg.grid_point[1] = j_lower + iSegment;
                    seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
                    seg.normal[0] = normal_x;
                    seg.normal[1] = normal_y;
                    seg.normal[2] = normal_z;
                    if(seg.length == 0.0)
                    {
                        continue;
                    }
                    lines[0].push_back(seg);
                }
                else if(iSegment == j_upper - j_lower)
                {
                    segment seg;
                    seg.start_point[0] = ((j_lower + iSegment) * dy - b) / a;
                    seg.start_point[1] = (j_lower + iSegment) * dy;
                    seg.end_point[0] = (i + 1) * dx; 
                    seg.end_point[1] = a * (i + 1) * dx + b;
                    seg.grid_point[0] = i;
                    seg.grid_point[1] = j_lower + iSegment;  
                    seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
                    seg.normal[0] = normal_x;
                    seg.normal[1] = normal_y;
                    seg.normal[2] = normal_z;
                    if(seg.length == 0.0)
                    {
                        continue;
                    }
                    lines[0].push_back(seg);
                }
                else
                {
                    segment seg;
                    seg.start_point[0] = ((j_lower + iSegment) * dy - b) / a;
                    seg.start_point[1] = (j_lower + iSegment) * dy;
                    seg.end_point[0] = ((j_lower + iSegment + 1) * dy - b) / a;
                    seg.end_point[1] = (j_lower + iSegment + 1) * dy;
                    seg.grid_point[0] = i;
                    seg.grid_point[1] = j_lower + iSegment;  
                    seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
                    seg.normal[0] = normal_x;
                    seg.normal[1] = normal_y;
                    seg.normal[2] = normal_z;
                    if(seg.length == 0.0)
                    {
                        continue;
                    }
                    lines[0].push_back(seg);
                }
            }
        }
    }

    // set the right surface of the left tile
    normal_x = 1.0;
    normal_y = 0.0;
    normal_z = 0.0;
    if(a * nx_left_tile * dx + b == ny_gap_top2 * dy)
    {
        nSegment = ny_gap_top2 - bottomWall_thickness;
        segment seg;
        seg.start_point[0] = nx_left_tile * dx;
        seg.start_point[1] = ny_gap_top2 * dy;
        seg.end_point[0] = nx_left_tile * dx;
        seg.end_point[1] = (ny_gap_top2 - 1) * dy;
        seg.grid_point[0] = nx_left_tile;
        seg.grid_point[1] = ny_gap_top2 - 1;
        seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
        seg.normal[0] = normal_x;
        seg.normal[1] = normal_y;
        seg.normal[2] = normal_z;
        if(seg.length == 0.0)
        {
            continue;
        }
        lines[1].push_back(seg);
        for(int iSegment = 1; iSegment < nSegment; iSegment++)
        {
            segment seg;
            seg.start_point[0] = nx_left_tile * dx;
            seg.start_point[1] = (ny_gap_top2 - iSegment) * dy;
            seg.end_point[0] = nx_left_tile * dx;
            seg.end_point[1] = (ny_gap_top2 - iSegment - 1) * dy;
            seg.grid_point[0] = nx_left_tile;
            seg.grid_point[1] = ny_gap_top2 - iSegment - 1;
            seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
            seg.normal[0] = normal_x;
            seg.normal[1] = normal_y;
            seg.normal[2] = normal_z;
            if(seg.length == 0.0)
            {
                continue;
            }
            lines[1].push_back(seg);
        }
    }
    else
    {
        nSegment = ny_gap_top2 - bottomWall_thickness + 1;
        segment seg;
        seg.start_point[0] = nx_left_tile * dx;
        seg.start_point[1] = a * nx_left_tile * dx + b;
        seg.end_point[0] = nx_left_tile * dx;
        seg.end_point[1] = (ny_gap_top2) * dy;
        seg.grid_point[0] = nx_left_tile;
        seg.grid_point[1] = ny_gap_top2;
        seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
        seg.normal[0] = normal_x;
        seg.normal[1] = normal_y;
        seg.normal[2] = normal_z;
        if(seg.length == 0.0)
        {
            continue;
        }
        lines[1].push_back(seg);
        for(int iSegment = 1; iSegment < nSegment; iSegment++)
        {
            segment seg;
            seg.start_point[0] = nx_left_tile * dx;
            seg.start_point[1] = (ny_gap_top2 - iSegment + 1) * dy;
            seg.end_point[0] = nx_left_tile * dx;
            seg.end_point[1] = (ny_gap_top2 - iSegment) * dy;
            seg.grid_point[0] = nx_left_tile;
            seg.grid_point[1] = ny_gap_top2 - iSegment;
            seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
            seg.normal[0] = normal_x;
            seg.normal[1] = normal_y;
            seg.normal[2] = normal_z;
            if(seg.length == 0.0)
            {
                continue;
            }
            lines[1].push_back(seg);
        }
    }

    // set the bottom surface of the gap
    normal_x = 0.0;
    normal_y = 1.0;
    normal_z = 0.0;
    for(int i = nx_left_tile; i < nx_right_tile; i++)
    {
        segment seg;
        seg.start_point[0] = i * dx;
        seg.start_point[1] = bottomWall_thickness * dy;
        seg.end_point[0] = (i + 1) * dx;
        seg.end_point[1] = bottomWall_thickness * dy;
        seg.grid_point[0] = i;
        seg.grid_point[1] = bottomWall_thickness - 1;
        seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
        seg.normal[0] = normal_x;
        seg.normal[1] = normal_y;
        seg.normal[2] = normal_z;
        if(seg.length == 0.0)
        {
            continue;
        }
        lines[2].push_back(seg);
    }

    // set the left surface of the right tile
    x1 = nx_right_tile * dx;
    y1 = ny_gap_top0 * dy;
    x2 = nx * dx;
    y2 = ny_gap_top1 * dy;
    a = (y2-y1)/(x2-x1);
    b = -x1*(y2-y1)/(x2-x1) + y1;
    normal_x = -1.0;
    normal_y = 0.0;
    normal_z = 0.0;
    for(int j = bottomWall_thickness; j < ny_gap_top0; j++)
    {
        segment seg;
        seg.start_point[0] = nx_right_tile * dx;
        seg.start_point[1] = j * dy;
        seg.end_point[0] = nx_right_tile * dx;
        seg.end_point[1] = (j + 1) * dy;
        seg.grid_point[0] = nx_right_tile - 1;
        seg.grid_point[1] = j;
        seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
        seg.normal[0] = normal_x;
        seg.normal[1] = normal_y;
        seg.normal[2] = normal_z;
        if(seg.length == 0.0)
        {
            continue;
        }
        lines[3].push_back(seg);
    }
    if(a * nx_right_tile * dx + b != ny_gap_top0 * dy)
    {
        segment seg;
        seg.start_point[0] = nx_right_tile * dx;
        seg.start_point[1] = ny_gap_top0 * dy;
        seg.end_point[0] = nx_right_tile * dx;
        seg.end_point[1] = a * nx_right_tile * dx + b;
        seg.grid_point[0] = nx_right_tile - 1;
        seg.grid_point[1] = ny_gap_top0;
        seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
        seg.normal[0] = normal_x;
        seg.normal[1] = normal_y;
        seg.normal[2] = normal_z;
        if(seg.length == 0.0)
        {
            continue;
        }
        lines[3].push_back(seg);
    }

    // set the top surface of the right tile
    normal_x = -bevel_depth * bevel_depth / nx_left_tile;
    normal_y = bevel_depth;
    normal_z = 0.0;
    normal_length = sqrt( normal_x * normal_x + normal_y * normal_y );
    normal_x /= normal_length;
    normal_y /= normal_length;
    for(int i = nx_right_tile; i < nx; i++)
    {
        y_lower = a * i * dx + b;
        y_upper = a * (i+1) * dx + b;
        j_lower = y_lower / dy;
        j_upper = y_upper / dy;
        if(j_upper == j_lower)
        {
            segment seg;
            seg.start_point[0] = i * dx;
            seg.start_point[1] = y_lower;
            seg.end_point[0] = (i+1) * dx;
            seg.end_point[1] = y_upper;
            seg.grid_point[0] = i;
            seg.grid_point[1] = j_lower;
            seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
            seg.normal[0] = normal_x;
            seg.normal[1] = normal_y;
            seg.normal[2] = normal_z;
            if(seg.length == 0.0)
            {
                continue;
            }
            lines[4].push_back(seg);
        }
        else
        {
            for(int iSegment = 0; iSegment < j_upper - j_lower + 1; iSegment++)
            {
                if(iSegment == 0)
                {
                    segment seg;
                    seg.start_point[0] = i * dx;
                    seg.start_point[1] = y_lower;
                    seg.end_point[0] = ((j_lower + iSegment + 1) * dy - b) / a; 
                    seg.end_point[1] = (j_lower + iSegment + 1) * dy;
                    seg.grid_point[0] = i;
                    seg.grid_point[1] = j_lower + iSegment;
                    seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
                    seg.normal[0] = normal_x;
                    seg.normal[1] = normal_y;
                    seg.normal[2] = normal_z;
                    if(seg.length == 0.0)
                    {
                        continue;
                    }
                    lines[4].push_back(seg);
                }
                else if(iSegment == j_upper - j_lower)
                {
                    segment seg;
                    seg.start_point[0] = ((j_lower + iSegment) * dy - b) / a;
                    seg.start_point[1] = (j_lower + iSegment) * dy;
                    seg.end_point[0] = (i + 1) * dx; 
                    seg.end_point[1] = a * (i + 1) * dx + b;
                    seg.grid_point[0] = i;
                    seg.grid_point[1] = j_lower + iSegment;  
                    seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
                    seg.normal[0] = normal_x;
                    seg.normal[1] = normal_y;
                    seg.normal[2] = normal_z;
                    if(seg.length == 0.0)
                    {
                        continue;
                    }
                    lines[4].push_back(seg);
                }
                else
                {
                    segment seg;
                    seg.start_point[0] = ((j_lower + iSegment) * dy - b) / a;
                    seg.start_point[1] = (j_lower + iSegment) * dy;
                    seg.end_point[0] = ((j_lower + iSegment + 1) * dy - b) / a;
                    seg.end_point[1] = (j_lower + iSegment + 1) * dy;
                    seg.grid_point[0] = i;
                    seg.grid_point[1] = j_lower + iSegment;  
                    seg.length = sqrt( pow((seg.end_point[0] - seg.start_point[0]), 2) + pow((seg.end_point[1] - seg.start_point[1]), 2) );
                    seg.normal[0] = normal_x;
                    seg.normal[1] = normal_y;
                    seg.normal[2] = normal_z;
                    if(seg.length == 0.0)
                    {
                        continue;
                    }
                    lines[4].push_back(seg);
                }
            }
        }
    }
}


void Grid2D::computeNcp()
{

    ncp=0;
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        if( bndr_global_2D[i][j]==0 || bndr_global_2D[i][j]==1
        || bndr_global_2D[i][j]==2 || bndr_global_2D[i][j]==8){
          ncp++;
          numcp_global_2D[i][j]=ncp-1;

        }
      }


}
