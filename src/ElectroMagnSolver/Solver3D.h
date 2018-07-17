#ifndef SOLVER3D_H
#define SOLVER3D_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver3D
//  --------------------------------------------------------------------------------------------------------------------
class Solver3D : public Solver
{

public:
    //! Creator for Solver
    Solver3D(PicParams &params) : Solver(params) 
    {
        nx_p = params.n_space[0]+1+2*params.oversize[0];
        nx_d = params.n_space[0]+2+2*params.oversize[0];
        ny_p = params.n_space[1]+1+2*params.oversize[1];
        ny_d = params.n_space[1]+2+2*params.oversize[1];
        nz_p = params.n_space[2]+1+2*params.oversize[2];
        nz_d = params.n_space[2]+2+2*params.oversize[2];

        dt_ov_dx = params.timestep / params.cell_length[0];
        dt_ov_dy = params.timestep / params.cell_length[1];
        dt_ov_dz = params.timestep / params.cell_length[2];
    };
    virtual ~Solver3D() {};

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields)=0;
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi)=0;
    virtual void solve_SLU(Field* rho, Field* phi){};
    virtual void finishSLU(){};
    virtual void initSLU_test(){};
    virtual void initSLU(){};
    virtual void solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez){};

protected:
    int nx_p;
    int nx_d;
    int ny_p;
    int ny_d;
    int nz_p;
    int nz_d;
    double dt_ov_dy;
    double dt_ov_dx;
    double dt_ov_dz;

};//END class

#endif
