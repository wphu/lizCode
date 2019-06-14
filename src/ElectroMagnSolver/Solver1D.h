#ifndef SOLVER1D_H
#define SOLVER1D_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver1D
//  --------------------------------------------------------------------------------------------------------------------
class Solver1D : public Solver
{

public:
    //! Creator for Solver
    Solver1D(PicParams &params) : Solver(params) 
    {
        nx_p = params.n_space[0]+1+2*params.oversize[0];
        nx_d = params.n_space[0]+2+2*params.oversize[0];

        dt_ov_dx = params.timestep / params.cell_length[0];

        bc_em_type_x  = params.bc_em_type_x;
        bc_em_value_x = params.bc_em_value_x;
    };
    virtual ~Solver1D() {};

    virtual void operator()(ElectroMagn* fields)=0;
    virtual void operator()(SmileiMPI* smpi, ElectroMagn* fields, Diagnostic* diag)=0;

protected:
    int nx_p;
    int nx_d;
    double dt_ov_dx;
    std::vector<std::string> bc_em_type_x;
    std::vector<double> bc_em_value_x;

};//END class

#endif

