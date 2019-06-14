/* ================================================================
umfpack is included in the SuiteSparse packgade, 
so install SuiteSparse
==================================================================*/
#ifdef umfpack

#ifndef EF_SOLVER3D_UMFPACK_H
#define EF_SOLVER3D_UMFPACK_H

#include "Solver3D.h"
#include "umfpack.h"
#include "Grid3D.h"
#include "Field.h"
#include "SmileiMPI_Cart3D.h"

class ElectroMagn;
//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver3D_UMFPACK : public Solver3D
{

public:
    //! Creator for EF_SOLVER3D_UMFPACK
    EF_Solver3D_UMFPACK(PicParams& params, Grid* grid, SmileiMPI* smpi);
    virtual ~EF_Solver3D_UMFPACK();

    //! Overloading of () operator
    virtual void operator()(ElectroMagn* fields);
    virtual void operator()(SmileiMPI* smpi, ElectroMagn* fields, Diagnostic* diag);
    void solve_UMFPACK(Field* rho, Field* phi);
    void finishUMFPACK();
    void initUMFPACK_test();
    void initUMFPACK();
    void solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez);
    //>>>SuperLU parameters


    //>>>geometry parameters
    int nx, ny, nz;
    double dx, dy, dz, dxx;

    Grid3D* grid3D;

protected:

    // for umfpack
    SuiteSparse_long             n;
    SuiteSparse_long             nnz;
    SuiteSparse_long             *Ap;
    SuiteSparse_long             *Ai;
    double          *Ax;
    double          *rhsb;
    double          *rhsx;
    double          *null;
    void            *Symbolic;
    void            *Numeric;
    double          Info[UMFPACK_INFO];
    double          Control[UMFPACK_CONTROL];
    SuiteSparse_long             status;






};//END class

#endif

#endif // for SuperLU_serial