/* ================================================================
MPI parallel SuperLU version: now only SuperLU_DIST-5.3.0 is used
==================================================================*/
#ifdef SuperLU_mpi

#ifndef EF_SOLVER3D_SLU_DIST_H
#define EF_SOLVER3D_SLU_DIST_H

#include "Solver3D.h"
#include "superlu_ddefs.h"
#include "Grid3D.h"
#include "Field.h"
#include "SmileiMPI_Cart3D.h"

class ElectroMagn;
//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver3D_SLU_DIST : public Solver3D
{

public:
    //! Creator for EF_SOLVER3D_SLU
    EF_Solver3D_SLU_DIST(PicParams& params, Grid* grid, SmileiMPI* smpi);
    virtual ~EF_Solver3D_SLU_DIST();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);
    void solve_SLU(Field* rho, Field* phi);
    void finishSLU();
    void initSLU();
    void solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez);
    //>>>SuperLU parameters


    //>>>geometry parameters
    int nx, ny, nz;
    double dx, dy, dz, dxx;

    Grid3D* grid3D;

protected:

    superlu_dist_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    gridinfo_t grid;
    double   *berr;
    double   *a, *b1, *xtrue;
    int_t    *asub, *xa;
    int_t    i, j, m, n, nnz;
    int_t    nprow, npcol;
    int      iam, info, ldb, ldx, nrhs;
    char     trans[1];
    char     **cpp, c;
    FILE *fp, *fopen();

    /* prototypes */
    /*
    extern void LUstructInit(const int_t, LUstruct_t *);
    extern void LUstructFree(LUstruct_t *);
    extern void Destroy_LU(int_t, gridinfo_t *, LUstruct_t *);
    */

};//END class

#endif

#endif // for SuperLU_mpi