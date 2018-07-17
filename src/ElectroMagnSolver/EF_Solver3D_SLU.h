/* ================================================================
now only support SuperLU-4.3 serial version
==================================================================*/
#ifdef SuperLU_serial

#ifndef EF_SOLVER3D_SLU_H
#define EF_SOLVER3D_SLU_H

#include "Solver3D.h"
#include "slu_ddefs.h"
#include "Grid3D.h"
#include "Field.h"
#include "SmileiMPI_Cart3D.h"

class ElectroMagn;
//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver3D_SLU : public Solver3D
{

public:
    //! Creator for EF_SOLVER3D_SLU
    EF_Solver3D_SLU(PicParams& params, Grid* grid, SmileiMPI* smpi);
    virtual ~EF_Solver3D_SLU();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);
    void solve_SLU(Field* rho, Field* phi);
    void finishSLU();
    void initSLU_test();
    void initSLU();
    void solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez);
    //>>>SuperLU parameters


    //>>>geometry parameters
    int nx, ny, nz;
    double dx, dy, dz, dxx;

    Grid3D* grid3D;

protected:

    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
    double         *a;
    int            *asub, *xa;
    int            *perm_c; /* column permutation vector */
    int            *perm_r; /* row permutations from partial pivoting */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            m, n, nnz;
    double         *rhsb, *rhsx, *xact;
    double         *R, *C;
    double         *ferr, *berr;
    double         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    //>>>end





};//END class

#endif

#endif // for SuperLU_serial