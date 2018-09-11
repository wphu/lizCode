/* ================================================================
solvling Ax = b for poisson equation, using KSP of PETSc package
ref: petsc-3.8.3/src/ksp/ksp/examples/tutorials/ex2.c
==================================================================*/
#ifdef petsc

#ifndef EF_SOLVER3D_PETSC_KSP_H
#define EF_SOLVER3D_PETSC_KSP_H

#include "Solver3D.h"
#include "Grid3D.h"
#include "Field.h"
#include "SmileiMPI_Cart3D.h"
#include <petscksp.h>

class ElectroMagn;
//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class EF_Solver3D_PETSc_KSP : public Solver3D
{

public:
    //! Creator for EF_SOLVER3D_PETSc_KSP
    EF_Solver3D_PETSc_KSP(PicParams& params, Grid* grid, SmileiMPI* smpi);
    virtual ~EF_Solver3D_PETSc_KSP();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);
    virtual void operator()( ElectroMagn* fields, SmileiMPI* smpi);
    void solve_PETSc_KSP(Field* rho, Field* phi);
    void finish_PETSc_KSP();
    void init_PETSc_KSP_test();
    void init_PETSc_KSP();
    void solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez);
    //>>>SuperLU parameters

    bool is_in_pestc_ksp_mpicomm()
    {
        if(petsc_ksp_mpi_rank < petsc_ksp_mpi_size)
        {
            return true;
        }
        else 
        {
            return false;
        }
    }

    // geometry parameters
    int nx, ny, nz;
    // ncp: number of calculated points, the dimension of matrix A is ncp x ncp
    int ncp;
    double dx, dy, dz, dxx;

    Grid3D* grid3D;

protected:
    double *rhsb;
    double *rhsx;
    int nnz;
    int petsc_ksp_mpi_rank;
    int petsc_ksp_mpi_size;
    double petsc_ksp_tolerance_rtol;

    // ======================== variables for PETSc_KSP ===================
    Vec            x,b,u;    /* approx solution, RHS, exact solution */
    Vec            x_seq;    /* global vector x */
    Mat            A;        /* linear system matrix */
    KSP            ksp;      /* linear solver context */
    PC             pc;       // PreCondition
    PetscRandom    rctx;     /* random number generator context */
    PetscReal      norm;     /* norm of solution error */
    PetscInt       I,J,Istart,Iend,its;
    PetscErrorCode ierr;
    PetscBool      flg = PETSC_FALSE;
    PetscScalar    v;
    PetscInt       *indices_global;
    VecScatter     ctx;







};//END class

#endif

#endif // petsc