#ifdef SuperLU_serial

#include "EF_Solver2D_SLU.h"

#include "ElectroMagn.h"
#include "Field2D.h"
#include <memory>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


EF_Solver2D_SLU::EF_Solver2D_SLU(PicParams& params, Grid* grid, SmileiMPI* smpi):
Solver2D(params)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);


    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dxy = dx * dy;

    grid2D = static_cast<Grid2D*>(grid);

    if(smpi2D->isMaster())
    {
        initSLU();
    }


}


EF_Solver2D_SLU::~EF_Solver2D_SLU()
{
}

void EF_Solver2D_SLU::operator() ( ElectroMagn* fields )
{

}



void EF_Solver2D_SLU::operator() ( ElectroMagn* fields , SmileiMPI* smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(fields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(fields->Ey_);


    Field2D* rho2D           = static_cast<Field2D*>(fields->rho_);
    Field2D* rho2D_global    = static_cast<Field2D*>(fields->rho_global);
    Field2D* phi2D_global    = static_cast<Field2D*>(fields->phi_global);
    Field2D* Ex2D_global    = static_cast<Field2D*>(fields->Ex_global);
    Field2D* Ey2D_global    = static_cast<Field2D*>(fields->Ey_global);

    smpi2D->barrier();
    smpi2D->gatherRho(rho2D_global, rho2D);

    if(smpi2D->isMaster()){
        solve_SLU(rho2D_global, phi2D_global);
        solve_Exy(phi2D_global, Ex2D_global, Ey2D_global);
    }

    smpi2D->barrier();
    smpi2D->scatterField(Ex2D_global, Ex2D);
    smpi2D->scatterField(Ey2D_global, Ey2D);
}


void EF_Solver2D_SLU::initSLU_test()
{

    SuperMatrix A, L, U, B;
    double   *a, *rhs;
    double   s, u, p, e, r, l;
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      nrhs, info, i, m, n, nnz, permc_spec;
    superlu_options_t options;
    SuperLUStat_t stat;
    /* Initialize matrix A. */
    m = n = 5;
    nnz = 12;
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
    s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
    a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
    a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
    asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
    asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
    asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
    xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* Create right-hand side matrix B. */
    nrhs = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    for (i = 0; i < m; ++i) rhs[i] = 1.0;
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    set_default_options(&options);
    options.ColPerm = NATURAL;

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Solve the linear system. */
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    printf("LU factorization: dgssvx() returns info %d\n", info);


}


void EF_Solver2D_SLU::initSLU()
{
    vector< vector<double> > val;
    vector< vector<int> >    row;

    int i,j,k,ii,ll,kk,v,hu,hd,hr,hl,i_ncp,i_nnz,nz_col,i_val;

    val.resize(grid2D->ncp);
    row.resize(grid2D->ncp);

    nnz=0;
    ii=0;
    v=0;
    nx=grid2D->nx;
    ny=grid2D->ny;

    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            // normal points in the calculation region
            if(grid2D->bndr_global_2D[i][j]==0) 
            {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i+1][j] - grid2D->numcp_global_2D[i][j];
                /*
                for(k=0; k<grid2D->ncp; k++)
                {
                    b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii+hr]=1.0; b[ii+1]=1.0;
                */

                nnz=nnz+5;
                
                val[ii].push_back(-4.0);
                row[ii].push_back(ii);
                val[ii-hl].push_back(1.0);
                row[ii-hl].push_back(ii);
                val[ii-1].push_back(1.0);
                row[ii-1].push_back(ii);
                val[ii+hr].push_back(1.0);
                row[ii+hr].push_back(ii);
                val[ii+1].push_back(1.0);
                row[ii+1].push_back(ii);

                /*
                for(k=0; k<grid2D->ncp; k++)
                {
                    if(b[k] != 0.0)
                    {
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                */

                ii++;
            }

            // Dirchlet boudnary points
            else if(grid2D->bndr_global_2D[i][j]==1) 
            {
                nnz++;

                val[ii].push_back(1.0);
                row[ii].push_back(ii);

                ii++;
            }

            // periodic boudnary points at lowwer boudary in y direction
            else if( grid2D->bndr_global_2D[i][j]==8 && j==0) 
            {
                hu = grid2D->numcp_global_2D[i][ny-1] - grid2D->numcp_global_2D[i][j];
                nnz=nnz+2;

                val[ii].push_back(1.0);
                row[ii].push_back(ii);
                val[ii+hu].push_back(-1.0);
                row[ii+hu].push_back(ii);

                ii++;
            }

            // periodic boudnary points at upper boudary in y direction
            else if ( grid2D->bndr_global_2D[i][j] == 8 && j == ny-1 ) 
            {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i+1][j] - grid2D->numcp_global_2D[i][j];
                hd = grid2D->numcp_global_2D[i][ny-1] - grid2D->numcp_global_2D[i][1];
                nnz=nnz+5;

                val[ii].push_back(-4.0);
                row[ii].push_back(ii);
                val[ii-hl].push_back(1.0);
                row[ii-hl].push_back(ii);
                val[ii-1].push_back(1.0);
                row[ii-1].push_back(ii);
                val[ii+hr].push_back(1.0);
                row[ii+hr].push_back(ii);
                val[ii-hd].push_back(1.0);
                row[ii-hd].push_back(ii);

                ii++;
            }

            // periodic boudnary points at left boudary in x direction
            else if( grid2D->bndr_global_2D[i][j]==8 && i==0) 
            {
                hr = grid2D->numcp_global_2D[nx-1][j] - grid2D->numcp_global_2D[i][j];
                nnz=nnz+2;

                val[ii].push_back(1.0);
                row[ii].push_back(ii);
                val[ii+hr].push_back(-1.0);
                row[ii+hr].push_back(ii);

                ii++;
            }

            // periodic boudnary points at right boudary in x direction
            else if ( grid2D->bndr_global_2D[i][j] == 8 && i == nx-1 ) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[1][j];
                nnz=nnz+5;

                val[ii].push_back(-4.0);
                row[ii].push_back(ii);
                val[ii-hl].push_back(1.0);
                row[ii-hl].push_back(ii);
                val[ii-1].push_back(1.0);
                row[ii-1].push_back(ii);
                val[ii-hr].push_back(1.0);
                row[ii-hr].push_back(ii);
                val[ii+1].push_back(1.0);
                row[ii+1].push_back(ii);

                ii++;
            }
        }
    }

    // convert the temp "val row col" to A (compressed column format, i.e. Harwell-Boeing format)
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(grid2D->ncp+1)) ) ABORT("Malloc fails for xa[].");

    i_ncp=0;
    i_nnz=0;
    nz_col=0;
    i_val=0;

    /* old algorithm
    // scan colomn (total ncp column)
    for( i_ncp=0; i_ncp<grid2D->ncp; i_ncp++){
        // is the first not zero element in the current column
        nz_col=0;
        for( i_nnz=0; i_nnz<nnz; i_nnz++)
        {
            if(col[i_nnz] == i_ncp)
            {
                a[i_val] = val[i_nnz];
                asub[i_val] = row[i_nnz];
                if ( nz_col == 0) 
                {
                    xa[i_ncp] = i_val;
                    nz_col = 1;
                }
                i_val++;
            }
        }
    }//>>>end scan
    */

    // new algorithm
    for(int i_col = 0; i_col < grid2D->ncp; i_col++)
    {
        nz_col=0;
        for(int i_row = 0; i_row < val[i_col].size(); i_row++)
        {
            a[i_val]    = val[i_col][i_row];
            asub[i_val] = row[i_col][i_row];
            if(nz_col == 0)
            {
                xa[i_col] = i_val;
                nz_col    = 1;
            }
            i_val++;
        }
    }
    xa[grid2D->ncp]=nnz;
    cout<<"A and b have been structrued"<<endl;

    // ================================ call superlu methods to factorize LU ==================
    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    equil = YES;
    u     = 1.0;
    trans = NOTRANS;

    set_default_options(&options);
    options.Equil = equil;
    options.DiagPivotThresh = u;
    options.Trans = trans;

    m = grid2D->ncp;
    n = grid2D->ncp;
    nrhs = 1;
    if ( !(rhsb = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsx = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);

    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for berr[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* ONLY PERFORM THE LU DECOMPOSITION */
    B.ncol = 0;  /* Indicate not to solve the system */

    cout<<"begin factorizing......"<<endl;

    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

    cout<<"end factorizing......"<<endl;

    printf("LU factorization: dgssvx() returns info %d\n", info);
    StatPrint(&stat);
    StatFree(&stat);
}


void EF_Solver2D_SLU::solve_SLU(Field* rho, Field* phi)
{

    Field2D* rho2D = static_cast<Field2D*>(rho);
    Field2D* phi2D = static_cast<Field2D*>(phi);

    //>>>convert Field2D rho to SuperLU right hand side matrix
    int ii;
    ii = 0;
    for ( int i=0; i<nx; i++)
    {
      for ( int j=0; j<ny; j++) {
        if ( grid2D->bndr_global_2D[i][j] == 0 ) {
          rhsb[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 1) {
          rhsb[ii] = grid2D->bndrVal_global_2D[i][j];
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == 0 || i == 0 )) {
          rhsb[ii] = 0.0;
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == ny-1 || i == nx-1 )) {
          rhsb[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else {
        }

      }
    }//>>>end convert



    /* ------------------------------------------------------------
       NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF A.
       ------------------------------------------------------------*/
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
    B.ncol = nrhs;  /* Set the number of right-hand side */

    /* Initialize the statistics variables. */
    StatInit(&stat);
    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

    //printf("Triangular solve: dgssvx() returns info %d\n", info);

    //>>>convert SuperLU solution X to Field2D phi
    ii=0;
    for ( int i=0; i<nx; i++)
    {
        for ( int j=0; j<ny; j++)
        {
          if ( grid2D->bndr_global_2D[i][j] == 0 || grid2D->bndr_global_2D[i][j] == 1
          || grid2D->bndr_global_2D[i][j] == 2 || grid2D->bndr_global_2D[i][j] == 8) {
            (*phi2D)(i,j) = rhsx[ii];
            ii++;
          }

          if(grid2D->bndr_global_2D[i][j] == 5) {
              (*phi2D)(i,j) = grid2D->bndrVal_global_2D[i][j];
          }

        }//>>>end convert
    }

    StatFree(&stat);


}


void EF_Solver2D_SLU::solve_Exy(Field* phi, Field* Ex, Field* Ey)
{
    Field2D* phi2D = static_cast<Field2D*>(phi);
    Field2D* Ex2D = static_cast<Field2D*>(Ex);
    Field2D* Ey2D = static_cast<Field2D*>(Ey);


    for(int j = 0; j < ny; j++)
    {
        for(int i = 1; i < nx-1; i++)
        {
            (*Ex2D)(i,j) = - ((*phi2D)(i+1,j) - (*phi2D)(i-1,j)) / (2.0*dx);
        }

        (*Ex2D)(0,j) = -(-3.0 * (*phi2D)(0,j) + 4.0 * (*phi2D)(1,j) - (*phi2D)(2,j)) / (2.0*dx);
        (*Ex2D)(nx-1,j) = -((*phi2D)(nx-3,j) - 4.0 * (*phi2D)(nx-2,j) + 3.0 * (*phi2D)(nx-1,j)) / (2.0*dx);
    }


    for(int i = 0; i < nx; i++)
    {
        for(int j = 1; j < ny-1; j++)
        {
            (*Ey2D)(i,j) = - ((*phi2D)(i,j+1) - (*phi2D)(i,j-1)) / (2.0*dy);
        }

        (*Ey2D)(i,0) = - (-3.0 * (*phi2D)(i,0) + 4.0 * (*phi2D)(i,1) - (*phi2D)(i,2)) / (2.0*dy);
        (*Ey2D)(i,ny-1) = - ((*phi2D)(i,ny-3) - 4.0 * (*phi2D)(i,ny-2) + 3.0 * (*phi2D)(i,ny-1)) / (2.0*dy);
    }


}


void EF_Solver2D_SLU::finishSLU()
{
    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork == 0 ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }

}

#endif // for SuperLU_type == serial