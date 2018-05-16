#ifdef SuperLU_mpi

#include "EF_Solver2D_SLU_DIST.h"

#include "ElectroMagn.h"
#include "Field2D.h"
#include <memory>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


EF_Solver2D_SLU_DIST::EF_Solver2D_SLU_DIST(PicParams& params, Grid* grid, SmileiMPI* smpi):
Solver2D(params)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);


    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dxy = dx * dy;

    nprow = params.number_of_procs[0];
    npcol = params.number_of_procs[1];

    // Now ensure that one compute node is used to solve SuperLU
    if(nprow * npcol > 20)
    {
        nprow = 4;
        npcol = 5;
    }

    grid2D = static_cast<Grid2D*>(grid);
    initSLU();
}


EF_Solver2D_SLU_DIST::~EF_Solver2D_SLU_DIST()
{
}

void EF_Solver2D_SLU_DIST::operator() ( ElectroMagn* fields )
{

}



void EF_Solver2D_SLU_DIST::operator() ( ElectroMagn* fields , SmileiMPI* smpi)
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
    smpi2D->gather_rho_all(rho2D_global, rho2D);

    solve_SLU(rho2D_global, phi2D_global);
    solve_Exy(phi2D_global, Ex2D_global, Ey2D_global);

    smpi2D->barrier();
    smpi2D->scatterField(Ex2D_global, Ex2D);
    smpi2D->scatterField(Ey2D_global, Ey2D);
}


void EF_Solver2D_SLU_DIST::initSLU()
{
    vector< vector<double> > val;
    vector< vector<int> >    row;

    int i,j,k,ii,ll,kk,v,hu,hd,hr,hl,i_ncp,i_nnz,nz_col,i_val;

    val.resize(grid2D->ncp*5);
    row.resize(grid2D->ncp*5);

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
    //a = new double[nnz];
    //asub = new int[nnz];
    //xa = new int[grid2D->ncp+1];
    if ( !(a = doubleMalloc_dist(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc_dist(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc_dist(grid2D->ncp+1)) ) ABORT("Malloc fails for xa[].");



    i_ncp=0;
    i_nnz=0;
    nz_col=0;
    i_val=0;

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
    MESSAGE("A and b have been structrued");

    xa[grid2D->ncp]=nnz;

    m = grid2D->ncp;
    n = grid2D->ncp;
    nrhs = 1;

    /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID. 
       ------------------------------------------------------------*/
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);

    /* Bail out if I do not belong in the grid. */
    iam = grid.iam;
    if ( iam >= nprow * npcol )
    {
        //cout<<"iam >= nprow * npcol: "<<iam<<endl;
        return;
    }

    /* Create compressed column matrix for A. */
    dCreate_CompCol_Matrix_dist(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* Generate the exact solution and compute the right-hand side. 
       The right-hand-side B and the result X are both stored in b1 */
    if ( !(b1 = doubleMalloc_dist(m * nrhs)) ) ABORT("Malloc fails for b1[]");
    if ( !(xtrue = doubleMalloc_dist(n*nrhs)) ) ABORT("Malloc fails for xtrue[]");
    *trans = 'N';
    ldx = n;
    ldb = m;
    dGenXtrue_dist(n, nrhs, xtrue, ldx);
    dFillRHS_dist(trans, nrhs, xtrue, ldx, &A, b1, ldb);

    if ( !(berr = doubleMalloc_dist(nrhs)) )
	ABORT("Malloc fails for berr[].");

    set_default_options_dist(&options);
    //options.IterRefine = NO;

    if (!iam) 
    {
        print_sp_ienv_dist(&options);
        print_options_dist(&options);
    }

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(n, &LUstruct);

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    MESSAGE("begin factorizing......");

    /* Call the linear equation solver: factorize and solve. */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b1, ldb, nrhs, &grid, &LUstruct, berr, &stat, &info);
    
    MESSAGE("end factorizing......");
    
    /* Check the accuracy of the solution. */
    if ( !iam ) 
    {
	    dinf_norm_error_dist(n, nrhs, b1, ldb, xtrue, ldx, &grid);
    }

    PStatPrint(&options, &stat, &grid);        /* Print the statistics. */
    PStatFree(&stat);

    MESSAGE("SuperLU_DIST init done!!!");
}


void EF_Solver2D_SLU_DIST::solve_SLU(Field* rho, Field* phi)
{
    if ( iam >= nprow * npcol )
    {
        //cout<<"iam >= nprow * npcol: "<<iam<<endl;
        return;
    }

    Field2D* rho2D = static_cast<Field2D*>(rho);
    Field2D* phi2D = static_cast<Field2D*>(phi);

    //>>>convert Field2D rho to SuperLU right hand side matrix
    int ii;
    ii = 0;
    for ( int i=0; i<nx; i++)
    {
      for ( int j=0; j<ny; j++) {
        if ( grid2D->bndr_global_2D[i][j] == 0 ) {
          b1[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 1) {
          b1[ii] = grid2D->bndrVal_global_2D[i][j];
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == 0 || i == 0 )) {
          b1[ii] = 0.0;
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == ny-1 || i == nx-1 )) {
          b1[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else {
        }

      }
    }//>>>end convert

    /* ------------------------------------------------------------
       NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF A.
       The right-hand-side B and the result X are both stored in b
       ------------------------------------------------------------*/
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
    PStatInit(&stat); /* Initialize the statistics variables. */

    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b1, ldb, nrhs, &grid, &LUstruct, berr, &stat, &info);

    //printf("Triangular solve: dgssvx() returns info %d\n", info);

    //>>>convert SuperLU solution X to Field2D phi
    ii=0;
    for ( int i=0; i<nx; i++)
    {
        for ( int j=0; j<ny; j++)
        {
          if ( grid2D->bndr_global_2D[i][j] == 0 || grid2D->bndr_global_2D[i][j] == 1
          || grid2D->bndr_global_2D[i][j] == 2 || grid2D->bndr_global_2D[i][j] == 8) {
            (*phi2D)(i,j) = b1[ii];
            ii++;
          }

          if(grid2D->bndr_global_2D[i][j] == 5) {
              (*phi2D)(i,j) = grid2D->bndrVal_global_2D[i][j];
          }

        }//>>>end convert
    }

    /* Check the accuracy of the solution. */
    /*
    if ( !iam ) 
    {
	    printf("Solve the system with a different B.\n");
	    dinf_norm_error_dist(n, nrhs, b1, ldb, xtrue, ldx, &grid);
    }
    */

    /* Print the statistics. */
    //PStatPrint(&options, &stat, &grid);

}


void EF_Solver2D_SLU_DIST::solve_Exy(Field* phi, Field* Ex, Field* Ey)
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




void EF_Solver2D_SLU_DIST::finishSLU()
{
   /* ------------------------------------------------------------
       DEALLOCATE STORAGE.
       ------------------------------------------------------------*/
    PStatFree(&stat);
    Destroy_CompCol_Matrix_dist(&A);
    Destroy_LU(n, &grid, &LUstruct);
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct);
    SUPERLU_FREE(b1);
    SUPERLU_FREE(xtrue);
    SUPERLU_FREE(berr);

    /* ------------------------------------------------------------
       RELEASE THE SUPERLU PROCESS GRID.
       ------------------------------------------------------------*/
    superlu_gridexit(&grid);

    /* ------------------------------------------------------------
       TERMINATES THE MPI EXECUTION ENVIRONMENT.
       ------------------------------------------------------------*/
    MPI_Finalize();
}

#endif // for SuperLU_mpi