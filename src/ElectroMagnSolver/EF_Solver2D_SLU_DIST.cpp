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

    grid2D = static_cast<Grid2D*>(grid);
    //grid2D = new Grid2D(params);
    if(smpi2D->isMaster()){
        //grid2D = static_cast<Grid2D*>(grid);
        initSLU();
    }
    //initSLU();
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
    smpi2D->gatherRho(rho2D_global, rho2D);

    if(smpi2D->isMaster()){
        solve_SLU(rho2D_global, phi2D_global);
        //phi2D_global->put_to(0.0);
        solve_Exy(phi2D_global, Ex2D_global, Ey2D_global);
    }

    //Ex2D_global->put_to(0.0);
    //Ey2D_global->put_to(0.0);


    smpi2D->barrier();
    smpi2D->scatterField(Ex2D_global, Ex2D);
    smpi2D->scatterField(Ey2D_global, Ey2D);
}


void EF_Solver2D_SLU_DIST::initSLU()
{
    double* val;
    double* b;
    int* row;
    int* col;

    //>>>structure the matrix A
    val = new double[grid2D->ncp*5];
    row = new int[grid2D->ncp*5];
    col = new int[grid2D->ncp*5];
    b   = new double[grid2D->ncp];

    for(int i = 0; i < grid2D->ncp*5; i++)
    {
        val[i] = 0.0;
    }

    //unique_ptr<double[]> val(new double[grid2D->ncp*5]);
    //unique_ptr<double[]> b(new double[grid2D->ncp]);
    //unique_ptr<int[]> row(new int[grid2D->ncp*5]);
    //unique_ptr<int[]> col(new int[grid2D->ncp*5]);

    int i,j,k,ii,ll,kk,v,hu,hd,hr,hl,i_ncp,i_nnz,nz_col,i_val;

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
            if(grid2D->bndr_global_2D[i][j]==0) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i+1][j] - grid2D->numcp_global_2D[i][j];
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii+hr]=1.0; b[ii+1]=1.0;
                nnz=nnz+5;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                      val[v] = b [k];
                      row[v] = ii;
                      col[v] = k;
                      v++;
                    }
                }
                ii++;
            }

            // Dirchlet boudnary points
            else if(grid2D->bndr_global_2D[i][j]==1) {
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=1.0;
                nnz++;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                ii++;
            }

            // periodic boudnary points
            else if( grid2D->bndr_global_2D[i][j]==8 && j==0) {
                hu = grid2D->numcp_global_2D[i][ny-1] - grid2D->numcp_global_2D[i][j];
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=1.0; b[ii+hu]=-1.0;
                nnz=nnz+2;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                ii++;
            }

            // periodic boudnary points
            else if ( grid2D->bndr_global_2D[i][j] == 8 && j == ny-1 ) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i+1][j] - grid2D->numcp_global_2D[i][j];
                hd = grid2D->numcp_global_2D[i][ny-1] - grid2D->numcp_global_2D[i][1];
                for(k=0; k<grid2D->ncp; k++){
                b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii+hr]=1.0; b[ii-hd]=1.0;
                nnz=nnz+5;
                for(k=0; k<grid2D->ncp; k++){
                if(b[k] != 0.0){
                    if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                  val[v] = b [k];
                  row[v] = ii;
                  col[v] = k;
                  v++;
                }
                }
                ii++;
            }

            // periodic boudnary points
            else if( grid2D->bndr_global_2D[i][j]==8 && i==0) {
                hr = grid2D->numcp_global_2D[nx-1][j] - grid2D->numcp_global_2D[i][j];
                for(k=0; k<grid2D->ncp; k++){
                    b[k]=0.0;
                }
                b[ii]=1.0; b[ii+hr]=-1.0;
                nnz=nnz+2;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                        val[v] = b [k];
                        row[v] = ii;
                        col[v] = k;
                        v++;
                    }
                }
                ii++;
            }

            // periodic boudnary points
            else if ( grid2D->bndr_global_2D[i][j] == 8 && i == nx-1 ) {
                hl = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[i-1][j];
                hr = grid2D->numcp_global_2D[i][j] - grid2D->numcp_global_2D[1][j];
                for(k=0; k<grid2D->ncp; k++){
                b[k]=0.0;
                }
                b[ii]=-4.0; b[ii-hl]=1.0; b[ii-1]=1.0; b[ii-hr]=1.0; b[ii+1]=1.0;
                nnz=nnz+5;
                for(k=0; k<grid2D->ncp; k++){
                    if(b[k] != 0.0){
                        if(v>=grid2D->ncp*5) cout<<"error"<<v<<endl;
                      val[v] = b [k];
                      row[v] = ii;
                      col[v] = k;
                      v++;
                    }
                }
                ii++;
            }


        }
    }

    //for(int i=0; i<grid2D->ncp*5;i++) cout<<i<<" "<<val[i]<<endl;


    //>>>convert the temp "val row col" to A (compressed column format, i.e. Harwell-Boeing format)
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

    //>>>scan colomn (total ncp column)
    for( i_ncp=0; i_ncp<grid2D->ncp; i_ncp++){
      //>>>is the first not zero element in the current column
      nz_col=0;
      for( i_nnz=0; i_nnz<nnz; i_nnz++){
        if(col[i_nnz] == i_ncp){
          a[i_val] = val[i_nnz];
          asub[i_val] = row[i_nnz];
          if ( nz_col == 0) {
            xa[i_ncp] = i_val;
            nz_col = 1;
          }
          i_val++;
        }
      }
    }//>>>end scan
    cout<<"adrrr "<<val[0]<<" "<<val[grid2D->ncp*5-2]<<" "<<val[grid2D->ncp*5-1]<<endl;
    delete[] val;
    delete[] b;
    delete[] row;
    delete[] col;

    xa[grid2D->ncp]=nnz;
    cout<<"ncp: "<<grid2D->ncp<<" "<<nnz<<endl;

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
        return;
    }

    /* Create compressed column matrix for A. */
    dCreate_CompCol_Matrix_dist(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* Generate the exact solution and compute the right-hand side. 
       The right-hand-side B and the result X are both stored in b */
    if ( !(b = doubleMalloc_dist(m * nrhs)) ) ABORT("Malloc fails for b[]");
    if ( !(xtrue = doubleMalloc_dist(n*nrhs)) ) ABORT("Malloc fails for xtrue[]");
    *trans = 'N';
    ldx = n;
    ldb = m;
    dGenXtrue_dist(n, nrhs, xtrue, ldx);
    dFillRHS_dist(trans, nrhs, xtrue, ldx, &A, b, ldb);

    set_default_options_dist(&options);

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

    /* Call the linear equation solver: factorize and solve. */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid, &LUstruct, berr, &stat, &info);

    /* Check the accuracy of the solution. */
    if ( !iam ) 
    {
	    dinf_norm_error_dist(n, nrhs, b, ldb, xtrue, ldx, &grid);
    }

    PStatPrint(&options, &stat, &grid);        /* Print the statistics. */
    PStatFree(&stat);
}


void EF_Solver2D_SLU_DIST::solve_SLU(Field* rho, Field* phi){

    Field2D* rho2D = static_cast<Field2D*>(rho);
    Field2D* phi2D = static_cast<Field2D*>(phi);

    //>>>convert Field2D rho to SuperLU right hand side matrix
    int ii;
    ii = 0;
    for ( int i=0; i<nx; i++)
    {
      for ( int j=0; j<ny; j++) {
        if ( grid2D->bndr_global_2D[i][j] == 0 ) {
          b[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 1) {
          b[ii] = grid2D->bndrVal_global_2D[i][j];
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == 0 || i == 0 )) {
          b[ii] = 0.0;
          ii++;
        }
        else if ( grid2D->bndr_global_2D[i][j] == 8 && ( j == ny-1 || i == nx-1 )) {
          b[ii] = - dxy * const_ephi0_inv * (*rho2D)(i,j);
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

    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid, &LUstruct, berr, &stat, &info);

    //printf("Triangular solve: dgssvx() returns info %d\n", info);

    //>>>convert SuperLU solution X to Field2D phi
    ii=0;
    for ( int i=0; i<nx; i++)
    {
        for ( int j=0; j<ny; j++)
        {
          if ( grid2D->bndr_global_2D[i][j] == 0 || grid2D->bndr_global_2D[i][j] == 1
          || grid2D->bndr_global_2D[i][j] == 2 || grid2D->bndr_global_2D[i][j] == 8) {
            (*phi2D)(i,j) = b[ii];
            ii++;
          }

          if(grid2D->bndr_global_2D[i][j] == 5) {
              (*phi2D)(i,j) = grid2D->bndrVal_global_2D[i][j];
          }

        }//>>>end convert
    }

   /* Check the accuracy of the solution. */
    if ( !iam ) 
    {
	    printf("Solve the system with a different B.\n");
	    dinf_norm_error_dist(n, nrhs, b, ldb, xtrue, ldx, &grid);
    }

    /* Print the statistics. */
    PStatPrint(&options, &stat, &grid);
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
    SUPERLU_FREE(b);
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