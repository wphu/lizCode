#ifdef SuperLU_mpi

#include "EF_Solver3D_SLU_DIST.h"

#include "ElectroMagn.h"
#include "Field3D.h"
#include <memory>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


EF_Solver3D_SLU_DIST::EF_Solver3D_SLU_DIST(PicParams& params, Grid* grid, SmileiMPI* smpi):
Solver3D(params)
{
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);


    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    dxx = dx * dx;

    nprow = params.number_of_procs[0];
    npcol = params.number_of_procs[1];

    // Now ensure that one compute node is used to solve SuperLU
    if(nprow * npcol > 20)
    {
        nprow = 2;
        npcol = 2;
    }

    grid3D = static_cast<Grid3D*>(grid);
    initSLU();
}


EF_Solver3D_SLU_DIST::~EF_Solver3D_SLU_DIST()
{
}

void EF_Solver3D_SLU_DIST::operator() ( ElectroMagn* fields )
{

}



void EF_Solver3D_SLU_DIST::operator() ( ElectroMagn* fields , SmileiMPI* smpi)
{
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);
    // Static-cast of the fields
    Field3D* Ex3D = static_cast<Field3D*>(fields->Ex_);
    Field3D* Ey3D = static_cast<Field3D*>(fields->Ey_);
    Field3D* Ez3D = static_cast<Field3D*>(fields->Ez_);

    Field3D* rho3D           = static_cast<Field3D*>(fields->rho_);
    Field3D* rho3D_global    = static_cast<Field3D*>(fields->rho_global);
    Field3D* phi3D_global    = static_cast<Field3D*>(fields->phi_global);
    Field3D* Ex3D_global    = static_cast<Field3D*>(fields->Ex_global);
    Field3D* Ey3D_global    = static_cast<Field3D*>(fields->Ey_global);
    Field3D* Ez3D_global    = static_cast<Field3D*>(fields->Ez_global);

    smpi3D->barrier();
    smpi3D->gather_rho_all(rho3D_global, rho3D);

    solve_SLU(rho3D_global, phi3D_global);
    solve_Exyz(phi3D_global, Ex3D_global, Ey3D_global, Ez3D_global);

    smpi3D->barrier();
    smpi3D->scatterField(Ex3D_global, Ex3D);
    smpi3D->scatterField(Ey3D_global, Ey3D);
    smpi3D->scatterField(Ez3D_global, Ez3D);
}


void EF_Solver3D_SLU_DIST::initSLU()
{
    vector< vector<double> > val;
    vector< vector<int> >    row;

    int i,j,k,ii,ll,kk,v,hu,hd,hr,hl,hi,ho,i_ncp,i_nnz,nz_col,i_val;

    val.resize(grid3D->ncp);
    row.resize(grid3D->ncp);

    nnz=0;
    ii=0;
    v=0;
    nx = grid3D->nx;
    ny = grid3D->ny;
    nz = grid3D->nz;

    MESSAGE("Begining structure temporary A ==============");

    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            for(k=0; k<nz; k++)
            {
                // normal points in the calculation region
                if(grid3D->bndr_global_3D(i,j,k)==0) 
                {
                    // order: west(hl), east(hr), north(hd), south(hu), bottom(hi), up(ho)
                    hl = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i-1,j,k);
                    hr = grid3D->numcp_global_3D(i+1,j,k) - grid3D->numcp_global_3D(i,j,k);
                    hd = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j-1,k);
                    hu = grid3D->numcp_global_3D(i,j+1,k) - grid3D->numcp_global_3D(i,j,k);
                    hi = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j,k-1);
                    ho = grid3D->numcp_global_3D(i,j,k+1) - grid3D->numcp_global_3D(i,j,k);

                    nnz = nnz + 7;
                    
                    val[ii].push_back(-6.0);
                    row[ii].push_back(ii);
                    val[ii-hl].push_back(1.0);
                    row[ii-hl].push_back(ii);
                    val[ii+hr].push_back(1.0);
                    row[ii+hr].push_back(ii);
                    val[ii-hd].push_back(1.0);
                    row[ii-hd].push_back(ii);
                    val[ii+hu].push_back(1.0);
                    row[ii+hu].push_back(ii);
                    val[ii-hi].push_back(1.0);
                    row[ii-hi].push_back(ii);
                    val[ii+ho].push_back(1.0);
                    row[ii+ho].push_back(ii);

                    ii++;
                }

                // Dirchlet boudnary points
                else if(grid3D->bndr_global_3D(i,j,k)==1) 
                {
                    nnz++;

                    val[ii].push_back(1.0);
                    row[ii].push_back(ii);

                    ii++;
                }

 
                // periodic boudnary points at left boudary in x direction
                else if( grid3D->bndr_global_3D(i,j,k) == 8 && i == 0) 
                {
                    hr = grid3D->numcp_global_3D(nx-1,j,k) - grid3D->numcp_global_3D(i,j,k);
                    nnz = nnz + 2;

                    val[ii].push_back(1.0);
                    row[ii].push_back(ii);
                    val[ii+hr].push_back(-1.0);
                    row[ii+hr].push_back(ii);

                    ii++;
                }

                // periodic boudnary points at right boudary in x direction
                else if ( grid3D->bndr_global_3D(i,j,k) == 8 && i == nx-1 ) {
                    if(j == 0 || j == ny - 1)
                    {
                        hl = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i-1,j,k);
                        hr = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(1,j,k);

                        nnz = nnz + 3;

                        val[ii].push_back(-2.0);
                        row[ii].push_back(ii);
                        val[ii-hl].push_back(1.0);
                        row[ii-hl].push_back(ii);
                        val[ii-hr].push_back(1.0);
                        row[ii-hr].push_back(ii);

                        ii++;
                    }
                    else
                    {
                        hl = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i-1,j,k);
                        hr = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(1,j,k);
                        hd = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j-1,k);
                        hu = grid3D->numcp_global_3D(i,j+1,k) - grid3D->numcp_global_3D(i,j,k);
                        hi = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j,k-1);
                        ho = grid3D->numcp_global_3D(i,j,k+1) - grid3D->numcp_global_3D(i,j,k);

                        nnz = nnz + 7;

                        val[ii].push_back(-6.0);
                        row[ii].push_back(ii);
                        val[ii-hl].push_back(1.0);
                        row[ii-hl].push_back(ii);
                        val[ii-hr].push_back(1.0);
                        row[ii-hr].push_back(ii);
                        val[ii-hd].push_back(1.0);
                        row[ii-hd].push_back(ii);
                        val[ii+hu].push_back(1.0);
                        row[ii+hu].push_back(ii);
                        val[ii-hi].push_back(1.0);
                        row[ii-hi].push_back(ii);
                        val[ii+ho].push_back(1.0);
                        row[ii+ho].push_back(ii);

                        ii++;
                    }
                }

               // periodic boudnary points at lowwer boudary in y direction
                else if( grid3D->bndr_global_3D(i,j,k)==8 && j==0) 
                {
                    hu = grid3D->numcp_global_3D(i,j+1,k) - grid3D->numcp_global_3D(i,j,k);

                    nnz = nnz + 2;

                    val[ii].push_back(1.0);
                    row[ii].push_back(ii);
                    val[ii+hu].push_back(-1.0);
                    row[ii+hu].push_back(ii);

                    ii++;
                }

                // periodic boudnary points at upper boudary in y direction
                else if ( grid3D->bndr_global_3D(i,j,k) == 8 && j == ny-1 ) 
                {
                    hl = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i-1,j,k);
                    hr = grid3D->numcp_global_3D(i+1,j,k) - grid3D->numcp_global_3D(i,j,k);
                    hd = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j-1,k);
                    hu = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,1,k);
                    hi = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j,k-1);
                    ho = grid3D->numcp_global_3D(i,j,k+1) - grid3D->numcp_global_3D(i,j,k);

                    nnz = nnz + 7;

                    val[ii].push_back(-6.0);
                    row[ii].push_back(ii);
                    val[ii-hl].push_back(1.0);
                    row[ii-hl].push_back(ii);
                    val[ii+hr].push_back(1.0);
                    row[ii+hr].push_back(ii);
                    val[ii-hd].push_back(1.0);
                    row[ii-hd].push_back(ii);
                    val[ii-hu].push_back(1.0);
                    row[ii-hu].push_back(ii);
                    val[ii-hi].push_back(1.0);
                    row[ii-hi].push_back(ii);
                    val[ii+ho].push_back(1.0);
                    row[ii+ho].push_back(ii);

                    ii++;
                }

               // periodic boudnary points at lowwer boudary in z direction
                else if( grid3D->bndr_global_3D(i,j,k)==8 && k==0) 
                {
                    ho = grid3D->numcp_global_3D(i,j,k+1) - grid3D->numcp_global_3D(i,j,k);
                    nnz = nnz + 2;

                    val[ii].push_back(1.0);
                    row[ii].push_back(ii);
                    val[ii+ho].push_back(-1.0);
                    row[ii+ho].push_back(ii);

                    ii++;
                }

                // periodic boudnary points at upper boudary in z direction
                else if ( grid3D->bndr_global_3D(i,j,k) == 8 && k == nz-1 ) 
                {
                    hl = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i-1,j,k);
                    hr = grid3D->numcp_global_3D(i+1,j,k) - grid3D->numcp_global_3D(i,j,k);
                    hd = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j-1,k);
                    hu = grid3D->numcp_global_3D(i,j+1,k) - grid3D->numcp_global_3D(i,j,k);
                    hi = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j,k-1);
                    ho = grid3D->numcp_global_3D(i,j,k) - grid3D->numcp_global_3D(i,j,1);

                    nnz = nnz + 7;

                    val[ii].push_back(-6.0);
                    row[ii].push_back(ii);
                    val[ii-hl].push_back(1.0);
                    row[ii-hl].push_back(ii);
                    val[ii+hr].push_back(1.0);
                    row[ii+hr].push_back(ii);
                    val[ii-hd].push_back(1.0);
                    row[ii-hd].push_back(ii);
                    val[ii-hu].push_back(1.0);
                    row[ii-hu].push_back(ii);
                    val[ii-hi].push_back(1.0);
                    row[ii-hi].push_back(ii);
                    val[ii-ho].push_back(1.0);
                    row[ii-ho].push_back(ii);

                    ii++;
                }
                else
                {
                    
                }

            }
        }

    }

    MESSAGE("Temporary A has been structrued");

    // convert the temp "val row col" to A (compressed column format, i.e. Harwell-Boeing format)
    //a = new double[nnz];
    //asub = new int[nnz];
    //xa = new int[grid3D->ncp+1];
    if ( !(a = doubleMalloc_dist(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc_dist(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc_dist(grid3D->ncp+1)) ) ABORT("Malloc fails for xa[].");



    i_ncp=0;
    i_nnz=0;
    nz_col=0;
    i_val=0;

    // new algorithm
    for(int i_col = 0; i_col < grid3D->ncp; i_col++)
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

    xa[grid3D->ncp]=nnz;

    m = grid3D->ncp;
    n = grid3D->ncp;
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


void EF_Solver3D_SLU_DIST::solve_SLU(Field* rho, Field* phi)
{
    if ( iam >= nprow * npcol )
    {
        //cout<<"iam >= nprow * npcol: "<<iam<<endl;
        return;
    }

    Field3D* rho3D = static_cast<Field3D*>(rho);
    Field3D* phi3D = static_cast<Field3D*>(phi);

    //>>>convert Field3D rho to SuperLU right hand side matrix
    int ii;
    ii = 0;
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++) 
        {
            for(int k = 0; k < nz; k++)
            {
                if(grid3D->bndr_global_3D(i,j,k) == 0 ) 
                {
                    b1[ii] = - dxx * const_ephi0_inv * (*rho3D)(i,j,k);
                    ii++;
                }
                else if(grid3D->bndr_global_3D(i,j,k) == 1) 
                {
                    b1[ii] = grid3D->bndrVal_global_3D(i,j,k);
                    ii++;
                }
                else if(grid3D->bndr_global_3D(i,j,k) == 8 && ( i == 0 || j == 0 || k == 0)) 
                {
                    b1[ii] = 0.0;
                    ii++;
                }
                else if(grid3D->bndr_global_3D(i,j,k) == 8 && ( i == nx - 1 || j == ny - 1 || k == nz - 1)) 
                {
                    b1[ii] = - dxx * const_ephi0_inv * (*rho3D)(i,j,k);
                    ii++;
                }
                else 
                {
                }
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

    //>>>convert SuperLU solution X to Field3D phi
    ii=0;
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                if (grid3D->bndr_global_3D(i,j,k) == 0 || grid3D->bndr_global_3D(i,j,k) == 1
                 || grid3D->bndr_global_3D(i,j,k) == 2 || grid3D->bndr_global_3D(i,j,k) == 8) 
                {
                    (*phi3D)(i,j,k) = b1[ii];
                    ii++;
                }

                if(grid3D->bndr_global_3D(i,j,k) == 5) 
                {
                    (*phi3D)(i,j,k) = grid3D->bndrVal_global_3D(i,j,k);
                }
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


void EF_Solver3D_SLU_DIST::solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez)
{
    Field3D* phi3D = static_cast<Field3D*>(phi);
    Field3D* Ex3D = static_cast<Field3D*>(Ex);
    Field3D* Ey3D = static_cast<Field3D*>(Ey);
    Field3D* Ez3D = static_cast<Field3D*>(Ez);

    for(int i = 1; i < nx - 1; i++)
    {
        for(int j = 1; j < ny - 1; j++)
        {
            for(int k = 1; k < nz - 1; k++)
            {
                (*Ex3D)(i,j,k) = - ((*phi3D)(i+1,j,k) - (*phi3D)(i-1,j,k)) / (2.0*dx);
                (*Ey3D)(i,j,k) = - ((*phi3D)(i,j+1,k) - (*phi3D)(i,j-1,k)) / (2.0*dx);
                (*Ez3D)(i,j,k) = - ((*phi3D)(i,j,k+1) - (*phi3D)(i,j,k-1)) / (2.0*dx);
            }
        }
    }

    for(int j = 1; j < ny - 1; j++)
    {
        for(int k = 1; k < nz - 1; k++)
        {
            (*Ex3D)(0,j,k) = -(-3.0 * (*phi3D)(0,j,k) + 4.0 * (*phi3D)(1,j,k) - (*phi3D)(2,j,k)) / (2.0*dx);
            (*Ex3D)(nx-1,j,k) = -((*phi3D)(nx-3,j,k) - 4.0 * (*phi3D)(nx-2,j,k) + 3.0 * (*phi3D)(nx-1,j,k)) / (2.0*dx);

            (*Ey3D)(0,j,k) = - ((*phi3D)(0,j+1,k) - (*phi3D)(0,j-1,k)) / (2.0*dx);
            (*Ey3D)(nx-1,j,k) = - ((*phi3D)(nx-1,j+1,k) - (*phi3D)(nx-1,j-1,k)) / (2.0*dx);
            (*Ez3D)(0,j,k) = - ((*phi3D)(0,j,k+1) - (*phi3D)(0,j,k-1)) / (2.0*dx);
            (*Ez3D)(nx-1,j,k) = - ((*phi3D)(nx-1,j,k+1) - (*phi3D)(nx-1,j,k-1)) / (2.0*dx);
        }
    }

    for(int i = 1; i < nx - 1; i++)
    {
        for(int k = 1; k < nz - 1; k++)
        {
            (*Ey3D)(i,0,k) = -(-3.0 * (*phi3D)(i,0,k) + 4.0 * (*phi3D)(i,1,k) - (*phi3D)(i,2,k)) / (2.0*dx);
            (*Ey3D)(i,ny-1,k) = -((*phi3D)(i,ny-3,k) - 4.0 * (*phi3D)(i,ny-2,k) + 3.0 * (*phi3D)(i,ny-1,k)) / (2.0*dx);

            (*Ex3D)(i,0,k) = - ((*phi3D)(i+1,0,k) - (*phi3D)(i-1,0,k)) / (2.0*dx);
            (*Ex3D)(i,ny-1,k) = - ((*phi3D)(i+1,ny-1,k) - (*phi3D)(i-1,ny-1,k)) / (2.0*dx);
            (*Ez3D)(i,0,k) = - ((*phi3D)(i,0,k+1) - (*phi3D)(i,0,k-1)) / (2.0*dx);
            (*Ez3D)(i,ny-1,k) = - ((*phi3D)(i,ny-1,k+1) - (*phi3D)(i,ny-1,k-1)) / (2.0*dx);
        }
    }

    for(int i = 1; i < nx - 1; i++)
    {
        for(int j = 1; j < ny - 1; j++)
        {
            (*Ez3D)(i,j,0) = -(-3.0 * (*phi3D)(i,j,0) + 4.0 * (*phi3D)(i,j,1) - (*phi3D)(i,j,2)) / (2.0*dx);
            (*Ez3D)(i,j,nz-1) = -((*phi3D)(i,j,nz-3) - 4.0 * (*phi3D)(i,j,nz-2) + 3.0 * (*phi3D)(i,j,nz-1)) / (2.0*dx);

            (*Ex3D)(i,j,0) = - ((*phi3D)(i+1,j,0) - (*phi3D)(i-1,j,0)) / (2.0*dx);
            (*Ex3D)(i,j,nz-1) = - ((*phi3D)(i+1,j,nz-1) - (*phi3D)(i-1,j,nz-1)) / (2.0*dx);
            (*Ey3D)(i,j,0) = - ((*phi3D)(i,j+1,0) - (*phi3D)(i,j-1,0)) / (2.0*dx);
            (*Ey3D)(i,j,nz-1) = - ((*phi3D)(i,j+1,nz-1) - (*phi3D)(i,j-1,nz-1)) / (2.0*dx);
        }
    }

}




void EF_Solver3D_SLU_DIST::finishSLU()
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