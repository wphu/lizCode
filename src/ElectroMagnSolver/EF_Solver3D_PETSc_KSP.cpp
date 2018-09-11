#ifdef petsc

#include "EF_Solver3D_PETSc_KSP.h"

#include "ElectroMagn.h"
#include "Field3D.h"
#include <memory>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


EF_Solver3D_PETSc_KSP::EF_Solver3D_PETSc_KSP(PicParams& params, Grid* grid, SmileiMPI* smpi):
Solver3D(params)
{
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);

    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    dxx = dx * dx;
    petsc_ksp_tolerance_rtol = params.petsc_ksp_tolerance_rtol;

    grid3D = static_cast<Grid3D*>(grid);

    MPI_Comm_rank( MPI_COMM_WORLD, &petsc_ksp_mpi_rank );

    // create sub communicator ofr petsc_ksp parallel solver using params.petsc_ksp_process_number
    petsc_ksp_mpi_size = params.petsc_ksp_process_number;
    MPI_Group group_world;
    MPI_Group petsc_ksp_group;
    MPI_Comm  petsc_ksp_comm;
    int *process_ranks;

    // make a list of processes in the new communicator
    process_ranks = new int[petsc_ksp_mpi_size];
    for(int i = 0; i < petsc_ksp_mpi_size; i++)
    {
        process_ranks[i] = i;
    }
    // get the group under MPI_COMM_WORLD
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    // create the new group
    MPI_Group_incl(group_world, petsc_ksp_mpi_size, process_ranks, &petsc_ksp_group);
    // create the new communicator
    MPI_Comm_create(MPI_COMM_WORLD, petsc_ksp_group, &petsc_ksp_comm);
    
    PETSC_COMM_WORLD = petsc_ksp_comm;

    if( is_in_pestc_ksp_mpicomm() )
    {
        init_PETSc_KSP();
    }

}


EF_Solver3D_PETSc_KSP::~EF_Solver3D_PETSc_KSP()
{
    VecScatterDestroy(&ctx);
}

void EF_Solver3D_PETSc_KSP::operator() ( ElectroMagn* fields )
{

}



void EF_Solver3D_PETSc_KSP::operator() ( ElectroMagn* fields , SmileiMPI* smpi)
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
    smpi3D->gatherAllRho(rho3D_global, rho3D);
    if( is_in_pestc_ksp_mpicomm() )
    {
        solve_PETSc_KSP(rho3D_global, phi3D_global);
        solve_Exyz(phi3D_global, Ex3D_global, Ey3D_global, Ez3D_global);
    }
    

    smpi3D->barrier();
    smpi3D->scatterField(Ex3D_global, Ex3D);
    smpi3D->scatterField(Ey3D_global, Ey3D);
    smpi3D->scatterField(Ez3D_global, Ez3D);
}


void EF_Solver3D_PETSc_KSP::init_PETSc_KSP_test()
{

}


void EF_Solver3D_PETSc_KSP::init_PETSc_KSP()
{
    vector< vector<double> > val;
    vector< vector<int> >    row;

    int i,j,k,ii,ll,kk,v,hu,hd,hr,hl,hi,ho,i_ncp,i_nnz,nz_col,i_val;

    ncp = grid3D->ncp;
    val.resize(grid3D->ncp);
    row.resize(grid3D->ncp);

    rhsb = new double[ncp];
    rhsx = new double[ncp];
    indices_global = new int[ncp];
    for(int i = 0; i < ncp; i++)
    {
        rhsb[i] = 0.0;
        indices_global[i] = i;
    }

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
                    hu = grid3D->numcp_global_3D(i,ny-1,k) - grid3D->numcp_global_3D(i,j,k);

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
                    ho = grid3D->numcp_global_3D(i,j,nz-1) - grid3D->numcp_global_3D(i,j,k);
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

    // =========================== PETSc KSP part ===================================
    char help[] =   "Solves a linear system in parallel with KSP.\n\
                            Input parameters include:\n\
                            -random_exact_sol : use a random exact solution vector\n\
                            -view_exact_sol   : write exact solution vector to stdout\n\
                            -m <mesh_x>       : number of mesh points in x-direction\n\
                            -n <mesh_n>       : number of mesh points in y-direction\n\n";
    


    int argc = 0;
    char **args;
    ierr = PetscInitialize(&argc,&args,(char*)0,help);
    if(ierr)
    {
        cout<<"PetscInitialize failed !!!"<<endl;
        ERROR("PetscInitialize failed !!!");
    }

    //MESSAGE("1111");
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Compute the matrix and right-hand-side vector that define
            the linear system, Ax = b.
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
        Create parallel matrix, specifying only its global dimensions.
        When using MatCreate(), the matrix format can be specified at
        runtime. Also, the parallel partitioning of the matrix is
        determined by PETSc at runtime.

        Performance tuning note:  For problems of substantial size,
        preallocation of matrix memory is crucial for attaining good
        performance. See the matrix chapter of the users manual for details.
    */
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRV(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,ncp,ncp);CHKERRV(ierr);
    ierr = MatSetFromOptions(A);CHKERRV(ierr);
    //ierr = MatSetUp(A);CHKERRV(ierr);
    ierr = MatMPIAIJSetPreallocation(A,7,NULL,7,NULL);CHKERRV(ierr);
    ierr = MatSeqAIJSetPreallocation(A,7,NULL);CHKERRV(ierr);
    //ierr = MatSeqSBAIJSetPreallocation(A,ncp,ncp,NULL);CHKERRV(ierr);
    //ierr = MatMPISBAIJSetPreallocation(A,ncp,ncp,NULL,ncp,NULL);CHKERRV(ierr);

    //MESSAGE("22222");
    /*
        Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors.  Determine which
        rows of the matrix are locally owned.
    */
    //ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRV(ierr);

    /*
        Set matrix elements for the 3-D, seven-point stencil in parallel.
        - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
        - Always specify global rows and columns of matrix entries.
        !!! But, now, each processor inserts all elements of global A, for easy ~~
    */
    //ierr = PetscLogStageRegister("Assembly", &stage);CHKERRV(ierr);
    //MESSAGE("ncp "<<ncp);
    for(int j = 0; j < ncp; j++)
    {
        for(int i = 0; i < row[j].size(); i++)
        {
            I = row[j][i];
            J = j;
            PetscScalar val_temp = val[j][i];
            ierr = MatSetValues(A,1,&I,1,&J,&val_temp,INSERT_VALUES);CHKERRV(ierr);
        }
    }

    //MESSAGE("4444");
    /*
        Assemble matrix, using the 2-step process:
        MatAssemblyBegin(), MatAssemblyEnd()
        Computations can be done while messages are in transition
        by placing code between these two statements.
    */
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
    //ierr = PetscLogStagePop();CHKERRV(ierr);

    /*
        Create parallel vectors.
        - We form 1 vector from scratch and then duplicate as needed.
        - When using VecCreate(), VecSetSizes and VecSetFromOptions()
        in this example, we specify only the
        vector's global dimension; the parallel partitioning is determined
        at runtime.
        - When solving a linear system, the vectors and matrices MUST
        be partitioned accordingly.  PETSc automatically generates
        appropriately partitioned matrices and vectors when MatCreate()
        and VecCreate() are used with the same communicator.
        - The user can alternatively specify the local vector and matrix
        dimensions when more sophisticated partitioning is needed
        (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
        below).
    */
    ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRV(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,ncp);CHKERRV(ierr);
    ierr = VecSetFromOptions(u);CHKERRV(ierr);
    ierr = VecDuplicate(u,&b);CHKERRV(ierr);
    ierr = VecDuplicate(b,&x);CHKERRV(ierr);  
    
    //MESSAGE("55555");
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
        Create linear solver context
    */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRV(ierr);

    /*
        Set operators. Here the matrix that defines the linear system
        also serves as the preconditioning matrix.
    */
    ierr = KSPSetOperators(ksp,A,A);CHKERRV(ierr);

    // set KSP solver type: KSPGMRES, KSPCG, KSPCR, KSPCGS, KSPBCGS
    // KSPCG, KSPCR do not give right solution
    ierr = KSPSetType(ksp, KSPGMRES); CHKERRV(ierr);

    // set precondition: PCGAMG, PCILU, PCJACOBI, PCBJACOBI, PCSOR, PCEISENSTAT, PCICC
    //                   PCASM, PCBDDC, PCKSP, PCCOMPOSITE
    // PCSOR is fastest
    //ierr = KSPGetPC(ksp, &pc); CHKERRV(ierr);
    //ierr = PCSetType(pc, PCSOR); CHKERRV(ierr);

    /*
        Set linear solver defaults for this problem (optional).
        - By extracting the KSP and PC contexts from the KSP context,
        we can then directly call any KSP and PC routines to set
        various options.
        - The following two statements are optional; all of these
        parameters could alternatively be specified at runtime via
        KSPSetFromOptions().  All of these defaults can be
        overridden at runtime, as indicated below.
    */
    ierr = KSPSetTolerances(ksp,petsc_ksp_tolerance_rtol,PETSC_DEFAULT,PETSC_DEFAULT,
                            PETSC_DEFAULT);CHKERRV(ierr);

    /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */
    ierr = KSPSetFromOptions(ksp);CHKERRV(ierr);

    VecScatterCreateToAll(x, &ctx, &x_seq);
    //MESSAGE("6666");
}


void EF_Solver3D_PETSc_KSP::solve_PETSc_KSP(Field* rho, Field* phi)
{

    Field3D* rho3D = static_cast<Field3D*>(rho);
    Field3D* phi3D = static_cast<Field3D*>(phi);

    // convert Field3D rho to rhsb
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
                    rhsb[ii] = - dxx * const_ephi0_inv * (*rho3D)(i,j,k);
                    ii++;
                }
                else if(grid3D->bndr_global_3D(i,j,k) == 1) 
                {
                    rhsb[ii] = grid3D->bndrVal_global_3D(i,j,k);
                    ii++;
                }
                else if(grid3D->bndr_global_3D(i,j,k) == 8 && ( i == 0 || j == 0 || k == 0)) 
                {
                    rhsb[ii] = 0.0;
                    ii++;
                }
                else if(grid3D->bndr_global_3D(i,j,k) == 8 && ( i == nx - 1 || j == ny - 1 || k == nz - 1)) 
                {
                    rhsb[ii] = - dxx * const_ephi0_inv * (*rho3D)(i,j,k);
                    ii++;
                }
                else 
                {
                }
            }
        }
    }// end convert

    // ======================= set right hand side: b ======================
    ierr = VecSetValues(b, ncp, indices_global, rhsb, INSERT_VALUES);
    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the linear system
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = KSPSolve(ksp,b,x);CHKERRV(ierr);

    // ================ gather x of all processors to global rhsx ==============
    VecScatterBegin(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD);
    VecGetArray(x_seq,&rhsx);

    /*
        Check the error
    */
    //ierr = VecAXPY(x,-1.0,u);CHKERRV(ierr);
    //ierr = VecNorm(x,NORM_2,&norm);CHKERRV(ierr);
    //ierr = KSPGetIterationNumber(ksp,&its);CHKERRV(ierr);

    /*
        Print convergence information.  PetscPrintf() produces a single
        print statement from all processes that share a communicator.
        An alternative is PetscFPrintf(), which prints to a file.
    */
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);CHKERRV(ierr);




    // convert rhsx to Field3D phi
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
                    (*phi3D)(i,j,k) = rhsx[ii];
                    ii++;
                }

                if(grid3D->bndr_global_3D(i,j,k) == 5) 
                {
                    (*phi3D)(i,j,k) = grid3D->bndrVal_global_3D(i,j,k);
                }
            }


        }// end convert
    }



}


void EF_Solver3D_PETSc_KSP::solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez)
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


void EF_Solver3D_PETSc_KSP::finish_PETSc_KSP()
{

}

#endif // for petsc