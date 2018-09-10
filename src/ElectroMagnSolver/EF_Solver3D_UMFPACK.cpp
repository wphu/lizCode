#ifdef umfpack

#include "EF_Solver3D_UMFPACK.h"

#include "ElectroMagn.h"
#include "Field3D.h"
#include <memory>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


EF_Solver3D_UMFPACK::EF_Solver3D_UMFPACK(PicParams& params, Grid* grid, SmileiMPI* smpi):
Solver3D(params)
{
    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);

    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    dxx = dx * dx;

    grid3D = static_cast<Grid3D*>(grid);

    if(smpi3D->isMaster())
    {
        initUMFPACK();
    }


}


EF_Solver3D_UMFPACK::~EF_Solver3D_UMFPACK()
{
}

void EF_Solver3D_UMFPACK::operator() ( ElectroMagn* fields )
{

}



void EF_Solver3D_UMFPACK::operator() ( ElectroMagn* fields , SmileiMPI* smpi)
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
    smpi3D->gatherRho(rho3D_global, rho3D);

    if(smpi3D->isMaster())
    {
        solve_UMFPACK(rho3D_global, phi3D_global);
        solve_Exyz(phi3D_global, Ex3D_global, Ey3D_global, Ez3D_global);
    }

    smpi3D->barrier();
    smpi3D->scatterField(Ex3D_global, Ex3D);
    smpi3D->scatterField(Ey3D_global, Ey3D);
    smpi3D->scatterField(Ez3D_global, Ez3D);
}


void EF_Solver3D_UMFPACK::initUMFPACK_test()
{

}


void EF_Solver3D_UMFPACK::initUMFPACK()
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

    MESSAGE("Temporary A has been structrued, n = "<<grid3D->ncp<<" nnz = "<<nnz);

    n = grid3D->ncp;
    Ap = new SuiteSparse_long[n + 1];
    Ai = new SuiteSparse_long[nnz];
    Ax = new double[nnz];
    rhsb  = new double[n];
    rhsx  = new double[n];
    for(int i = 0; i < n; i++)
    {
        rhsx[i] = 500.0;
    }

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
            Ax[i_val]    = val[i_col][i_row];
            Ai[i_val] = row[i_col][i_row];
            if(nz_col == 0)
            {
                Ap[i_col] = i_val;
                nz_col    = 1;
            }
            i_val++;
        }
    }
    Ap[grid3D->ncp]=nnz;
    vector< vector<double> >().swap(val);
    vector< vector<int> >().swap(row);
    cout<<"A and b have been structrued"<<endl;

    
    cout<<"begin factorizing......"<<endl;

    //get the default control parameters
    umfpack_dl_defaults (Control) ;

    status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status) ;
    printf ("\nSymbolic factorization of A: ") ;
    (void) umfpack_dl_report_symbolic (Symbolic, Control) ;

    status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status) ;
    // print the numeric factorization
    printf ("\nNumeric factorization of A: ") ;
    (void) umfpack_dl_report_numeric (Numeric, Control) ;

    cout<<"end factorizing......"<<endl;

}


void EF_Solver3D_UMFPACK::solve_UMFPACK(Field* rho, Field* phi)
{

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
    }//>>>end convert


    status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, rhsx, rhsb, Numeric, Control, Info);
    umfpack_dl_report_info(Control, Info) ;
    umfpack_dl_report_status(Control, status) ;

   //>>>convert  solution X to Field3D phi
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


        }//>>>end convert
    }



}


void EF_Solver3D_UMFPACK::solve_Exyz(Field* phi, Field* Ex, Field* Ey, Field* Ez)
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


void EF_Solver3D_UMFPACK::finishUMFPACK()
{
    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);
}

#endif // for SuperLU_type == serial