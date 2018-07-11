#include "Diagnostic3D.h"
#include "PSI3D.h"
#include "SmileiMPI_Cart3D.h"

#include <algorithm>

Diagnostic3D::Diagnostic3D(PicParams& params, SmileiMPI* smpi, Grid* grid, ElectroMagn* EMfields, vector<PSI*>& vecPSI) :
Diagnostic(params)
{
    dims_global.resize(3);
    for(int i = 0; i < 3; i++)
    {
        dims_global[i] = params.n_space_global[i] + 1;
    }

    Grid3D* grid3D = static_cast<Grid3D*>(grid);

    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    dx_inv_   = 1.0/params.cell_length[0];
    dy_inv_   = 1.0/params.cell_length[1];
    dz_inv_   = 1.0/params.cell_length[2];

    SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);
    i_domain_begin = smpi3D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi3D->getCellStartingGlobalIndex(1);
    k_domain_begin = smpi3D->getCellStartingGlobalIndex(2);

    n_species = params.species_param.size();
    particleFlux.resize(n_species);
    heatFlux.resize(n_species);
    particleFlux_global.resize(n_species);
    heatFlux_global.resize(n_species);

    for(int i_species = 0; i_species < n_species; i_species++)
    {
        particleFlux[i_species]         = new Field3D(dims_global, ("particleFlux"          + params.species_param[i_species].species_type).c_str());
        heatFlux[i_species]             = new Field3D(dims_global, ("heatFlux"              + params.species_param[i_species].species_type).c_str());
        particleFlux_global[i_species]  = new Field3D(dims_global, ("particleFlux_global"   + params.species_param[i_species].species_type).c_str());
        heatFlux_global[i_species]      = new Field3D(dims_global, ("heatFlux_global"       + params.species_param[i_species].species_type).c_str());
    }
}


void Diagnostic3D::run( SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime )
{
    Species *s1;
    Particles *p1;
	Particles *psi_particles;
    bool has_find;
    bool is_in_wall;
    int iLine_cross, iSegment_cross;
	double v_square, v_magnitude, energy;
	double mass_ov_2;
	double wlt0, wlt;			// weight * cell_length / time for calculating flux
    double angle;
	int iAngle;
	double flux_temp;
    double length1, length2, length12;
    unsigned int ic, jc, kc;
    
    Grid3D* grid3D = static_cast<Grid3D*>(grid);

	// reset diagnostic parameters to zero
	if( ((itime - 1) % step_dump) == 0 ) 
    {
		for(int i_species = 0; i_species < n_species; i_species++)
		{
            particleFlux[i_species]->put_to(0.0);
            heatFlux[i_species]->put_to(0.0);
		}
	}


    // absorb particles which hit wall, and calcualte particle flux, heat flux, and average angles
    for(int i_species = 0; i_species < n_species; i_species++)
    {
        s1 = vecSpecies[i_species];
        p1 = &(s1->particles);
        s1->indexes_of_particles_to_absorb.clear();
        mass_ov_2 = 0.5 * s1->species_param.mass;

        for (int ibin = 0 ; ibin < (unsigned int)s1->bmin.size() ; ibin++) 
        {
            for (int iPart=(unsigned int)s1->bmin[ibin] ; iPart<(unsigned int)s1->bmax[ibin]; iPart++ ) 
            {
                ic = p1->position(0, iPart) * dx_inv_;
                jc = p1->position(1, iPart) * dy_inv_;
                kc = p1->position(2, iPart) * dz_inv_;
                //cout<<"iPart"<<ic<<" "<<jc<<" "<<kc<<endl;
                //vector<unsigned int> grid_dims = grid3D->iswall_global_3D.get_dims();
                //cout<<"dims "<<grid_dims[0]<<" "<<grid_dims[1]<<" "<<grid_dims[2]<<endl;
                if(grid3D->iswall_global_3D(ic,jc,kc) == 1 && grid3D->iswall_global_3D(ic,jc,kc+1) == 1 && grid3D->iswall_global_3D(ic,jc+1,kc) == 1 && grid3D->iswall_global_3D(ic,jc+1,kc+1) == 1
                && grid3D->iswall_global_3D(ic+1,jc,kc) == 1 && grid3D->iswall_global_3D(ic+1,jc,kc+1) == 1 && grid3D->iswall_global_3D(ic+1,jc+1,kc) == 1 && grid3D->iswall_global_3D(ic+1,jc+1,kc+1) == 1)
                {   
                    //if(iLine_cross != 4) { cout<<"iLine_cross = "<<iLine_cross<<endl; }
                    //cout<<"has find "<<endl;                   
                    if( (itime % step_dump) > (step_dump - step_ave) || (itime % step_dump) == 0 )
                    {
                        v_square = pow(p1->momentum(0,iPart), 2) + pow(p1->momentum(1,iPart), 2) + pow(p1->momentum(2,iPart), 2);
                        v_magnitude = sqrt(v_square);
                        (*particleFlux[i_species])(ic,jc,kc) += 1.0;
                        (*heatFlux[i_species])    (ic,jc,kc) += mass_ov_2 * v_square;
                    }
                    s1->indexes_of_particles_to_absorb.push_back(iPart);
                    
                }
            }//iPart
        }// ibin
        

        // copy PSI particles to psi_particles, because after MPi particle exchanging
        // the PSI particles will be erased
        for(int iPart=0; iPart<s1->indexes_of_particles_to_absorb.size(); iPart++)
        {
            int iPart_psi = s1->indexes_of_particles_to_absorb[iPart];
            p1->cp_particle(iPart_psi, s1->psi_particles);
        }
        s1->erase_particles_from_bins(s1->indexes_of_particles_to_absorb);

    }


    // MPI gather diagnostic parameters to master
    if( (itime % step_dump) == 0 ) 
    {
		for(int i_species = 0; i_species < n_species; i_species++)
		{
            s1 = vecSpecies[i_species];
			wlt0 = s1->species_param.weight * dx * dy / (timestep * step_ave);
		}
	}
 
}