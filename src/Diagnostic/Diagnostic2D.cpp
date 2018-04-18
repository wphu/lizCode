#include "Diagnostic2D.h"
#include "Field2D.h"
#include "PSI2D.h"
#include "SmileiMPI_Cart2D.h"

#include <algorithm>

Diagnostic2D::Diagnostic2D(PicParams& params, Grid* grid, SmileiMPI* smpi, ElectroMagn* EMfields, vector<PSI*>& vecPSI) :
Diagnostic(params)
{
    Grid2D* grid2D = static_cast<Grid2D*>(grid);

    dx_inv_   = 1.0/params.cell_length[0];
    dy_inv_   = 1.0/params.cell_length[1];

    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
    j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);

    // init particleFlux
    particleFlux.resize(params.species_param.size());
    heatFlux.resize(params.species_param.size());
    averageAngle.resize(params.species_param.size());
    particleFlux_global.resize(params.species_param.size());
    heatFlux_global.resize(params.species_param.size());
    averageAngle_global.resize(params.species_param.size());
    for (unsigned int iSpec=0 ; iSpec<params.species_param.size() ; iSpec++) 
    {
        particleFlux[iSpec].resize(grid2D->lines.size());
        heatFlux[iSpec].resize(grid2D->lines.size());
        averageAngle[iSpec].resize(grid2D->lines.size());
        
        particleFlux_global[iSpec].resize(grid2D->lines.size());
        heatFlux_global[iSpec].resize(grid2D->lines.size());
        averageAngle_global[iSpec].resize(grid2D->lines.size());
        for(int iLine = 0; iLine < grid2D->lines.size(); iLine++)
        {
            particleFlux[iSpec][iLine].resize(grid2D->lines[iLine].size());
            heatFlux[iSpec][iLine].resize(grid2D->lines[iLine].size());
            averageAngle[iSpec][iLine].resize(grid2D->lines[iLine].size());
            
            particleFlux_global[iSpec][iLine].resize(grid2D->lines[iLine].size());
            heatFlux_global[iSpec][iLine].resize(grid2D->lines[iLine].size());
            averageAngle_global[iSpec][iLine].resize(grid2D->lines[iLine].size());
        }
    }

    /*
    psiRate.resize(vecPSI.size());
    psiRate_global.resize(vecPSI.size());
    for(unsigned int ipsi=0; ipsi<vecPSI.size(); ipsi++)
    {
        psiRate[ipsi] = new Field2D(dim_global, "psiRate");
        psiRate_global[ipsi] = new Field2D(dim_global, "psiRate_global");
    }
    */
}


void Diagnostic2D::run( SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime )
{
    Species *s1;
    Particles *p1;
	Particles *psi_particles;
    bool has_find;
    int iLine_cross, iSegment_cross;
	double v_square, v_magnitude, energy;
	double mass_ov_2;
	double wlt;			// weight * cell_length / time for calculating flux
    double angle;
	int iAngle;
	double flux_temp;
    
    Grid2D* grid2D = static_cast<Grid2D*>(grid);

	// reset diagnostic parameters to zero
	if( ((itime - 1) % step_dump) == 0 ) 
    {
		for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
		{
            for(int iLine = 0; iLine < particleFlux[iSpec].size(); iLine++)
            {
                for(int iSegment = 0; iSegment < particleFlux[iSpec][iLine].size(); iSegment++)
                {
                    particleFlux[iSpec][iLine][iSegment] = 0.0;
                    heatFlux[iSpec][iLine][iSegment] = 0.0;
                    averageAngle[iSpec][iLine][iSegment] = 0.0;
                }

            }

		}
        /*
        for(int iPsi = 0; iPsi < vecPSI.size();  iPsi++)
        {
            psiRate[iPsi]->put_to(0.0);
        }
        */
	}


    // absorb particles which hit wall, and calcualte particle flux, heat flux, and average angles
    for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
    {
        s1 = vecSpecies[iSpec];
        p1 = &(s1->particles);
        s1->indexes_of_particles_to_absorb.clear();
        mass_ov_2 = 0.5 * s1->species_param.mass;
        for (int ibin = 0 ; ibin < (unsigned int)s1->bmin.size() ; ibin++) 
        {
            for (int iPart=(unsigned int)s1->bmin[ibin] ; iPart<(unsigned int)s1->bmax[ibin]; iPart++ ) 
            {
                has_find = find_cross_segment(grid2D, p1, iPart, &iLine_cross, &iSegment_cross);
                if( has_find )
                {
                    s1->indexes_of_particles_to_absorb.push_back(iPart);
                    if( (itime % step_dump) > (step_dump - step_ave) || (itime % step_dump) == 0 )
                    {
                        v_square = pow(p1->momentum(0,iPart), 2) + pow(p1->momentum(1,iPart), 2) + pow(p1->momentum(2,iPart), 2);
                        v_magnitude = sqrt(v_square);
                        angle = 0.0;
                        particleFlux[iSpec][iLine_cross][iSegment_cross] += 1.0;
                        heatFlux[iSpec][iLine_cross][iSegment_cross] += mass_ov_2 * v_square;
                        averageAngle[iSpec][iLine_cross][iSegment_cross] += angle;
                    }

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
    if( (itime % step_dump) == 0 ) {
		for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
		{
            for(int iLine = 0; iLine < particleFlux[iSpec].size(); iLine++)
            {
                smpi->reduce_sum_double( &(particleFlux[iSpec][iLine][0]), &(particleFlux_global[iSpec][iLine][0]), particleFlux[iSpec][iLine].size() );
                smpi->reduce_sum_double( &(heatFlux[iSpec][iLine][0]), &(heatFlux_global[iSpec][iLine][0]), heatFlux[iSpec][iLine].size() );
                smpi->reduce_sum_double( &(averageAngle[iSpec][iLine][0]), &(averageAngle_global[iSpec][iLine][0]), averageAngle[iSpec][iLine].size() );
            }
		}
	}
    
}


bool Diagnostic2D::find_cross_segment(Grid2D *grid2D, Particles *particles, int iPart, int *iLine_cross, int *iSegment_cross)
{
    bool has_find = false;
    double xpn, ypn;
    int ic, jc, ic0, jc0, ic1, jc1;
    double pos_new[2];
    double pos_old[2];
    vector<int> vecSegment;

    //Locate particle on the primal grid & calculate the projection coefficients
    xpn = particles->position(0, iPart) * dx_inv_;  // normalized distance to the first node
    ic  = floor(xpn);                   // index of the central node

    ypn = particles->position(1, iPart) * dy_inv_;  // normalized distance to the first node
    jc   = floor(ypn);                  // index of the central node



    pos_new[0] = particles->position(0, iPart);
    pos_new[1] = particles->position(1, iPart);
    pos_old[0] = particles->position_old(0, iPart);
    pos_old[1] = particles->position_old(1, iPart);

    if( grid2D->iswall_global_2D[ic][jc] == 1 || grid2D->iswall_global_2D[ic+1][jc] == 1 || grid2D->iswall_global_2D[ic+1][jc+1] == 1
        ||grid2D->iswall_global_2D[ic][jc+1] == 1 ) 
    {
        ic0 = ic;
        jc0 = jc;
        ic1 = particles->position_old(0, iPart) * dx_inv_;
        jc1 = particles->position_old(0, iPart) * dy_inv_;
        if(ic0 > ic1)
        {
            int ic_temp = ic0;
            ic0 = ic1;
            ic1 = ic_temp;
        }
        if(jc0 > ic1)
        {
            int jc_temp = jc0;
            jc0 = jc1;
            jc1 = jc0;
        }
        
        for(int iLine = 0; iLine < grid2D->lines.size(); iLine++)
        {
            // find segments which the particle crosses
            for(int iSegment = 0; iSegment < grid2D->lines[iLine].size(); iSegment++)
            {
                if( grid2D->lines[iLine][iSegment].grid_point[0] >= ic0 && grid2D->lines[iLine][iSegment].grid_point[0] <= ic1
                    && grid2D->lines[iLine][iSegment].grid_point[1] >= jc0 && grid2D->lines[iLine][iSegment].grid_point[1] <= jc1 )
                {
                    vecSegment.push_back(iSegment);
                }
            }
            // determine if the segment cross the particle trajectory
            for(int i = 0; i < vecSegment.size(); i++)
            {
                if( is_cross( grid2D->lines[iLine][vecSegment[i]].start_point, grid2D->lines[iLine][vecSegment[i]].end_point, pos_new, pos_old )
                {
                    *iLine_cross = iLine;
                    *iSegment_cross = vecSegment[i];
                    has_find = true;
                    return has_find;
                }
            }
        }
    }
    return has_find;
}

// ref: https://www.cnblogs.com/wuwangchuxin0924/p/6218494.html
bool Diagnostic2D::is_cross(double start_point[], double end_point[], double pos_new[], double pos_old[])
{
    if(!( min(start_point[0],end_point[0]) <= max(pos_new[0],pos_old[0]) && min(pos_new[1],pos_old[1]) <= max(start_point[1],end_point[1])
       && min(pos_new[0], pos_old[0]) <= max(start_point[0],end_point[0]) && min(start_point[1],end_point[1]) <= max(pos_new[1],pow_old[1]) ))
    {
        return false;
    }
    double u, v, w, z;
    u = (pos_new[0] - start_point[0]) * (end_point[1] - start_point[1]) - (end_point[0] - start_point[0]) * (pos_new[1] - start_point[1]);
    v = (pos_old[0] - start_point[0]) * (end_point[1] - start_point[1]) - (end_point[0] - start_point[0]) * (pos_old[1] - start_point[1]);
    w = (start_point[0] - pos_new[0]) * (pos_old[1] - pos_new[1])       - (pos_old[0]   - pos_new[0])     * (start_point[1] - pos_new[1]);
    z = (end_point[0] - pos_new[0])   * (pos_old[1] - pos_new[1])       - (pos_old[0]   - pos_new[0])     * (end_point[1] - pos_new[1]);
    return( u * v <= 0.000000000001 && w * z <= 0.000000000001 );
}