#include "Diagnostic2D.h"
#include "Field2D.h"
#include "PSI2D.h"
#include "SmileiMPI_Cart2D.h"

Diagnostic2D::Diagnostic2D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields, vector<PSI*>& vecPSI) :
Diagnostic(params)
{
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
    for (unsigned int ispec=0 ; ispec<params.species_param.size() ; ispec++) 
    {
        particleFlux[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        heatFlux[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        averageAngle[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        
        particleFlux_global[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        heatFlux_global[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        averageAngle_global[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
    }

    psiRate.resize(vecPSI.size());
    psiRate_global.resize(vecPSI.size());
    for(unsigned int ipsi=0; ipsi<vecPSI.size(); ipsi++)
    {
        psiRate[ipsi] = new Field2D(dim_global, "psiRate");
        psiRate_global[ipsi] = new Field2D(dim_global, "psiRate_global");
    }
}


void Diagnostic2D::run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime )
{
    Species *s1;
	Particles *psi_particles;
	double v_square, v_magnitude, energy;
	double mass_ov_2;
	double wlt;			// weight * cell_length / time for calculating flux
    double angle;
	int iAngle;
	double flux_temp;
    double xpn, ypn;
    int i, j, ic, jc;

	// reset diagnostic parameters to zero
	if( ((itime - 1) % step_dump) == 0 ) 
    {
		for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
		{
            particleFlux[iSpec]->put_to(0.0);
            heatFlux[iSpec]->put_to(0.0);
            averageAngle[iSpec]->put_to(0.0);
		}
        for(int iPsi = 0; iPsi < vecPSI.size();  iPsi++)
        {
            psiRate[iPsi]->put_to(0.0);
        }
	}

    // calculate particle flux, heat flux, average angle, and psiRate (like: deposition rate and sputtering rate)
    if( (itime % step_dump) > (step_dump - step_ave) || (itime % step_dump) == 0 )
    {
        for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
        {
            s1 = vecSpecies[iSpec];
            psi_particles = &(s1->psi_particles);
            mass_ov_2 = 0.5 * s1->species_param.mass;
            Field2D* particleFlux_temp =  static_cast<Field2D*>(particleFlux[iSpec]);
            Field2D* heatFlux_temp =  static_cast<Field2D*>(heatFlux[iSpec]);
            Field2D* averageAngle_temp =  static_cast<Field2D*>(averageAngle[iSpec]);
            for(int iPart = 0; iPart < psi_particles->size(); iPart++)
            {
                    xpn = psi_particles->position(0, iPart) * dx_inv_;  // normalized distance to the first node
                    ic  = floor(xpn);                   // index of the central node

                    ypn = psi_particles->position(1, iPart) * dy_inv_;  // normalized distance to the first node
                    jc   = floor(ypn);                  // index of the central node

                    i = ic-i_domain_begin; // index of first point for projection in x
                    j = jc-j_domain_begin; // index of first point for projection in y

                    v_square = psi_particles->momentum(0,iPart) * psi_particles->momentum(0,iPart) + psi_particles->momentum(1,iPart) * psi_particles->momentum(1,iPart) + psi_particles->momentum(2,iPart) * psi_particles->momentum(2,iPart);
                    v_magnitude = sqrt(v_square);
                    angle = 0.0;

                    (*particleFlux_temp)(i,j) += 1.0;
                    (*heatFlux_temp)(i,j) += mass_ov_2 * v_square;
                    (*averageAngle_temp)(i,j) += angle;
            }
        }
    }

    // MPI gather diagnostic parameters to master
    if( (itime % step_dump) == 0 ) {
		for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
		{
            smpi->reduce_sum_field(particleFlux[iSpec], particleFlux_global[iSpec]);
            smpi->reduce_sum_field(heatFlux[iSpec], heatFlux_global[iSpec]);
            smpi->reduce_sum_field(averageAngle[iSpec], averageAngle_global[iSpec]);
		}
	}
    
}