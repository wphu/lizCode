#ifndef DIAGNOSTIC1D_H
#define DIAGNOSTIC1D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"
#include "Field1D.h"
#include "PSI1D.h"
#include "Grid.h"
#include "PSI1D.h"

class Collisions;

class Diagnostic1D : public Diagnostic {

public:

    Diagnostic1D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields, vector<Species*>& vecSpecies, vector<Collisions*> &vecCollisions, vector<PSI*>& vecPSI);
    virtual ~Diagnostic1D();

    //run the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime ) ;

    //calculate velocity and temperature of each species
	void calVT(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime);

    //calculate total energy(particles and electric field)
    void calTotalEnergy(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime);

    //parameters for calculate energy distribution of incident ions on the target
    int n_energy;
    double energy_max;

    Field1D *ptclNum1D;

    vector<int> particle_number;                        //particle_number[i_spec]
    vector<double> total_particle_energy;               //total_paritcle_energy[ispec]
    double total_electric_field_energy;

	vector<double> particle_flux_left;                  //particle_flux_left[i_spec]
    vector<double> particle_flux_right;                 //particle_flux_right[i_spec]
	vector<double> heat_flux_left;                      //heat_flux_left[i_spec]
    vector<double> heat_flux_right;                     //heat_flux_right[i_spec]

    vector< vector<double> > angle_distribution_left;   //angle_distribution_left[i_spec][i_angle]
    vector< vector<double> > angle_distribution_right;  //angle_distribution_right[i_spec][i_angle]
    vector< vector<double> > energy_distribution_left;  //energy_distribution_left[i_spec][i_energy]
    vector< vector<double> > energy_distribution_right; //energy_distribution_right[i_spec][i_energy]

    //psi_rate has different meaning for different psi process, like physical sputtering rate, reflection rate
    vector<double> psi_rate_left;                            //psi_rate_left[i_psi]
    vector<double> psi_rate_right;                           //psi_rate_right[i_psi]

    //radiative_energy_collision[iCollision][iBin]
    vector< vector<double> > radiative_energy_collision;


protected :

    //Inverse of the spatial step 1/dx
    double dx_inv_;
    //parameters to project macroscopic velocity and temperature
    int index_domain_begin;

};

#endif
