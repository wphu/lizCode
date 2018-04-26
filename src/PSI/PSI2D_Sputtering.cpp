#include "PSI2D_Sputtering.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "Diagnostic2D.h"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI2D_Sputtering::PSI2D_Sputtering(
    PicParams& params,
    SmileiMPI* smpi
    ):
PSI2D(params, smpi)
{



}

PSI2D_Sputtering::~PSI2D_Sputtering()
{

}

// initialize
void PSI2D_Sputtering::init(vector<Species*>& vecSpecies)
{
    Species   *s1, *s2;
    Particles *p1, *p2;


    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->psi_particles);
    p2 = &(s2->psi_particles);

    double an1 = s1->species_param.atomic_number;
    double am1 = s1->species_param.atomic_mass;
    double an2 = s2->species_param.atomic_number;
    double am2 = s2->species_param.atomic_mass;
    double es = s2->species_param.surface_binding_energy;
    double n = s2->species_param.density_solid;

    sputtering = new PhysicalSputtering_EmpiricalFormula(an1, am1, an2, am2, es, n);
}

// Calculates the PSI2D for a given Collisions object
void PSI2D_Sputtering::performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime)
{
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double energy_incident;
    double v_square, v_magnitude;
    double momentum[3];
    // sputtering probability
    double pSput;
    int iDim;
    bool has_find;
    int iLine_cross, iSegment_cross;
    Species   *s1, *s2;
    Particles *p1, *p2;

    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->psi_particles);
    p2 = &(s2->psi_particles);

    Diagnostic2D *diag2D = static_cast<Diagnostic2D*>(diag);
    Grid2D *grid2D = static_cast<Grid2D*>(grid);

    iDim = 0;
    int nPartEmit = 0;
    int nPart = p1->size();
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        has_find = diag2D->find_cross_segment(grid2D, p1, iPart, iLine_cross, iSegment_cross);
        
        momentum[0] = p1->momentum(0,iPart);
        momentum[1] = p1->momentum(1,iPart);
        momentum[2] = p1->momentum(2,iPart);
        v_square = pow(momentum[0], 2) + pow(momentum[1], 2) + pow(momentum[2], 2);
        energy_incident = 0.5 * s1->species_param.mass * v_square;
        theta = angle_2vectors(momentum, grid2D->lines[iLine_cross][iSegment_cross].normal);
        theta *= ( 180.0 / params.const_pi );
        
        pSput = sputtering->phy_sput_yield(theta, energy_incident / const_e);

        diag2D->psiRate[n_PSI][iLine_cross][iSegment_cross] += pSput;

        // add sputtered particle if pSput > ran_p
        double ran_p = (double)rand() / RAND_MAX;
        if( pSput > ran_p ) 
        {
            nPartEmit++;
        }


    };

    s2->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s2);
    new_particles.clear();
}
