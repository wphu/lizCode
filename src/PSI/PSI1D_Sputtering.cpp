#include "PSI1D_Sputtering.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "Diagnostic1D.h"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI1D_Sputtering::PSI1D_Sputtering(
    PicParams& params,
    SmileiMPI* smpi,
    vector<Species*>& vecSpecies,
    int n_psi_in,
    unsigned int psi_species1,
    unsigned int psi_species2,
    bool psi_is_self_consistent,
    string psiPosition,
    double emitTemperature
):
PSI1D(params, smpi)
{
    species1    = psi_species1;
    species2    = psi_species2;
    psiPos      = psiPosition;
    emitTemp    = emitTemperature;
    const_e     = params.const_e;
    n_psi       = n_psi_in;
    is_self_consistent = psi_is_self_consistent;
    
    count_of_particles_to_insert_s2.resize(params.n_space[0]);

    new_particles.initialize(0, params);

    init(vecSpecies);
}

PSI1D_Sputtering::~PSI1D_Sputtering()
{

}

// initialize
void PSI1D_Sputtering::init(vector<Species*>& vecSpecies)
{
    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];

    double an1 = s1->species_param.atomic_number;
    double am1 = s1->species_param.atomic_mass;
    double an2 = s2->species_param.atomic_number;
    double am2 = s2->species_param.atomic_mass;
    double es = s2->species_param.surface_binding_energy;
    double density = s2->species_param.density_solid;

    sputtering = new PhysicalSputtering_EmpiricalFormula(an1, am1, an2, am2, es, density);
}


// Calculates the PSI1D for a given Collisions object
void PSI1D_Sputtering::performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime)
{
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double ke;
    double v_square, v_magnitude;
    // sputtering probability
    double pSput;
    int iDim;

    p1 = &(s1->psi_particles);

    Diagnostic1D *diag1D = static_cast<Diagnostic1D*>(diag);

    for(int i = 0; i < count_of_particles_to_insert_s2.size(); i++)
    {
        count_of_particles_to_insert_s2[i] = 0;
    }

    iDim = 0;
    nPartEmit = 0;
    int nPart = p1->size();
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        if(p1->position(iDim,iPart) < smpi->getDomainLocalMin(iDim) && psiPos == "left") 
        {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            theta = acos(abs( p1->momentum(0,iPart) ) / sqrt( v_square ));
            theta *= ( 180.0 / params.const_pi );
            ke = 0.5 * s1->species_param.mass * v_square;
            pSput = sputtering->phy_sput_yield( theta, ke/const_e );
            double ran_p = (double)rand() / RAND_MAX;
            if( pSput > ran_p ) 
            {
                nPartEmit++;
            }
            diag1D->psi_rate_left[n_psi] += pSput;
        }
        else if(p1->position(iDim,iPart) > smpi->getDomainLocalMax(iDim) && psiPos == "right") 
        {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            theta = acos(abs( p1->momentum(0,iPart) ) / sqrt( v_square ));
            theta *= ( 180.0 / params.const_pi );
            ke = 0.5 * s1->species_param.mass * v_square;
            pSput = sputtering->phy_sput_yield( theta, ke/const_e );
            double ran_p = (double)rand() / RAND_MAX;
            if( pSput > ran_p ) 
            {
                nPartEmit++;
            }
            diag1D->psi_rate_right[n_psi] += pSput;
        }
    };

    if( smpi->isWestern() || smpi->isEastern() )
    {
        emit(params, vecSpecies);
        s2->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s2);
        new_particles.clear();
    }

}


void PSI1D_Sputtering::emit(PicParams& params, vector<Species*>& vecSpecies)
{
    new_particles.create_particles(nPartEmit);
    if(psiPos == "left")
    {
        count_of_particles_to_insert_s2.front() = nPartEmit;
        for(int iPart=0; iPart < nPartEmit; iPart++)
        {
            new_particles.position(0,iPart)=(((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
            new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

            double ran;
            do 
            {
                ran = (double)rand() / RAND_MAX;
            }
            while (ran == 0.0);
            // initialize using the Maxwell distribution function in x-direction
            double psm = sqrt(2.0 * const_e * emitTemp / s2->species_param.mass) * sqrt(-log(ran));
            double theta = M_PI*(double)rand() / RAND_MAX;
            double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
            new_particles.momentum(0,iPart) = abs( psm*sin(theta)*cos(phi) );
            new_particles.momentum(1,iPart) = 0.0;
            new_particles.momentum(2,iPart) = 0.0;

            new_particles.al_imp(0,iPart) = 0.0;
            new_particles.al_imp(1,iPart) = 0.0;
            new_particles.al_imp(2,iPart) = 0.0;
            new_particles.au_imp(0,iPart) = 0.0;
            new_particles.au_imp(1,iPart) = 0.0;
            new_particles.au_imp(2,iPart) = 0.0;
        }
    }
    else if(psiPos == "right")
    {
        
        count_of_particles_to_insert_s2.back() = nPartEmit;
        for(int iPart=0; iPart<nPartEmit; iPart++)
        {
           new_particles.position(0,iPart) = params.cell_length[0] * params.n_space_global[0] - (((double)rand() / RAND_MAX)) * params.cell_length[0] * posOffset;
           new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

           double ran;
           do 
           {
               ran = (double)rand() / RAND_MAX;
           }
           while (ran == 0.0);

           // initialize using the Maxwell distribution function in x-direction
           double psm = sqrt(2.0 * const_e * emitTemp / s2->species_param.mass) * sqrt(-log(ran));
           double theta = M_PI*(double)rand() / RAND_MAX;
           double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
           new_particles.momentum(0,iPart) = -abs( psm*sin(theta)*cos(phi) );
           new_particles.momentum(1,iPart) = 0.0;
           new_particles.momentum(2,iPart) = 0.0;

           new_particles.al_imp(0,iPart) = 0.0;
           new_particles.al_imp(1,iPart) = 0.0;
           new_particles.al_imp(2,iPart) = 0.0;
           new_particles.au_imp(0,iPart) = 0.0;
           new_particles.au_imp(1,iPart) = 0.0;
           new_particles.au_imp(2,iPart) = 0.0;
       }
    }
    else {
        ERROR("no such psiPos: " << psiPos);
    }
}
