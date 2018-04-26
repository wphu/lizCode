#include "PSI1D_RefDep.h"
#include "SmileiMPI.h"
#include "Field2D.h"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI1D_RefDep::PSI1D_RefDep(
    PicParams& params,
    SmileiMPI* smpi,
    vector<Species*>& vecSpecies,
    unsigned int psi_species1,
    unsigned int psi_species2,
    string psiPosition,
    double emitTemperature
):
PSI1D(params, smpi)
{
    species1 = psi_species1;
    species2 = psi_species2;
    psiPos = psiPosition;
    emitTemp = emitTemperature;
    const_e = params.const_e;
}

PSI1D_RefDep::~PSI1D_RefDep()
{

}

// initialize
void PSI1D_RefDep::init(vector<Species*>& vecSpecies)
{
    // init Backscattering_EmpiricalFormula
    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->psi_particles);
    p2 = &(s2->psi_particles);

    int nz1 = s1->species_param.atomic_number;
    int m1 = s1->species_param.atomic_mass;
    int ne = s2->species_param.ne;
    vector<int> nz2 = s2->species_param.nz2;
    vector<int> nw = s2->species_param.nw;

    backscattering = new Backscatterin_EmpiricalFormula(nz1, m1, ne, nz2, nw);
}

// Calculates the PSI1D for a given Collisions object
void PSI1D_RefDep::performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime)
{
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double ke;
    // Backscattering number and energy cofficients
    double rn, re;
    double v_square, v_magnitude;
    // sputtering probability
    double pSput;
    int iDim;

    iDim = 0;
    int nPartEmit = 0;
    int nPart = p1->size();
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        if( p1->position(iDim,iPart) < smpi->getDomainLocalMin(iDim) || p1->position(iDim,iPart) > smpi->getDomainLocalMax(iDim) ) {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            theta = abs( p1->momentum(0,iPart) ) / sqrt( v_square );
            theta *= ( 180.0 / params.const_pi );
            ke = 0.5 * s1->species_param.mass * v_square;
            //ke *= params.norm_temperature;
            backscattering->scatter( rn, re, theta, ke/const_e );
            double ran_p = (double)rand() / RAND_MAX;
            if( rn > ran_p ) {
                emitTemp = re * ke;
                new_particles.create_particle();
                if(psiPos == "left")
                {
                    new_particles.position(0,iPart)=(((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
                    new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

                    double ran;
                    do {
                        ran = (double)rand() / RAND_MAX;
                    }
                    while (ran == 0.0);

                    // initialize using the Maxwell distribution function in x-direction
                    double psm = sqrt(2.0 * emitTemp / s1->species_param.mass) * sqrt(-log(ran));
                    double theta = M_PI*(double)rand() / RAND_MAX;
                    double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
                    new_particles.momentum(0,iPart) = abs( psm*sin(theta)*cos(phi) );
                }
                else if(psiPos == "right")
                {
                    new_particles.position(0,iPart)=params.cell_length[0]*params.n_space_global[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
                    new_particles.position_old(0,iPart) = new_particles.position(0,iPart);

                    double ran;
                    do {
                        ran = (double)rand() / RAND_MAX;
                    }
                    while (ran == 0.0);
                    // initialize using the Maxwell distribution function in x-direction
                    double psm = sqrt(2.0 * emitTemp / s1->species_param.mass) * sqrt(-log(ran));
                    double theta = M_PI*(double)rand() / RAND_MAX;
                    double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
                    new_particles.momentum(0,iPart) = -abs( psm*sin(theta)*cos(phi) );
                }
                else
                {
                    return;
                }
                new_particles.momentum(1,nPartEmit) = 0.0;
                new_particles.momentum(2,nPartEmit) = 0.0;

                new_particles.al_imp(0,iPart) = 0.0;
                new_particles.al_imp(1,iPart) = 0.0;
                new_particles.al_imp(2,iPart) = 0.0;
                new_particles.au_imp(0,iPart) = 0.0;
                new_particles.au_imp(1,iPart) = 0.0;
                new_particles.au_imp(2,iPart) = 0.0;

                new_particles.weight(iPart) = weight_const;
                new_particles.charge(iPart) = s1->species_param.charge;
                nPartEmit++;
            }
        }
    };

    if( smpi->isWestern() || smpi->isEastern() ) {
        emit(params, vecSpecies);
        s2->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s2);
        new_particles.clear();
    };
    
}

void PSI1D_RefDep::emit(PicParams& params, vector<Species*>& vecSpecies)
{
    if(psiPos == "left")
    {
        count_of_particles_to_insert_s2.front() = nPartEmit;
    }
    else if(psiPos == "right")
    {
        count_of_particles_to_insert_s2.back() = nPartEmit;
    }
    else
    {
        ERROR("no such emitPos: " << psiPos);
    }
}


