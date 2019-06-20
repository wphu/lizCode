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
    int n_psi_in,
    unsigned int psi_species1,
    unsigned int psi_species2,
    unsigned int psi_species3,
    bool psi_is_self_consistent,
    string psiPosition,
    double emitTemperature
):
PSI1D(params, smpi)
{
    species1    = psi_species1;
    species2    = psi_species2;
    species3    = psi_species3;
    psiPos      = psiPosition;
    emitTemp    = emitTemperature;
    const_e     = params.const_e;
    n_psi       = n_psi_in;  
    is_self_consistent = psi_is_self_consistent;

    count_of_particles_to_insert_s3.resize(params.n_space[0]);

    new_particles.initialize(0, params);

    init(vecSpecies);
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
    s3 = vecSpecies[species3];

    int an1 = s1->species_param.atomic_number;
    int am1 = s1->species_param.atomic_mass;
    int ne2 = s2->species_param.ne2;
    vector<int> an2_vector = s2->species_param.an2_vector;
    vector<int> nw2_vector = s2->species_param.nw2_vector;

    backscattering = new Backscatterin_EmpiricalFormula(an1, am1, ne2, an2_vector, nw2_vector);
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
    int iDim;
    int nPart;
    int i_particle_new;

    Diagnostic1D *diag1D = static_cast<Diagnostic1D*>(diag);

    p1 = &(s1->psi_particles);

    for(int i = 0; i < count_of_particles_to_insert_s3.size(); i++)
    {
        count_of_particles_to_insert_s3[i] = 0;
    }

    iDim = 0;
    nPartEmit = 0;
    nPart = p1->size();
    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        if(p1->position(iDim,iPart) < smpi->getDomainLocalMin(iDim) && psiPos == "left") 
        {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            theta = acos(abs( p1->momentum(0,iPart) ) / sqrt( v_square ));
            theta *= ( 180.0 / params.const_pi );
            ke = 0.5 * s1->species_param.mass * v_square;
            ke /= const_e;
            backscattering->scatter( rn, re, theta, ke);

            if(rn > 1.0)
            {
                while(rn > 1.0)
                {
                    backscattering->scatter(rn, re, theta, ke);
                    cout<<"rn > 1:  "<<rn<<"  "<<re<<"  "<<theta<<"  "<<ke<<endl;
                }
            }

            double ran_p = (double)rand() / RAND_MAX;
            if( rn > ran_p ) 
            {
                if(rn > 1.0)
                {
                    cout<<"rn > 1:  "<<rn<<"  "<<theta<<"  "<<ke<<endl;
                }

                emitTemp = re * ke;
                new_particles.create_particle();
                i_particle_new = new_particles.size() - 1;

                new_particles.position(0,i_particle_new)     = (((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
                new_particles.position_old(0,i_particle_new) = p1->position(0,iPart);

                double ran;
                do 
                {
                    ran = (double)rand() / RAND_MAX;
                }
                while (ran == 0.0);

                // initialize using the Maxwell distribution function in x-direction
                double psm = sqrt(2.0 * const_e * emitTemp / s3->species_param.mass) * sqrt(-log(ran));
                double theta = M_PI*(double)rand() / RAND_MAX;
                double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
                new_particles.momentum(0,i_particle_new) = abs( psm*sin(theta)*cos(phi) );

                new_particles.momentum(1,i_particle_new) = 0.0;
                new_particles.momentum(2,i_particle_new) = 0.0;

                new_particles.al_imp(0,i_particle_new) = 0.0;
                new_particles.al_imp(1,i_particle_new) = 0.0;
                new_particles.al_imp(2,i_particle_new) = 0.0;
                new_particles.au_imp(0,i_particle_new) = 0.0;
                new_particles.au_imp(1,i_particle_new) = 0.0;
                new_particles.au_imp(2,i_particle_new) = 0.0;
                nPartEmit++;
            }   
            diag1D->psi_rate_left[n_psi] += rn;         
        }
        else if(p1->position(iDim,iPart) > smpi->getDomainLocalMax(iDim) && psiPos == "right") 
        {
            v_square = pow(p1->momentum(0,iPart),2) + pow(p1->momentum(1,iPart),2) + pow(p1->momentum(2,iPart),2);
            theta = acos(abs( p1->momentum(0,iPart) ) / sqrt( v_square ));
            theta *= ( 180.0 / params.const_pi );
            ke = 0.5 * s1->species_param.mass * v_square;
            ke /= const_e;
            backscattering->scatter(rn, re, theta, ke);
            if(rn > 1.0)
            {
                while(rn > 1.0)
                {
                    backscattering->scatter(rn, re, theta, ke);
                    cout<<"rn > 1:  "<<rn<<"  "<<re<<"  "<<theta<<"  "<<ke<<endl;
                }
            }
            
            double ran_p = (double)rand() / RAND_MAX;

            if( rn > ran_p ) 
            {
                emitTemp = re * ke;
                new_particles.create_particle();
                i_particle_new = new_particles.size() - 1;

                new_particles.position(0,i_particle_new) = params.cell_length[0]*params.n_space_global[0] - (((double)rand() / RAND_MAX))*params.cell_length[0]*posOffset;
                new_particles.position_old(0,i_particle_new) = p1->position(0,iPart);

                double ran;
                do 
                {
                    ran = (double)rand() / RAND_MAX;
                }
                while (ran == 0.0);
                // initialize using the Maxwell distribution function in x-direction
                double psm = sqrt(2.0 * const_e * emitTemp / s3->species_param.mass) * sqrt(-log(ran));
                double theta = M_PI*(double)rand() / RAND_MAX;
                double phi   = 2.0 * M_PI*(double)rand() / RAND_MAX;
                new_particles.momentum(0,i_particle_new) = -abs( psm*sin(theta)*cos(phi) );

                new_particles.momentum(1,i_particle_new) = 0.0;
                new_particles.momentum(2,i_particle_new) = 0.0;

                new_particles.al_imp(0,i_particle_new) = 0.0;
                new_particles.al_imp(1,i_particle_new) = 0.0;
                new_particles.al_imp(2,i_particle_new) = 0.0;
                new_particles.au_imp(0,i_particle_new) = 0.0;
                new_particles.au_imp(1,i_particle_new) = 0.0;
                new_particles.au_imp(2,i_particle_new) = 0.0;
                nPartEmit++;
            }
            diag1D->psi_rate_right[n_psi] += rn;
        }
    };

    
    if( smpi->isWestern() || smpi->isEastern() ) 
    {
        emit(params, vecSpecies);
        s3->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s3);
        new_particles.clear();
    };
    
}

void PSI1D_RefDep::emit(PicParams& params, vector<Species*>& vecSpecies)
{
    if(psiPos == "left")
    {
        count_of_particles_to_insert_s3.front() = nPartEmit;
    }
    else if(psiPos == "right")
    {
        count_of_particles_to_insert_s3.back() = nPartEmit;
    }
    else
    {
        ERROR("no such emitPos: " << psiPos);
    }
}


