#include "PSI2D_RefDep.h"
#include "SmileiMPI.h"
#include "Field2D.h"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PSI2D_RefDep::PSI2D_RefDep(
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
PSI2D(params, smpi)
{
    species1    = psi_species1;
    species2    = psi_species2;
    species3    = psi_species3;
    psiPos      = psiPosition;
    emitTemp    = emitTemperature;
    n_psi       = n_psi_in;
    is_self_consistent = psi_is_self_consistent;

    const_e = params.const_e;

    init(vecSpecies);
}

PSI2D_RefDep::~PSI2D_RefDep()
{

}

// initialize
void PSI2D_RefDep::init(vector<Species*>& vecSpecies)
{
    // init Backscattering_EmpiricalFormula
    s1 = vecSpecies[species1];
    s2 = vecSpecies[species2];
    p1 = &(s1->psi_particles);
    p2 = &(s2->psi_particles);

    int an1 = s1->species_param.atomic_number;
    int am1 = s1->species_param.atomic_mass;
    int ne2 = s2->species_param.ne2;
    vector<int> an2_vector = s2->species_param.an2_vector;
    vector<int> nw2_vector = s2->species_param.nw2_vector;

    backscattering = new Backscatterin_EmpiricalFormula(an1, am1, ne2, an2_vector, nw2_vector);
}

// Calculates the PSI2D for a given Collisions object
void PSI2D_RefDep::performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime)
{
    // rn: number coefficient, re: energy coefficient
    double rn, re;
    // the angle of particle velocity with the surface normal
    double theta;
    // kinetic energy_ion
    double energy_incident;
    double v_square, v_magnitude;
    double momentum[3];
    int iDim;
    bool has_find;
    bool is_in_wall;
    int iLine_cross, iSegment_cross;
    Species   *s1, *s2, *s3;
    Particles *p1, *p2, *p3;

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
        has_find = diag2D->find_cross_segment(grid2D, p1, iPart, iLine_cross, iSegment_cross, is_in_wall);
        
        momentum[0] = p1->momentum(0,iPart);
        momentum[1] = p1->momentum(1,iPart);
        momentum[2] = p1->momentum(2,iPart);
        v_square = pow(momentum[0], 2) + pow(momentum[1], 2) + pow(momentum[2], 2);
        energy_incident = 0.5 * s1->species_param.mass * v_square;
        theta = angle_2vectors(momentum, grid2D->lines[iLine_cross][iSegment_cross].normal);
        theta *= ( 180.0 / params.const_pi );
        
        backscattering->scatter(rn, re, theta, energy_incident / const_e);

        diag2D->psiRate[n_psi][iLine_cross][iSegment_cross] += rn;

        // add reflected particle if rn > ran_p
        double ran_p = (double)rand() / RAND_MAX;
        if( rn > ran_p ) 
        {
            nPartEmit++;
        }
    };

    s2->insert_particles_to_bins(new_particles, count_of_particles_to_insert_s2);
    new_particles.clear();

}
