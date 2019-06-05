#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "PicParams.h"
#include "SmileiMPI.h"
#include "PSI.h"
#include "Grid.h"

#include <iostream>
#include <vector>

using namespace std;

class Diagnostic {

public :

    Diagnostic(PicParams &params, SmileiMPI* smpi, vector<Species*>& vecSpecies, vector<Collisions*> &vecCollisions, vector<PSI*>& vecPSI);
    virtual ~Diagnostic() {};

    //run the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int timestep ) {};

    const int n_species;
    const int n_collision;
    const int n_psi;
    const int n_dim_field;

    //n_space (from params) always 3D
    const std::vector<int> n_space;
    const std::vector<int> n_space_global;
    std::vector<int> dim;
    std::vector<int> dim_global;


protected :

    // pi * 0.5
    double pi_ov_2;
    int step_dump;
    int step_ave;
    double timestep;
    double const_e;

    // Phi is the angle between the magnetic field and the y-direction
    double sinPhi, cosPhi;

    vector<double> sim_length;

    vector<unsigned int> oversize;
};

#endif
