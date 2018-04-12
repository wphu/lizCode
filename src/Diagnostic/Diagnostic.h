#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "PicParams.h"
#include "SmileiMPI.h"
#include "PSI.h"

#include <iostream>
#include <vector>

using namespace std;

class Diagnostic {

public :

    Diagnostic(PicParams &params);
    virtual ~Diagnostic() {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int timestep ) {};

    const unsigned int n_species;

    //! nDim_field (from params)
    const unsigned int n_dim_field;

    //! n_space (from params) always 3D
    const std::vector<unsigned int> n_space;
    const std::vector<unsigned int> n_space_global;
    std::vector<unsigned int> dim;
    std::vector<unsigned int> dim_global;


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
