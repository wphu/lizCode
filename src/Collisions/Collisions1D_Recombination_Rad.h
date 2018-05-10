// ref: Kinetic Modelling of the Plasma Recombination, Contrib. Plasma Phys. 56, No. 6-8, 698 – 704 (2016) / DOI 10.1002/ctpp.201611004

#ifndef COLLISIONS1D_RECOMBINATION_RAD_H
#define COLLISIONS1D_RECOMBINATION_RAD_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"

using namespace std;

class Particles;

class Collisions1D_Recombination_Rad : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Recombination_Rad(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
        unsigned int n_col,
        vector<unsigned int> sg1,  // electron
        vector<unsigned int> sg2,  // ion
        vector<unsigned int> sg3,  // atom
        string CS_fileName);

    ~Collisions1D_Recombination_Rad();


    double cross_section(double ke);
    void calculate_scatter_velocity(double v_magnitude, double mass1, double mass2,
                                    vector<double>& momentum_unit, vector<double>& momentum_temp);


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    // get the maximum value of crossSection*velocity
    double maxCV(Particles* particles, double eMass);

private:
    //>the ionization threshold energy
    double energy_ionization_threshold;

};


#endif
