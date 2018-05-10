/* Three-Body (TB) recombination colllsion
   ref: Kinetic Modelling of the Plasma Recombination, Contrib. Plasma Phys. 56, No. 6-8, 698 – 704 (2016) / DOI 10.1002/ctpp.201611004
   note: Now the collision is only specified for hydrogen isotopes, if it is applied to other element, the corresponding 
         parameters should be modified properly.
*/
#ifndef COLLISIONS1D_RECOMBINATION_TB_H
#define COLLISIONS1D_RECOMBINATION_TB_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"

using namespace std;

class Particles;

class Collisions1D_Recombination_TB : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Recombination_TB(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
        unsigned int n_col,
        vector<unsigned int> sg1,  // electron
        vector<unsigned int> sg2,  // ion
        vector<unsigned int> sg3,  // atom
        string CS_fileName);
    ~Collisions1D_Recombination_TB();


    // ke1 and ke2 are the energys of the two electrons
    double cross_section(double ke1, double ke2);
    void calculate_scatter_velocity(double v_magnitude, double mass1, double mass2,
                                    vector<double>& momentum_unit, vector<double>& momentum_temp);


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    // get the maximum value of crossSection*velocity
    double maxCV(Particles* particles, double eMass);

private:

    double energy_ionization_threshold;

    // the max principle quantum number of the product atom to calculate recombination collision cross section
    int nmax;

};


#endif
