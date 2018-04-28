/*
PSI class
*/

#ifndef PSI_H
#define PSI_H

#include <vector>
#include <string>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Species.h"
#include "Diagnostic.h"
#include "Grid.h"

class Diagnostic;

using namespace std;
class PSI
{

public:
    //! Constructor for PSI between two species
    PSI(PicParams& params, SmileiMPI* smpi)
    {
        const_e = params.const_e;
        const_pi = params.const_pi;
        count_of_particles_to_insert_s2.resize(params.n_space[0]);
        for(int i = 0; i < count_of_particles_to_insert_s2.size(); i++)
        {
            count_of_particles_to_insert_s2[i] = 0;
        }
        posOffset = 0.1;

    };
    virtual ~PSI(){};

    void setRelPsi(PSI* relevantPsi)
    {
        relPsi = relevantPsi;
    }

    //! Method called in the main smilei loop to apply PSI at each timestep
    virtual void performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime){};

    // emit particles
    void emit(PicParams&, vector<Species*>&){};

    // calculate angle between two 3D vectors, unit is degree
    double angle_2vectors(double v1[], double v2[])
    {
        double mag1 = sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
        double mag2 = sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
        double dot_product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        if(mag1 * mag2 == 0.0)
        {
            ERROR("angle_2vectors equal to zero!");
        }
        double angle = acos(abs(dot_product) / (mag1 * mag2));
        return 180.0 * angle / const_pi;
    };

    //! Identification number of the PSI object
    int n_PSI;

    // PSI position : only left and right for 1D case
    string psiPos;

    // emit kind, regular or fieldEmit for injection PSI
    string emitKind;

    // relevant PSI, emitting number of two species may be relevant
    // such as nPartEmit(A) = relCoff * nPartEmit(B)
    PSI *relPsi;
    string relSpecies;

    // position offset of injected or sputtered particles
    double posOffset;
    // the energy/temperature of the new particles
    double emitTemp;
    double weight_const;

    // if is_self_consistent is false, only calculate the psiRate (like sputtering rate or reflection rate)
    // not create new particles
    bool is_self_consistent;

    //! Group of the species numbers that are associated for PSI.
    //> actually, each species gourp only contains one species for PSI
    //> for PSI_Injection, only species1 is used;
    //> for sputtering and secondary electron emission, species1 is the incident particle.
    unsigned int species1, species2;

    Species   *s1, *s2;
    Particles *p1, *p2;

    Particles new_particles;
    vector<int> count_of_particles_to_insert_s1;
    vector<int> count_of_particles_to_insert_s2;
    double const_e;
    double const_pi;


private:

};

#endif
