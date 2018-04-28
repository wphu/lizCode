/*
PSI2D class
*/

#ifndef PSI2D_H
#define PSI2D_H

#include <vector>
#include "PSI.h"
#include "Diagnostic2D.h"
#include "Grid2D.h"

class PSI2D : public PSI
{

public:
    //! Constructor for PSI between two species
    PSI2D(PicParams& params, SmileiMPI* smpi) : PSI(params, smpi){};
    virtual ~PSI2D(){};




    //! Method called in the main smilei loop to apply PSI at each timestep
    virtual void performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime){};

    // calculate position coordinate of mirror reflection
    void cal_mirror_reflection(double start_point[], double end_point[], double position_old[], double position_new[])
    {
        double v0[2], v1[2], v2[3];
        v0[0] = start_point[0] - position_old[0];
        v0[1] = start_point[1] - position_old[1];
        v1[0] = 0.5 * (end_point[0] - start_point[0]);
        v1[1] = 0.5 * (end_point[1] - start_point[1]);
        v2[0] = 2.0 * (v0[0] + v1[0]);
        v2[1] = 2.0 * (v0[1] + v1[1]);
        position_new[0] = position_old[0] + v2[0];
        position_new[1] = position_old[1] + v2[1];
    }

    void cal_velocity(double normal[], double energy_s2, double momentum[])
    {
        double vt = sqrt(2.0 * energy_s2 * const_e / s2->species_param.mass);
        double x1;
        double x2;
        // velocity in the coordinate, whose y direction is parallel to the segment surface normal
        double momentum_normal[3];

        do {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        momentum_normal[0] = vt * sqrt( -log(x1) ) * sin(2.0 * const_pi * x2);

        do {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        // the velocity should point out the segment surface
        momentum_normal[1] = abs(vt * sqrt( -log(x1) ) * sin(2.0 * const_pi * x2));

        do {
            x1 = (double)rand() / RAND_MAX;
        }
        while (x1 == 0.0);
        x2 = (double)rand() / RAND_MAX;
        momentum_normal[2] = vt * sqrt( -log(x1) ) * sin(2.0 * const_pi * x2);

        momentum[0] = momentum_normal[0] * normal[0] + momentum_normal[1] * normal[1];
        momentum[1] = momentum_normal[0] * normal[1] + momentum_normal[1] * normal[0];
        momentum[2] = momentum_normal[2];
    }



private:



};


#endif
