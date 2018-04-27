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



private:



};


#endif
