
#ifndef PSI2D_SEE_H
#define PSI2D_SEE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI2D.h"

class PSI2D_SEE : public PSI2D
{

public:
    //! Constructor for Collisions between two species
    PSI2D_SEE(PicParams& params, SmileiMPI* smpi);
    ~PSI2D_SEE();



    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime);

private:


};


#endif
