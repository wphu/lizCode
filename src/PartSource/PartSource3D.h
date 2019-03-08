/*
PartSource3D class
*/

#ifndef PARTSOURCE3D_H
#define PARTSOURCE3D_H

#include <vector>
#include "PartSource.h"

class PartSource3D : public PartSource
{

public:
    //! Constructor for PSI between two species
    PartSource3D(PicParams& params, SmileiMPI* smpi) : PartSource(params, smpi) {};
    virtual ~PartSource3D(){};

    // sputtered particle number
    int nPartEmit;

    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void emitLoad(PicParams&,std::vector<Species*>&,int, ElectroMagn*){};

private:

};


#endif
