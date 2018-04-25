/* ==============================================================
Subroutine to evaluate empirical formulas for number and energy
Backscattering coefficients of light ions incident on Elemental
and Compound targets
Ref: Subroutines for some plasma surface interaction processes:
     hpysical sputtering, chemical erosion, radiation enhanced
     sublimation, backscattering and thermal evaporation.
================================================================*/
#ifndef PSI1D_BACKSCATTERING_H
#define PSI1D_BACKSCATTERING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "Backscattering.h"

class PSI1D_Backscattering : public PSI1D
{
public:
    PSI1D_Backscattering(
        PicParams& params,
        SmileiMPI* smpi,
        unsigned int psi_species1,
        unsigned int psi_species2,
        string psiPosition,
        double emitTemperature
    );

    ~PSI1D_Backscattering();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int, ElectroMagn*);

    // emit particles
    void emit(PicParams&, vector<Species*>&);

    Backscattering *backscattering;


private:


};


#endif
