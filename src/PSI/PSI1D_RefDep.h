/* ==============================================================
Subroutine to evaluate empirical formulas for number and energy
Backscattering coefficients of light ions incident on Elemental
and Compound targets
Ref: Subroutines for some plasma surface interaction processes:
     hpysical sputtering, chemical erosion, radiation enhanced
     sublimation, backscattering and thermal evaporation.
================================================================*/
#ifndef PSI1D_REFDEP_H
#define PSI1D_REFDEP_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "Backscatterin_EmpiricalFormula.h"

class PSI1D_RefDep : public PSI1D
{
public:
    PSI1D_RefDep(
        PicParams& params,
        SmileiMPI* smpi,
        vector<Species*>& vecSpecies,
        unsigned int psi_species1,
        unsigned int psi_species2,
        bool psi_is_self_consistent,
        string psiPosition,
        double emitTemperature
    );

    ~PSI1D_RefDep();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime);

    // emit particles
    void emit(PicParams&, vector<Species*>&);

    void init(std::vector<Species*>&);
    
    Backscattering *backscattering;

private:


};


#endif
