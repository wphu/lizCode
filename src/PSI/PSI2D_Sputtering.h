
#ifndef PSI2D_SPUTTERING_H
#define PSI2D_SPUTTERING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI2D.h"
#include "PhysicalSputtering_EmpiricalFormula.h"

class PSI2D_Sputtering : public PSI2D
{

public:
    //! Constructor
    PSI2D_Sputtering(       
        PicParams& params,
        SmileiMPI* smpi,
        vector<Species*>& vecSpecies,
        unsigned int psi_species1,
        unsigned int psi_species2,
        string psiPosition,
        double emitTemperature);
    ~PSI2D_Sputtering();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams& params, SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime);

    void init(std::vector<Species*>&);
    
    Sputtering *sputtering;
private:


};


#endif
