#ifndef DIAGNOSTIC2D_H
#define DIAGNOSTIC2D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"


class Diagnostic2D : public Diagnostic {

public :

    Diagnostic2D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields);
    virtual ~Diagnostic2D() {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime );

    std::vector<Field*> particleFlux;
    std::vector<Field*> heatFlux;
    std::vector<Field*> averageAngle;
    std::vector<Field*> psiRate;   //sputteringRate, depositionRate, and so on

    std::vector<Field*> particleFlux_global;
    std::vector<Field*> heatFlux_global;
    std::vector<Field*> aveAngle_global;
    std::vector<Field*> psiRate_global;

    // calculate velocity and temperature of each species
    // not implemented
	void calVT(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime){};

protected :


};

#endif
