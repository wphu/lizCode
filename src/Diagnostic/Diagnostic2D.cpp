#include "Diagnostic2D.h"
#include "Field2D.h"
#include "PSI2D.h"

Diagnostic2D::Diagnostic2D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields, vector<PSI*>& vecPSI) :
Diagnostic(params)
{
    // init particleFlux
    particleFlux.resize(params.species_param.size());
    averageAngle.resize(params.species_param.size());
    heatFlux.resize(params.species_param.size());
    particleFlux_global.resize(params.species_param.size());
    averageAngle_global.resize(params.species_param.size());
    heatFlux_global.resize(params.species_param.size());
    for (unsigned int ispec=0 ; ispec<params.species_param.size() ; ispec++) 
    {
        particleFlux[ispec] = new Field2D(dim_, params.species_param[ispec].species_type);
        heatFlux[ispec] = new Field2D(dim_, params.species_param[ispec].species_type);
        averageAngle[ispec] = new Field2D(dim_, params.species_param[ispec].species_type);
        
        particleFlux_global[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        heatFlux_global[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
        averageAngle_global[ispec] = new Field2D(dim_global, params.species_param[ispec].species_type);
    }

    psiRate.resize(vecPSI.size());
    psiRate_global.resize(vecPSI.size());
    for(unsigned int ipsi=0; ipsi<vecPSI.size(); ipsi++)
    {
        psiRate[ipsi] = new Field2D(dim_, "psiRate");
        psiRate_global[ipsi] = new Field2D(dim_global, "psiRate_global");
    }
}


void Diagnostic2D::run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime )
{
    
}