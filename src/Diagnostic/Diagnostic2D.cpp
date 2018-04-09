#include "Diagnostic2D.h"

Diagnostic2D::Diagnostic2D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields, vector<PSI*>& vecPSI) :
Diagnostic(params)
{
    // init particleFlux
    particleFlux.resize(ispec<params.species_param.size());
    particleFlux_global.resize(ispec<params.species_param.size());
    for (unsigned int ispec=0 ; ispec<params.species_param.size() ; ispec++) 
    {
        particleFlux = new Field2D(dim_, params.species_param.species_type);
        heatFlux = new Field2D(dim_, params.species_param.species_type);
        averageAngle = new Field2D(dim_, params.species_param.species_type);
        
        particleFlux_global = new Field2D(dim_global, params.species_param.species_type);
        heatFlux_global = new Field2D(dim__global, params.species_param.species_type);
        averageAngle_global = new Field2D(dim__global, params.species_param.species_type);
    }

    psiRate.resize(vecPSI.size());
    for(unsigned int ipsi=0; ipsi<vecPSI.size(); ipsi++)
    {
        psiRate = new Field2D(dim_, "psiRate");
        psiRate_global = new Field2D(dim_global, "psiRate_global");
    }
}


void Diagnostic2D::run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime )
{
    
}