#ifndef DIAGNOSTICFACTORY_H
#define DIAGNOSTICFACTORY_H

#include "Diagnostic1D.h"
#include "Diagnostic2D.h"
#include "Diagnostic3D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Create appropriate IO environment for the geometry
//! \param params : Parameters
//! \param smpi : MPI environment
//  --------------------------------------------------------------------------------------------------------------------

class DiagnosticFactory {
public:

    static Diagnostic* create(PicParams& params, SmileiMPI* smpi, Grid* grid, ElectroMagn* EMfields, vector<Species*>& vecSpecies, vector<Collisions*> &vecCollisions, vector<PSI*>& vecPSI) {
        Diagnostic* diag = NULL;
        if ( params.geometry == "1d3v" ) {
            diag = new Diagnostic1D(params, smpi, EMfields, vecSpecies, vecCollisions, vecPSI);
        }
        else if ( params.geometry == "2d3v" ) {
            diag = new Diagnostic2D(params, smpi, grid, EMfields, vecSpecies, vecCollisions, vecPSI);
        }
        else if ( params.geometry == "3d3v" ) {
            diag = new Diagnostic3D(params, smpi, grid, EMfields, vecSpecies, vecCollisions, vecPSI);
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }

        return diag;
    } // END createGlobalDiagnostics


};

#endif
