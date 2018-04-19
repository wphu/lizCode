#ifndef DIAGNOSTIC2D_H
#define DIAGNOSTIC2D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"
#include "PSI2D.h"
#include "Grid2D.h"
#include "Particles.h"

class Field;
class PSI;

class Diagnostic2D : public Diagnostic {

public :

    Diagnostic2D(PicParams& params, SmileiMPI* smpi, Grid* grid, ElectroMagn* EMfields, vector<PSI*>& vecPSI);
    virtual ~Diagnostic2D() {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime );

    // find the segment which particle cross
    bool find_cross_segment(Grid2D* grid2D, Particles *particles, int iPart, int *iLine_cross, int *iSegment_cross);
    
    // determine if a particle crosses a segment
    bool is_cross(double start_point[], double end_point[], double pos_new[], double pos_old[]);

    int n_species;
    int n_line;
    int n_segment_total;
    vector<int> n_segments;

    vector< vector< vector<double> > > particleFlux;
    vector< vector< vector<double> > > heatFlux;
    vector< vector< vector<double> > > averageAngle;
    vector< vector< vector<double> > > psiRate;   //sputteringRate, depositionRate, and so on

    vector< vector< vector<double> > > particleFlux_global;
    vector< vector< vector<double> > > heatFlux_global;
    vector< vector< vector<double> > > averageAngle_global;
    vector< vector< vector<double> > > psiRate_global;
    

    // calculate velocity and temperature of each species
    // not implemented
	void calVT(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime){};

protected :
    double dx_inv_, dy_inv_;
    int i_domain_begin, j_domain_begin;


};

#endif
