#include "PartSource3D_Load.h"
#include "SmileiMPI_Cart3D.h"
#include "Field3D.h"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
PartSource3D_Load::PartSource3D_Load(
    PicParams& params,
    SmileiMPI* smpi,
    unsigned int load_species1,
    vector<double> mean_vel,
    string load_kind,
    int    load_number,
    int    every_time,
    double load_dn,
    double load_density,
    double load_temperature,
    double load_Pos_start,
    double load_Pos_end,
    double load_Pos_Ystart,
    double load_Pos_Yend,
    double load_Pos_Zstart,
    double load_Pos_Zend,
    int    load_step_update
):
PartSource3D (params, smpi)
{
    dt_ov_dx = params.timestep / params.cell_length[0];
    dt = params.timestep;

    species1        = load_species1;
    mean_velocity   = mean_vel;
    loadKind        = load_kind;
    loadNumber      = load_number;
    everyTime       = every_time;
    loadDn          = load_dn;
    loadDensity     = load_density;
    loadTemperature = load_temperature;
    loadPos_start   = load_Pos_start;
    loadPos_end     = load_Pos_end;
    loadPos_Ystart  = load_Pos_Ystart;
    loadPos_Yend    = load_Pos_Yend;
    loadPos_Zstart  = load_Pos_Zstart;
    loadPos_Zend    = load_Pos_Zend;
    step_update     = load_step_update;

    YZArea = 1.0;

    // define the x-direction range of source region
    // the MPI domain is not in the source region
    if(smpi->getDomainLocalMax(0) <= loadPos_start || smpi->getDomainLocalMin(0) >= loadPos_end) {
        loadPos_start = 0.0;
        loadBin_start = 0;
        loadPos_end = 0.0;
        loadBin_end = 0;

        loadPos_Ystart = 0.0;
        loadBin_Ystart = 0;
        loadPos_Yend = 0.0;
        loadBin_Yend = 0;

        loadPos_Zstart = 0.0;
        loadBin_Zstart = 0;
        loadPos_Zend = 0.0;
        loadBin_Zend = 0;
    }
    // the left end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(0) <= loadPos_start && smpi->getDomainLocalMax(0) > loadPos_start
    && smpi->getDomainLocalMax(0) <= loadPos_end) {
        double temp = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        loadBin_start = temp + 0.001;
        if(loadBin_start < 0) { loadBin_start = 0; }
        loadPos_start = smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        loadPos_end = smpi->getDomainLocalMax(0);
        loadBin_end = params.n_space[0] - 1;
    }
    // the right end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(0) < loadPos_end && smpi->getDomainLocalMax(0) >= loadPos_end
    && smpi->getDomainLocalMin(0) >= loadPos_start) {
        loadPos_start = smpi->getDomainLocalMin(0);
        loadBin_start = 0;
        loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        if(loadBin_end > params.n_space[0] - 1) { loadBin_end = params.n_space[0] - 1; }
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_end*params.cell_length[0] )) );
        loadPos_end = smpi->getDomainLocalMin(0) + (loadBin_end+1)*params.cell_length[0];
    }
    // the whole MPI domain is in the source region
    else if(smpi->getDomainLocalMin(0) >= loadPos_start && smpi->getDomainLocalMax(0) <= loadPos_end) {
        loadPos_start = smpi->getDomainLocalMin(0);
        loadBin_start = 0;
        loadPos_end = smpi->getDomainLocalMax(0);
        loadBin_end = params.n_space[0] - 1;
    }
    // the whole MPI domain source region is in the bin
    else if(smpi->getDomainLocalMin(0) < loadPos_start && smpi->getDomainLocalMax(0) > loadPos_end) {
        //loadBin_start = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        //loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0] )) );
        double temp = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        loadBin_start = temp + 0.001;
        if(loadBin_start < 0) { loadBin_start = 0; }
        loadPos_start = smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0];
        loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        if(loadBin_end > params.n_space[0] - 1) { loadBin_end = params.n_space[0] - 1; }
        loadPos_end = smpi->getDomainLocalMin(0) + (loadBin_end+1)*params.cell_length[0];

    }

    // define the y-direction range of source region
    // the MPI domain is not in the source region
    if(smpi->getDomainLocalMax(1) <= loadPos_Ystart || smpi->getDomainLocalMin(1) >= loadPos_Yend) {
        loadPos_start = 0.0;
        loadBin_start = 0;
        loadPos_end = 0.0;
        loadBin_end = 0;

        loadPos_Ystart = 0.0;
        loadBin_Ystart = 0;
        loadPos_Yend = 0.0;
        loadBin_Yend = 0;

        loadPos_Zstart = 0.0;
        loadBin_Zstart = 0;
        loadPos_Zend = 0.0;
        loadBin_Zend = 0;
    }
    // the down end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(1) <= loadPos_Ystart && smpi->getDomainLocalMax(1) > loadPos_Ystart
    && smpi->getDomainLocalMax(1) <= loadPos_Yend) {
        double temp = (loadPos_Ystart - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        loadBin_Ystart = temp + 0.001;
        if(loadBin_Ystart < 0) { loadBin_Ystart = 0; }
        loadPos_Ystart = smpi->getDomainLocalMin(1) + loadBin_Ystart*params.cell_length[1];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        loadPos_Yend = smpi->getDomainLocalMax(1);
        loadBin_Yend = params.n_space[1] - 1;
    }
    // the up end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(1) < loadPos_Yend && smpi->getDomainLocalMax(1) >= loadPos_Yend
    && smpi->getDomainLocalMin(1) >= loadPos_Ystart) {
        loadPos_Ystart = smpi->getDomainLocalMin(1);
        loadBin_Ystart = 0;
        loadBin_Yend = (loadPos_Yend - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        if(loadBin_Yend > params.n_space[1] - 1) { loadBin_Yend = params.n_space[1] - 1; }
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_end*params.cell_length[0] )) );
        loadPos_Yend = smpi->getDomainLocalMin(1) + (loadBin_Yend+1)*params.cell_length[1];
    }
    // the whole MPI domain is in the source region
    else if(smpi->getDomainLocalMin(1) >= loadPos_Ystart && smpi->getDomainLocalMax(1) <= loadPos_Yend) {
        loadPos_Ystart = smpi->getDomainLocalMin(1);
        loadBin_Ystart = 0;
        loadPos_Yend = smpi->getDomainLocalMax(1);
        loadBin_Yend = params.n_space[1] - 1;
    }
    // the whole MPI domain source region is in the bin
    else if(smpi->getDomainLocalMin(1) < loadPos_Ystart && smpi->getDomainLocalMax(1) > loadPos_Yend) {
        //loadBin_start = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        //loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0] )) );
        double temp = (loadPos_Ystart - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        loadBin_Ystart = temp + 0.001;
        if(loadBin_Ystart < 0) { loadBin_Ystart = 0; }
        loadPos_Ystart = smpi->getDomainLocalMin(1) + loadBin_Ystart*params.cell_length[1];
        loadBin_Yend = (loadPos_Yend - smpi->getDomainLocalMin(1)) / params.cell_length[1];
        if(loadBin_Yend > params.n_space[1] - 1) { loadBin_Yend = params.n_space[1] - 1; }
        loadPos_Yend = smpi->getDomainLocalMin(1) + (loadBin_Yend+1)*params.cell_length[1];

    }

   // define the z-direction range of source region
    // the MPI domain is not in the source region
    if(smpi->getDomainLocalMax(2) <= loadPos_Zstart || smpi->getDomainLocalMin(2) >= loadPos_Zend) {
        loadPos_start = 0.0;
        loadBin_start = 0;
        loadPos_end = 0.0;
        loadBin_end = 0;

        loadPos_Ystart = 0.0;
        loadBin_Ystart = 0;
        loadPos_Yend = 0.0;
        loadBin_Yend = 0;

        loadPos_Zstart = 0.0;
        loadBin_Zstart = 0;
        loadPos_Zend = 0.0;
        loadBin_Zend = 0;
    }
    // the down end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(2) <= loadPos_Zstart && smpi->getDomainLocalMax(2) > loadPos_Zstart
    && smpi->getDomainLocalMax(2) <= loadPos_Zend) {
        double temp = (loadPos_Zstart - smpi->getDomainLocalMin(2)) / params.cell_length[2];
        loadBin_Zstart = temp + 0.001;
        if(loadBin_Zstart < 0) { loadBin_Zstart = 0; }
        loadPos_Zstart = smpi->getDomainLocalMin(2) + loadBin_Zstart*params.cell_length[2];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        loadPos_Zend = smpi->getDomainLocalMax(2);
        loadBin_Zend = params.n_space[2] - 1;
    }
    // the up end of source region is in the MPI domain
    else if(smpi->getDomainLocalMin(2) < loadPos_Zend && smpi->getDomainLocalMax(2) >= loadPos_Zend
    && smpi->getDomainLocalMin(2) >= loadPos_Zstart) {
        loadPos_Zstart = smpi->getDomainLocalMin(2);
        loadBin_Zstart = 0;
        loadBin_Zend = (loadPos_Zend - smpi->getDomainLocalMin(2)) / params.cell_length[2];
        if(loadBin_Zend > params.n_space[2] - 1) { loadBin_Zend = params.n_space[2] - 1; }
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_end*params.cell_length[0] )) );
        loadPos_Zend = smpi->getDomainLocalMin(2) + (loadBin_Zend+1)*params.cell_length[2];
    }
    // the whole MPI domain is in the source region
    else if(smpi->getDomainLocalMin(2) >= loadPos_Zstart && smpi->getDomainLocalMax(2) <= loadPos_Zend) {
        loadPos_Zstart = smpi->getDomainLocalMin(2);
        loadBin_Zstart = 0;
        loadPos_Zend = smpi->getDomainLocalMax(2);
        loadBin_Zend = params.n_space[2] - 1;
    }
    // the whole MPI domain source region is in the bin
    else if(smpi->getDomainLocalMin(2) < loadPos_Zstart && smpi->getDomainLocalMax(2) > loadPos_Zend) {
        //loadBin_start = (loadPos_start - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_start] /=
        //( params.cell_length[0] / (( smpi->getDomainLocalMin(0) + (loadBin_start+1)*params.cell_length[0] ) - loadPos_start) );
        //loadBin_end = (loadPos_end - smpi->getDomainLocalMin(0)) / params.cell_length[0];
        //numPart_in_each_bin[loadBin_end] /=
        //( params.cell_length[0] / (loadPos_end - ( smpi->getDomainLocalMin(0) + loadBin_start*params.cell_length[0] )) );
        double temp = (loadPos_Zstart - smpi->getDomainLocalMin(2)) / params.cell_length[2];
        loadBin_Zstart = temp + 0.001;
        if(loadBin_Zstart < 0) { loadBin_Zstart = 0; }
        loadPos_Zstart = smpi->getDomainLocalMin(2) + loadBin_Zstart*params.cell_length[2];
        loadBin_Zend = (loadPos_Zend - smpi->getDomainLocalMin(2)) / params.cell_length[2];
        if(loadBin_Zend > params.n_space[2] - 1) { loadBin_Zend = params.n_space[2] - 1; }
        loadPos_Zend = smpi->getDomainLocalMin(2) + (loadBin_Zend+1)*params.cell_length[2];

    }



    // load particles by number_density and Temperature
    numPart_in_each_cell = loadDensity / params.species_param[species1].weight;
    if(loadKind == "nT") {
        for(int ibin = 0; ibin < numPart_in_each_bin.size(); ibin++)
        {
            numPart_in_each_bin[ibin] = (loadBin_Yend - loadBin_Ystart + 1) * (loadBin_Zend - loadBin_Zstart + 1) * numPart_in_each_cell;
            //cout<<"numPart_in_each_bin "<<numPart_in_each_bin[ibin]<<endl;
        }
    }
    else if(loadKind == "dn")
    {
        loadStep = 1.0 + loadNumber * params.species_param[species1].weight / (loadDn * params.timestep);
        loadRem = loadDn * loadStep * params.timestep / params.species_param[species1].weight - loadNumber;
        loadRemTot = 0.0;
        MESSAGE("loadStep = "<<loadStep);
    }

    DEBUG(0, "x loadBin_position: "<<loadBin_start<<" "<<loadBin_end);
    DEBUG(0, "y loadBin_position: "<<loadBin_Ystart<<" "<<loadBin_Yend);
    DEBUG(0, "z loadBin_position: "<<loadBin_Zstart<<" "<<loadBin_Zend);
}

PartSource3D_Load::~PartSource3D_Load()
{

}



// Calculates the PartSource for a given PartSource object
void PartSource3D_Load::emitLoad(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime, ElectroMagn* fields)
{
    Species   *s1;
    Particles *p1;
    int iPart, nPart;
    vector <double> cell_length;            // cell_length for initialize position, maybe not equal to real cell lenghth
    vector<double> max_jutt_cumul;
    double *indexes=new double[params.nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];

    if(everyTime == 0 && itime > 1) { return; }
    if(itime % step_update != 0){ return; }
    if(loadKind == "nT" && loadBin_end != loadBin_start && loadBin_Yend != loadBin_Ystart && loadBin_Zend != loadBin_Zstart)
    {
        s1 = vecSpecies[species1];
        p1 = &(s1->particles);

        cell_length.resize(params.nDim_particle);
        cell_length[0] = params.cell_length[0];
        cell_length[1] = params.cell_length[1];
        cell_length[2] = params.cell_length[2];
        max_jutt_cumul.resize(0);
        temp[0] = loadTemperature;
        temp[1] = loadTemperature;
        temp[2] = loadTemperature;
        vel[0] = mean_velocity[0];
        vel[1] = mean_velocity[1];
        vel[2] = mean_velocity[2];

        indexes_of_particles_to_erase.clear();
        for(int ibin = 0; ibin < count_of_particles_to_insert.size(); ibin++ )
        {
            count_of_particles_to_insert[ibin] = 0;
        }

        // erase unnecessary paritcles in source region
        for(int ibin=loadBin_start; ibin<=loadBin_end; ibin++)
        {
            for(int iPart = s1->bmin[ibin]; iPart < s1->bmax[ibin]; iPart++)
            {
                if( p1->position(1,iPart) > loadPos_Ystart && p1->position(1,iPart) < loadPos_Yend 
                 && p1->position(2,iPart) > loadPos_Zstart && p1->position(2,iPart) < loadPos_Zend ) 
                {
                    indexes_of_particles_to_erase.push_back(iPart);
                }
            }
            count_of_particles_to_insert[ibin] = numPart_in_each_bin[ibin];
            //cout<<count_of_particles_to_insert[ibin]<<endl;
        }
        s1->erase_particles_from_bins(indexes_of_particles_to_erase);

        // insert lost particles into bins
        new_particles.clear();
        for(int ibin = loadBin_start; ibin <= loadBin_end; ibin++ )
        {
            new_particles.create_particles(count_of_particles_to_insert[ibin]);
        }
        s1->insert_particles_to_bins(new_particles, count_of_particles_to_insert);


        // re-initialize paritcles in source region
        for(int ibin = loadBin_start; ibin <= loadBin_end; ibin++)
        {
            // Because bin is only in x-direction, y-direction has no bin
            // The "Bin" in loadBin_Ystart and loadBin_Yend is just to be consistent with the x-direction
            for(int j = loadBin_Ystart; j <= loadBin_Yend; j++)
            {
                for(int k = loadBin_Zstart; k <= loadBin_Zend; k++)
                {
                    //cout<<"i,j,k "<<ibin<<" "<<j<<" "<<k<<endl;
                    //cout<<"i,j,k "<<loadBin_end<<" "<<loadBin_Yend<<" "<<loadBin_Zend<<endl;
                    iPart = s1->bmax[ibin] - (loadBin_Yend - j) * (loadBin_Zend - loadBin_Zstart + 1) * numPart_in_each_cell - (loadBin_Zend - k + 1) * numPart_in_each_cell;
                    nPart = numPart_in_each_cell;
                    indexes[0] = smpi->getDomainLocalMin(0) + ibin*params.cell_length[0];
                    indexes[1] = smpi->getDomainLocalMin(1) + j   *params.cell_length[1];
                    indexes[2] = smpi->getDomainLocalMin(2) + k   *params.cell_length[2];

                    s1->initPosition(nPart, iPart, indexes, params.nDim_particle,
                                cell_length, s1->species_param.initPosition_type);

                    s1->initMomentum(nPart,iPart, temp, vel,
                                s1->species_param.initMomentum_type, max_jutt_cumul, params);

                    s1->initWeight_constant(nPart, species1, iPart, s1->species_param.weight);
                    s1->initCharge(nPart, species1, iPart, s1->species_param.charge);
                }
            }
        }
    }
    else if(loadKind == "dn" && itime%loadStep == 0 && loadBin_end != loadBin_start)
    {

    }

    delete [] indexes;
    delete [] temp;
    delete [] vel;
}
