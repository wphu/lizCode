/*
 * SmileiIO_Cart1D.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO_Cart1D.h"

#include <sstream>

#include "PicParams.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"
#include "ElectroMagn.h"
#include "Species.h"


using namespace std;

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies)
: SmileiIO( params, smpi )
{
    initVDF(params, smpi, fields, vecSpecies);
    reloadP(params, smpi, vecSpecies);

    if( smpi->isMaster() )
    {
        //create data patterns
        createFieldsPattern(params, fields);
        //at present, not calculate VDF in the core, VDF can be calculated from particles in the "restore" directory
        //createPartsPattern(params, fields, vecSpecies);

    }

}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
}

//> create hdf5 data hierarchical structure: datespace, dateset and so on
void SmileiIO_Cart1D::createFieldsPattern( PicParams& params, ElectroMagn* fields )
{
    fieldsGroup.dims_global[2] = params.n_space_global[0] + 1;
    fieldsGroup.dims_global[1] = 1;
    fieldsGroup.dims_global[0] = 1;

    fieldsGroup.ndims_[0] = fieldsGroup.dims_global[0];
    fieldsGroup.ndims_[1] = fieldsGroup.dims_global[1];
    fieldsGroup.ndims_[2] = fieldsGroup.dims_global[2];


    fieldsGroup.offset[0] = 0;
    fieldsGroup.offset[1] = 0;
    fieldsGroup.offset[2] = 0;


    fieldsGroup.stride[0] = 1;
    fieldsGroup.stride[1] = 1;
    fieldsGroup.stride[2] = 1;


    fieldsGroup.block[0] = 1;
    fieldsGroup.block[1] = 1;
    fieldsGroup.block[2] = 1;


    //For attribute
    fieldsGroup.aDims = 3;

    createFieldsGroup(fields);

} //END createPattern

//Create particles h5 file pattern
void SmileiIO_Cart1D::createDiagsPattern( PicParams& params, Diagnostic* diag)
{
 

}


void SmileiIO_Cart1D::createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies )
{

    //For particles, size ofdims_global should be 5: dims_global[nx][ny][nz][nvelocity][ntime]
    //But to be simple, the size is set 4, nz dimension is deleted.
    ptclsGroup.dims_global[2] = vx_dim;
    ptclsGroup.dims_global[1] = params.n_space_global[0];
    ptclsGroup.dims_global[0] = 1;

    ptclsGroup.ndims_[0] = ptclsGroup.dims_global[0];
    ptclsGroup.ndims_[1] = ptclsGroup.dims_global[1];
    ptclsGroup.ndims_[2] = ptclsGroup.dims_global[2];

    ptclsGroup.offset[0] = 0;
    ptclsGroup.offset[1] = 0;
    ptclsGroup.offset[2] = 0;

    ptclsGroup.stride[0] = 1;
    ptclsGroup.stride[1] = 1;
    ptclsGroup.stride[2] = 1;

    ptclsGroup.block[0] = 1;
    ptclsGroup.block[1] = 1;
    ptclsGroup.block[2] = 1;

    ptclsGroup.aDims = 3;

    createPartsGroup(vecSpecies);


}



void SmileiIO_Cart1D::initVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
{
    vx_dim = 200;

    vector<unsigned int> dims_VDF;
    vector<unsigned int> dims_VDF_global;
    dims_VDF.resize(4);
    dims_VDF_global.resize(4);

    dims_VDF[0] = params.n_space[0];
    dims_VDF[1] = 1;
    dims_VDF[2] = 1;
    dims_VDF[3] = vx_dim;

    dims_VDF_global[0] = params.n_space_global[0];
    dims_VDF_global[1] = 1;
    dims_VDF_global[2] = 1;
    dims_VDF_global[3] = vx_dim;

    for(int isp=0; isp<vecSpecies.size(); isp++)
    {
        vx_VDF.push_back(new Array4D());
        vx_VDF[isp]->allocateDims(dims_VDF);

        vx_VDF_global.push_back(new Array4D());
        vx_VDF_global[isp]->allocateDims(dims_VDF_global);

        vx_VDF_tot_global.push_back(new Array4D());
        vx_VDF_tot_global[isp]->allocateDims(dims_VDF_global);
    }

}



void SmileiIO_Cart1D::calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
{
    Species *s;
    Particles *p;


    for(int isp=0; isp<vx_VDF.size(); isp++)
    {
        s = vecSpecies[isp];
        p = &(s->particles);

        vxMax = 3*sqrt(2.0 * s->species_param.thermT[0] * params.const_e / s->species_param.mass);
        //WARNING("thermalVelocity" <<  s->species_param.thermT[0] );
        vxMin = -vxMax;
        vx_d = (vxMax - vxMin) / vx_dim;
        int vx_dim2 = vx_dim / 2;

        vx_VDF[isp]->put_to(0.0);
        for(int ibin = 0; ibin < ( s->bmin.size() ); ibin++)
        {
            for(int iPart = s->bmin[ibin]; iPart < s->bmax[ibin]; iPart++)
            {
                int ivx = p->momentum(0,iPart) / vx_d + vx_dim2;
                if( ivx < 0 ) {
                    ivx = 0;
                }
                if( ivx >= vx_dim ) {
                    ivx = vx_dim - 1;
                }
                (*vx_VDF[isp])(ibin,0,0,ivx) += 1.0;
            }
        }
        smpi->gatherVDF(vx_VDF_global[isp], vx_VDF[isp]);

        if( (itime % (params.dump_step + 1)) == 0 )
        {
            vx_VDF_tot_global[isp]->put_to(0.0);
        }

        for (int ibin = 0; ibin < vx_VDF_global[isp]->dims_[0]; ibin++)
        {
            for (int ivx = 0; ivx < vx_VDF_global[isp]->dims_[3]; ivx++)
            {
                (*vx_VDF_tot_global[isp])(0,0,0,ivx) += (*vx_VDF_global[isp])(ibin,0,0,ivx);
            }

        }

        if( (itime % (params.dump_step)) == 0 )
        {
            for (int ivx = 0; ivx < vx_VDF_global[isp]->dims_[3]; ivx++)
            {
                (*vx_VDF_tot_global[isp])(0,0,0,ivx) /= params.dump_step;
            }
        }

    }
}



//write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO_Cart1D::write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime)
{
    const char* h5_name;
    int n_dims_data = 3;
    int iDiag;

    hid_t       group_id;
    hid_t       dataspace_id;
    hid_t       memspace_id;
    hid_t       dataset_id;
    herr_t      status;

    hsize_t     dims_global[3];
    hsize_t     count[3];             //size of subset in the file
    hsize_t     offset[3];            //subset offset in the file
    hsize_t     stride[3];
    hsize_t     block[3];

    //set stride and block as default
    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;

    block[0] = 1;
    block[1] = 1;
    block[2] = 1;


    Diagnostic1D* diag1D = static_cast<Diagnostic1D*>(diag);
    if(params.is_calVDF)
    {
        calVDF( params, smpi, fields, vecSpecies, itime);
    }

    if(itime % params.dump_step == 0 && smpi->isMaster())
    {
        string step_output_string;

        step_output = itime / params.dump_step;
        step_output_string = to_string(step_output);
        
        if(step_output_max >= 10 && step_output_max <100)
        {
            if(step_output < 10)
            {
                step_output_string = "0" + step_output_string;
            }
        }
        else if(step_output_max >= 100 && step_output_max <1000)
        {
            if(step_output < 10)
            {
                step_output_string = "00" + step_output_string;
            }
            else if(step_output<100)
            {
                step_output_string = "0" + step_output_string;
            }
        }
        else if(step_output_max >= 1000)
        {
            WARNING("step_output_max is too large, please change the code in SmileiIO_Cart1D.cpp");
        }

        data_file_name = "data/data" + step_output_string + ".h5";
        data_file_id = H5Fcreate( data_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        //============= write attributes, n_dim: number of dimension ======================
        hid_t attrs_dataspace_id, attrs_id;
        int n_dim = 1;
        hsize_t attrs_dims[1];
        attrs_dims[0] = 1;
        attrs_dataspace_id = H5Screate_simple(1, attrs_dims, NULL);
        attrs_id           = H5Acreate2(data_file_id, "n_dim", H5T_STD_I32BE, attrs_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attrs_id, H5T_NATIVE_INT, &n_dim);
        H5Sclose(attrs_dataspace_id);
        H5Aclose(attrs_id);

        //=============write fields============================================
        fieldsGroup.group_id = H5Gcreate(data_file_id, "/Fields", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        for(int i = 0; i < fieldsGroup.dataset_stringName.size(); i++)
        {
            fieldsGroup.dataspace_id = H5Screate_simple(n_dims_data, fieldsGroup.dims_global, NULL);
            h5_name = fieldsGroup.dataset_stringName[i].c_str();
            fieldsGroup.dataset_id[i] = H5Dcreate2(fieldsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, fieldsGroup.dataspace_id,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            fieldsGroup.status = H5Dwrite(fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fieldsGroup.dataset_data[i]);
            fieldsGroup.status = H5Sclose(fieldsGroup.dataspace_id);
            fieldsGroup.status = H5Dclose(fieldsGroup.dataset_id[i]);
        }
        fieldsGroup.status = H5Gclose(fieldsGroup.group_id);


        //============== write Diagnostic: particle_flux, heat_flux and angle_distribution ============
        //============== particle_number, particle_energy, radiative_energy_collision      ============
        group_id = H5Gcreate(data_file_id, "/Diagnostic", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

        //total_electric_field_energy
        dims_global[0] = 1;
        dims_global[1] = 1;
        dims_global[2] = 1;

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);
 
        dataset_id   = H5Dcreate2(group_id, "total_electric_field_energy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->total_electric_field_energy));
        status       = H5Dclose(dataset_id);
        status       = H5Sclose(dataspace_id);

        //particle_number, total_particle_energy, particle_flux_left, particle_flux_right, heat_flux_left, heat_flux_right
        dims_global[0] = 1;
        dims_global[1] = 1;
        dims_global[2] = diag1D->n_species;

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);
 
        dataset_id   = H5Dcreate2(group_id, "particle_number", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->particle_number[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "total_particle_energy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->total_particle_energy[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "particle_flux_left", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->particle_flux_left[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "particle_flux_right", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->particle_flux_right[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "heat_flux_left", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->heat_flux_left[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "heat_flux_right", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->heat_flux_right[0]));
        status       = H5Dclose(dataset_id);

        status       = H5Sclose(dataspace_id);    

        //angle_distribution
        dims_global[0] = 1;
        dims_global[1] = diag1D->n_species;
        dims_global[2] = 90;

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);

        dataset_id = H5Dcreate2(group_id, "angle_distribution_left", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(int ispec = 0; ispec < diag1D->n_species; ispec++)
        {
            offset[0] = 0;
            offset[1] = ispec;
            offset[2] = 0;
            count[0]  = 1;
            count[1]  = 1;
            count[2]  = 90;

            memspace_id = H5Screate_simple(n_dims_data, count, NULL);
            status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
            status      = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, &(diag1D->angle_distribution_left[ispec][0])) ;
            status      = H5Sclose(memspace_id);
        }
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate2(group_id, "angle_distribution_right", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(int ispec = 0; ispec < diag1D->n_species; ispec++)
        {
            offset[0] = 0;
            offset[1] = ispec;
            offset[2] = 0;
            count[0]  = 1;
            count[1]  = 1;
            count[2]  = 90;

            memspace_id = H5Screate_simple(n_dims_data, count, NULL);
            status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
            status      = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, &(diag1D->angle_distribution_right[ispec][0]));
            status      = H5Sclose(memspace_id);
        }
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        //energy_distribution
        dims_global[0] = 1;
        dims_global[1] = diag1D->n_species;
        dims_global[2] = diag1D->n_energy;

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);

        dataset_id = H5Dcreate2(group_id, "energy_distribution_left", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(int ispec = 0; ispec < diag1D->n_species; ispec++)
        {
            offset[0] = 0;
            offset[1] = ispec;
            offset[2] = 0;
            count[0]  = 1;
            count[1]  = 1;
            count[2]  = diag1D->n_energy;

            memspace_id = H5Screate_simple(n_dims_data, count, NULL);
            status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
            status      = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, &(diag1D->energy_distribution_left[ispec][0]));
            status      = H5Sclose(memspace_id);
        }
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate2(group_id, "energy_distribution_right", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(int ispec = 0; ispec < diag1D->n_species; ispec++)
        {
            offset[0] = 0;
            offset[1] = ispec;
            offset[2] = 0;
            count[0]  = 1;
            count[1]  = 1;
            count[2]  = diag1D->n_energy;

            memspace_id = H5Screate_simple(n_dims_data, count, NULL);
            status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
            status      = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, &(diag1D->energy_distribution_right[ispec][0]));
            status      = H5Sclose(memspace_id);
        }
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        //psi_rate
        dims_global[0] = 1;
        dims_global[1] = 1;
        dims_global[2] = diag1D->n_psi;

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);
 
        dataset_id   = H5Dcreate2(group_id, "psi_rate_left", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->psi_rate_left[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "psi_rate_right", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->psi_rate_right[0]));
        status       = H5Dclose(dataset_id);
        status       = H5Sclose(dataspace_id);

        //radiative_energy_collision
        dims_global[0] = 1;
        dims_global[1] = diag1D->n_collision;
        dims_global[2] = diag1D->n_space_global[0];

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);

        dataset_id = H5Dcreate2(group_id, "radiative_energy_collision", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(int i_collision = 0; i_collision < diag1D->n_collision; i_collision++)
        {
            offset[0] = 0;
            offset[1] = i_collision;
            offset[2] = 0;
            count[0]  = 1;
            count[1]  = 1;
            count[2]  = diag1D->n_space_global[0];

            memspace_id = H5Screate_simple(n_dims_data, count, NULL);
            status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
            status      = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, &(diag1D->radiative_energy_collision[i_collision][0]));
            status      = H5Sclose(memspace_id);
        }
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        //sigma_left, sigma_right
        dims_global[0] = 1;
        dims_global[1] = 1;
        dims_global[2] = diag1D->n_species;

        dataspace_id = H5Screate_simple(n_dims_data, dims_global, NULL);

        dataset_id   = H5Dcreate2(group_id, "sigma_left", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->sigma_left[0]));
        status       = H5Dclose(dataset_id);

        dataset_id   = H5Dcreate2(group_id, "sigma_right", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status       = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(diag1D->sigma_right[0]));
        status       = H5Dclose(dataset_id);

        status       = H5Sclose(dataspace_id); 

        //close dignostic group
        status = H5Gclose(group_id);

        //close hdf5 file
        status = H5Fclose(data_file_id);


        /* to be removed to Diagnostic Class
        //==============write particle velocity distribution function=========================
        ptclsGroup.group_id = H5Gopen(data_file_id, "/VDF", H5P_DEFAULT);
        for(int i = 0; i < ptclsGroup.dataset_stringName.size(); i++)
        {
            h5_name = ptclsGroup.dataset_stringName[i].c_str();
            ptclsGroup.dataset_id[i] = H5Dopen( ptclsGroup.group_id, h5_name, H5P_DEFAULT);
        }

        ptclsGroup.count[0]  = ptclsGroup.dims_global[0];
        ptclsGroup.count[1]  = ptclsGroup.dims_global[1];
        ptclsGroup.count[2]  = ptclsGroup.dims_global[2];
        for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
        {
            ptclsGroup.memspace_id = H5Screate_simple (n_dims_data, ptclsGroup.count, NULL);
            ptclsGroup.dataspace_id = H5Dget_space (ptclsGroup.dataset_id[i]);
            ptclsGroup.status = H5Sselect_hyperslab (ptclsGroup.dataspace_id, H5S_SELECT_SET, ptclsGroup.offset,
                                                ptclsGroup.stride, ptclsGroup.count, ptclsGroup.block);
            ptclsGroup.status = H5Dwrite (ptclsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, ptclsGroup.memspace_id,
                                    ptclsGroup.dataspace_id, H5P_DEFAULT, ptclsGroup.dataset_data[i]);

            ptclsGroup.status = H5Sclose (ptclsGroup.memspace_id);
            ptclsGroup.status = H5Sclose (ptclsGroup.dataspace_id);
        }

        //close dataset, group, and so on
        for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
        {
            ptclsGroup.status = H5Dclose( ptclsGroup.dataset_id[i] );
        }
        ptclsGroup.status = H5Gclose( ptclsGroup.group_id );
        */

    }


}
