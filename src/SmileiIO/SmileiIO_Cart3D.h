/*
 * SmileIO_Cart3D.h
 *
 *  Created on: 3 juil. 2013
 */
#ifndef SMILEIO_CART3D_H
#define SMILEIO_CART3D_H

#include <string>
#include <vector>

#include "SmileiIO.h"
#include "Diagnostic3D.h"
#include "Grid3D.h"
#include "Segment.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart3D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart3D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart3D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies);
    //! Destructor for SmileiIO
    ~SmileiIO_Cart3D();

    virtual void write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);
    virtual void read( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime){};

    //! Build memory and file space for // HDF5 write/read
    void createFieldsPattern( PicParams& params, ElectroMagn* fields );
    // Create particles h5 file pattern
    void createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies );

    // Create particles h5 file pattern
    void createDiagsPattern( PicParams& params, Diagnostic* diag);

    void initVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies );
    // calculate velocity distribution function
    void calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies);

    // write grid to grid.h5 file
    virtual void writeGrid(Grid* grid);

    // read grid from grid.h5 file
    virtual void readGrid(Grid* grid);
private:



};

#endif /* SMILEIO_CART3D_H_ */
