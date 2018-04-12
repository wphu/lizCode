/*
 * SmileIO_Cart2D.h
 *
 *  Created on: 3 juil. 2013
 */
#ifndef SMILEIO_CART2D_H
#define SMILEIO_CART2D_H

#include <string>
#include <vector>

#include "SmileiIO.h"
#include "Diagnostic2D.h"
#include "Grid2D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart2D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart2D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart2D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag );
    //! Destructor for SmileiIO
    ~SmileiIO_Cart2D();

    virtual void write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    //! Build memory and file space for // HDF5 write/read
    void createFieldsPattern( PicParams& params, ElectroMagn* fields );
    // Create particles h5 file pattern
    void createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies );

    // Create particles h5 file pattern
    void createDiagsPattern( PicParams& params, Diagnostic2D* diag2D );

    void initVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies );
    // calculate velocity distribution function
    void calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies);

    // write grid to grid.h5 file
    virtual void writeGrid(Grid* grid);
private:



};

#endif /* SMILEIO_CART2D_H_ */
