#ifndef INTERPOLATOR3D_H
#define INTERPOLATOR3D_H


#include "Interpolator.h"
#include "SmileiMPI_Cart3D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 3D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D : public Interpolator
{
public:
    Interpolator3D(PicParams&params, SmileiMPI*smpi): Interpolator(params, smpi) {
        SmileiMPI_Cart3D* smpi3D = static_cast<SmileiMPI_Cart3D*>(smpi);
        i_domain_begin = smpi3D->getCellStartingGlobalIndex(0);
        j_domain_begin = smpi3D->getCellStartingGlobalIndex(1);
        k_domain_begin = smpi3D->getCellStartingGlobalIndex(2);
    };

    virtual ~Interpolator3D() {};
    virtual void mv_win(unsigned int shift) {i_domain_begin += shift;}
    virtual void setMvWinLimits(unsigned int shift) {i_domain_begin = shift;}

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    double dy_inv_;
    double dz_inv_;
    int i_domain_begin;
    int j_domain_begin;
    int k_domain_begin;
};

#endif
