#ifndef INTERPOLATOR3D0ORDER_H
#define INTERPOLATOR3D0ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D0Order : public Interpolator3D
{

public:
    Interpolator3D0Order(PicParams&, SmileiMPI*);
    ~Interpolator3D0Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);
    inline double compute( double* coeffx, double* coeffy, double* coeffz, Field3D* f, int idx, int idy, int idz) {
	double interp_res(0.);
	for (int iloc=0 ; iloc<2 ; iloc++) 
    {
	    for (int jloc=0 ; jloc<2 ; jloc++) 
        {
            for(int kloc = 0; kloc < 2; kloc++)
            {
                interp_res += coeffx[iloc] * coeffy[jloc] * coeffz[kloc] * (*f)(idx+iloc, idy+jloc, idz+kloc);
            }
	    }
	}
	return interp_res;
    };

private:
    // Last prim index computed
    int ip_, jp_, kp_;
    // Interpolation coefficient on Prim grid
    double coeffxp_[2], coeffyp_[2], coeffzp_[2];
    // Interpolation coefficient on Dual grid
    double coeffxd_[2], coeffyd_[2], coeffzd_[2];


};//END class

#endif
