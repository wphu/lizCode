#include "Interpolator3D0Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D0Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D0Order::Interpolator3D0Order(PicParams &params, SmileiMPI *smpi) : Interpolator3D(params, smpi)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    dz_inv_ = 1.0/params.cell_length[2];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D0Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{
    // Static cast of the electromagnetic fields
    Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
    Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
    Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
    Field3D* Bx3D = static_cast<Field3D*>(EMfields->Bx_);
    Field3D* By3D = static_cast<Field3D*>(EMfields->By_);
    Field3D* Bz3D = static_cast<Field3D*>(EMfields->Bz_);


    // Normalized particle position
    double xpn = particles.position(0, ipart)*dx_inv_;
    double ypn = particles.position(1, ipart)*dy_inv_;
    double zpn = particles.position(2, ipart)*dz_inv_;


    // Indexes of the central nodes
    ip_ = floor(xpn);
    jp_ = floor(ypn);
    kp_ = floor(zpn);


    // Declaration and calculation of the coefficient for interpolation
    double delta;

    delta   = xpn - (double)ip_;
    int iloc;
    if(delta < 0.5)
    {
        iloc = 0;
    }
    else
    {
        iloc = 1;
    }

    delta   = ypn - (double)jp_;
    int jloc;
    if(delta < 0.5)
    {
        jloc = 0;
    }
    else
    {
        jloc = 1;
    }

    delta   = zpn - (double)kp_;
    int kloc;
    if(delta < 0.5)
    {
        kloc = 0;
    }
    else
    {
        kloc = 1;
    }

    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    kp_ = kp_ - k_domain_begin;



    (*ELoc).x =  (*Ex3D)(ip_+iloc, jp_+jloc, kp_+kloc);
    (*ELoc).y =  (*Ey3D)(ip_+iloc, jp_+jloc, kp_+kloc);
    (*ELoc).z =  (*Ez3D)(ip_+iloc, jp_+jloc, kp_+kloc);

    (*BLoc).x =  (*Bx3D)(ip_+iloc, jp_+jloc, kp_+kloc);
    (*BLoc).y =  (*By3D)(ip_+iloc, jp_+jloc, kp_+kloc);
    (*BLoc).z =  (*Bz3D)(ip_+iloc, jp_+jloc, kp_+kloc);


} // END Interpolator3D0Order

void Interpolator3D0Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{

}
