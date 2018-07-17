#include "Interpolator3D1Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D1Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D1Order::Interpolator3D1Order(PicParams &params, SmileiMPI *smpi) : Interpolator3D(params, smpi)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    dz_inv_ = 1.0/params.cell_length[2];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
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
    coeffxp_[0] = 1.0 - delta;
    coeffxp_[1] = delta;


    delta   = ypn - (double)jp_;
    coeffyp_[0] = 1.0 - delta;
    coeffyp_[1] = delta;

    delta   = zpn - (double)kp_;
    coeffzp_[0] = 1.0 - delta;
    coeffzp_[1] = delta;

    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    kp_ = kp_ - k_domain_begin;


    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    (*ELoc).x =  compute( coeffxp_, coeffyp_, coeffzp_, Ex3D, ip_, jp_, kp_);

    // -------------------------
    // Interpolation of Ey^(p,d)
    // -------------------------
    (*ELoc).y = compute( coeffxp_, coeffyp_, coeffzp_, Ey3D, ip_, jp_, kp_);

    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    (*ELoc).z = compute( coeffxp_, coeffyp_, coeffzp_, Ez3D, ip_, jp_, kp_);

    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    (*BLoc).x = compute( coeffxp_, coeffyp_, coeffzp_, Bx3D, ip_, jp_, kp_);

    // -------------------------
    // Interpolation of By^(d,p)
    // -------------------------
    (*BLoc).y = compute( coeffxp_, coeffyp_, coeffzp_, By3D, ip_, jp_, kp_);

    // -------------------------
    // Interpolation of Bz^(d,d)
    // -------------------------
    (*BLoc).z = compute( coeffxp_, coeffyp_, coeffzp_, Bz3D, ip_, jp_, kp_);

} // END Interpolator3D1Order

void Interpolator3D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{

}
