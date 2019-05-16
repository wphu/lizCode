#include "Collisions1D_Recombination_TB.h"
#include "SmileiMPI.h"
#include "Field2D.h"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_Recombination_TB::Collisions1D_Recombination_TB(PicParams& params, vector<Species*>& vecSpecies, SmileiMPI* smpi,
    unsigned int n_col,
    vector<unsigned int> sg1,
    vector<unsigned int> sg2,
    vector<unsigned int> sg3,
    string CS_fileName)
: Collisions1D(params)
{
    n_collisions    = n_col;

    // reaction: e + e + H+ = e + H + hv
    // species_group1: e
    // species_group2: H+
    // species_group3: H
    species_group1  = sg1;
    species_group2  = sg2;
    species_group3  = sg3;
    //crossSection_fileName = CS_fileName;

    // Calculate total number of bins
    nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;

    //readCrossSection();


    // only for hydrogen isotope
    nmax = 20;

    // calculate E_bound, formular (14)
    double Ry = 13.6 * const_e;
    double numerator = 0;
    double denominator = 0;
    for(int i = 1; i <= nmax; i++)
    {
        numerator += i * i * i * i;
    }
    for(int i = 1; i <= nmax; i++)
    {
        denominator += i * i * i * i * i * i;
    }
    E_bound = Ry * numerator / denominator;

}

Collisions1D_Recombination_TB::~Collisions1D_Recombination_TB()
{

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_Recombination_TB::collide(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime)
{
    vector<unsigned int> *sg1, *sg2, *sg3;

    vector<int> index1, index2;
    vector<int> n1, n2;
    vector<double> density1, density2;
    double n1_max, n2_max;
    vector<double> momentum_unit(3, 0.0), momentum_temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2, bmin3, bmax3;
    unsigned int npairs; // number of pairs of macro-particles
    unsigned int i11, i12, i2, i3;
    Species   *s1, *s2, *s3;
    Particles *p1, *p2, *p3;
    double m1, m2, m3, m12, W1, W2, W3;

    double sigma_cr, sigma_cr_max, ke11, ke12, ke11_primary, ke_secondary,
           ran, P_collision;
    double v11_square, v11_magnitude, v11_magnitude_primary, v12_square, v12_magnitude;
    double ke12_post, v12_post;

    int iBin_global;
    double ke_radiative;

    Diagnostic1D *diag1D = static_cast<Diagnostic1D*>(diag);

    sg1 = &species_group1;
    sg2 = &species_group2;
    sg3 = &species_group3;

    // electons                         ions                            atoms
    s1 = vecSpecies[(*sg1)[0]];      s2 = vecSpecies[(*sg2)[0]];        s3 = vecSpecies[(*sg3)[0]];
    p1 = &(s1->particles);           p2 = &(s2->particles);             p3 = &(s3->particles);
    m1 = s1->species_param.mass;     m2 = s2->species_param.mass;       m3 = s3->species_param.mass;
    W1 = p1->weight(0);              W2 = p2->weight(0);                //W3 = p3->weight(0);
    bmin1 = s1->bmin;                bmin2 = s2->bmin;                  bmin3 = s3->bmin;
    bmax1 = s1->bmax;                bmax2 = s2->bmax;                  bmax3 = s3->bmax;

    
    count_of_particles_to_insert_s3.resize(nbins);
    count_of_particles_to_erase_s1.resize(nbins);
    count_of_particles_to_erase_s2.resize(nbins);
    for(int ibin=0; ibin<nbins; ibin++)
    {
        count_of_particles_to_insert_s3[ibin] = 0;
        count_of_particles_to_erase_s1[ibin] = 0;
        count_of_particles_to_erase_s2[ibin] = 0;
    }

    n1.resize(nbins);
    density1.resize(nbins);
    n1_max = 0.0;
    for(unsigned int ibin=0 ; ibin<nbins ; ibin++)
    {
        n1[ibin] = bmax1[ibin] - bmin1[ibin];
        density1[ibin] = n1[ibin] * W1;
    }

    n2.resize(nbins);
    density2.resize(nbins);
    n2_max = 0.0;
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++)
    {
        n2[ibin] = bmax2[ibin] - bmin2[ibin];
        density2[ibin] = n2[ibin] * W2;
        if( density2[ibin] > n2_max ) { n2_max = density2[ibin]; };
    }

    totNCollision = 0;
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {
        //MESSAGE("nbins000"<<"  "<<ibin<<"  "<<bmin2[ibin]<<" "<<bmax2[ibin]);
        // calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize( n1[ibin] );

        for(int iPart = 0; iPart < n1[ibin]; iPart++)
        {
            index1[iPart] = bmin1[ibin] + iPart;
        }
        random_shuffle(index1.begin(), index1.end());

        // calculate the particle number of species2 in each cell, and the indexs of particles in the cell
        index2.resize( n2[ibin] );

        for(int iPart = 0; iPart < n2[ibin]; iPart++)
        {
            index2[iPart] = bmin2[ibin] + iPart;
        }
        random_shuffle(index2.begin(), index2.end());

        npairs = n1[ibin] / 2;
        if(n1[ibin] % 2 == 1)
        {
            npairs++;
        }
        if(npairs > n2[ibin])
        {
            npairs = n2[ibin];
        }

        for(int i = 0; i < npairs; i++)
        {
            i11 = index1[2*i];
            i12 = index1[2*i+1];

            // Collision  only erase the i11 particle
            if(n1[ibin] % 2 == 1 && i == npairs - 1)
            {
                i11 = index1[2*i];
                i12 = index2[2*(i-1)+1];
            }

            v11_square = pow(p1->momentum(0,i11),2) + pow(p1->momentum(1,i11),2) + pow(p1->momentum(2,i11),2);
            v11_magnitude = sqrt(v11_square);
            // kinetic energy of i11 electron
            ke11 = 0.5 * m1 * v12_square;

            v12_square = pow(p1->momentum(0,i12),2) + pow(p1->momentum(1,i12),2) + pow(p1->momentum(2,i12),2);
            v12_magnitude = sqrt(v12_square);
            // kinetic energy of i12 electron
            ke12 = 0.5 * m1 * v12_square;

            // post-collision energy of i11 electron
            ke12_post = ke11 + ke12 + E_bound;
            v12_post = sqrt( 2.0 * ke12_post / m1 );

            ke_radiative = -E_bound;

            P_collision = 1.0 - exp( sqrt(v11_magnitude * v12_magnitude) * cross_section(ke11, ke12)
                          * sqrt(density1[ibin] * density2[ibin]) * timestep );

            // Generate a random number between 0 and 1
            double ran_p = (double)rand() / RAND_MAX;
            if(ran_p < P_collision){
                // erase i11 electron and ion
                count_of_particles_to_erase_s1[ibin]++;
                indexes_of_particles_to_erase_s1.push_back(i11);
                count_of_particles_to_erase_s2[ibin]++;
                indexes_of_particles_to_erase_s2.push_back(i2);

                // Calculate the scatter velocity of i12 electron
                momentum_unit[0] = p1->momentum(0,i12) / v12_magnitude;
                momentum_unit[1] = p1->momentum(1,i12) / v12_magnitude;
                momentum_unit[2] = p1->momentum(2,i12) / v12_magnitude;
                calculate_scatter_velocity(v12_post, m1, m2, momentum_unit, momentum_temp);
                p1->momentum(0,i12) = momentum_temp[0];
                p1->momentum(1,i12) = momentum_temp[1];
                p1->momentum(2,i12) = momentum_temp[2];

                // Copy the particle of species2 to species3, and change the charge
                p2->cp_particle(i2, new_particles3);
                count_of_particles_to_insert_s3[ibin]++;
                idNew = new_particles3.size() - 1;
                new_particles3.charge(idNew) = s3->species_param.charge;
                totNCollision++;

                iBin_global = smpi->getDomainLocalMin(0) / params.cell_length[0] + ibin;
                diag1D->radiative_energy_collision[n_collisions][iBin_global] += ke_radiative;
            }
        }

    } // end loop on bins

    s1->erase_particles_from_bins(indexes_of_particles_to_erase_s1);
    indexes_of_particles_to_erase_s1.clear();

    s2->erase_particles_from_bins(indexes_of_particles_to_erase_s2);
    indexes_of_particles_to_erase_s2.clear();

    s3->insert_particles_to_bins(new_particles3, count_of_particles_to_insert_s3);
    new_particles3.clear();
}

double Collisions1D_Recombination_TB::maxCV(Particles* particles, double eMass)
{
    int nPart = particles->size();
    double v_square;
    double v_magnitude;
    double ke;
    double maxCrossSectionV = 0.0;
    double crossSectionV;

    for(unsigned int iPart = 0; iPart < nPart; iPart++)
    {
        v_square = particles->momentum(0,iPart) * particles->momentum(0,iPart) +
                          particles->momentum(1,iPart) * particles->momentum(1,iPart) +
                          particles->momentum(2,iPart) * particles->momentum(2,iPart);
        v_magnitude = sqrt(v_square);
        // ke is energy (eV)
        ke = 0.5 * eMass * v_square / const_e;
        crossSectionV = v_magnitude * interpCrossSection( ke );
        if(crossSectionV > maxCrossSectionV) {maxCrossSectionV = crossSectionV;};
    }
    return maxCrossSectionV;
}

//>the method is eqution (11) from the ref: a Monte Carlo collision model for the particle in cell method: applications to
//>argon and oxygen discharges.
//>and the code is transformed from C.F. Sang's fortran code
void Collisions1D_Recombination_TB::calculate_scatter_velocity( double v_magnitude, double mass1, double mass2,
                                                                vector<double>& momentum_unit, vector<double>& momentum_temp)
{
    double up1, up2, up3;
    double r11, r12, r13, r21, r22, r23, r31, r32, r33;
    double mag;

    double ra = (double)rand() / RAND_MAX;
    double costheta = 1.0 - 2.0 * ra;
    double sintheta = sqrt(1.0 - abs(costheta * costheta) );

    ra = (double)rand() / RAND_MAX;
    double pi = 3.1415926;
    double phi = 2.0 * pi * ra;
    double cosphi = cos(phi);
    double sinphi = sin(phi);

    double ve=v_magnitude*sqrt(1.0-2.0*mass1*(1.0-costheta)/mass2);

    r13 = momentum_unit[0];
    r23 = momentum_unit[1];
    r33 = momentum_unit[2];
    if(r33 == 1.0 ){
        up1= 0.;
        up2= 1.;
        up3= 0.;
    }
    else{
        up1= 0.;
        up2= 0.;
        up3= 1.;
    }

    r12 = r23 * up3 - r33 * up2;
    r22 = r33 * up1 - r13 * up3;
    r32 = r13 * up2 - r23 * up1;
    mag = sqrt(r12 * r12 + r22 * r22 + r32 * r32);
    r12 = r12 / mag;
    r22 = r22 / mag;
    r32 = r32 / mag;
    r11 = r22 * r33 - r32 * r23;
    r21 = r32 * r13 - r12 * r33;
    r31 = r12 * r23 - r22 * r13;
    momentum_temp[0] = ve * (r11 * sintheta * cosphi + r12 * sintheta * sinphi + r13 * costheta);
    momentum_temp[1] = ve * (r21 * sintheta * cosphi + r22 * sintheta * sinphi + r23 * costheta);
    momentum_temp[2] = ve * (r31 * sintheta * cosphi + r32 * sintheta * sinphi + r33 * costheta);


}



double Collisions1D_Recombination_TB::cross_section(double ke1, double ke2)
{
    double cs = 0.0;
    for(int i = 0; i <= nmax; i++)
    {
        cs += cross_section_each(ke1, ke2, i);
    }
    return cs;
}

double Collisions1D_Recombination_TB::cross_section_sub(double E1, double E2, int n)
{
    double A = 0.3440;
    double a[4]  = {-0.014353, 0.75206, -0.29548, 0.056884};
    double Ry = 13.6 * const_e;
    double En = Ry / (n*n);
    double e1 = 1.0 + E1 / En;
    double e2 = 1.0 + E2 / En;
    double e12 = e1 - e2;

    double cs = 1.0 / (e1 * e1) + 1.0 / (e12 * e12) - 1.0 / (e1 * e12);
    // this formula below may be wrong in the ref
    // cs = 1.0 / (e1 * e1) + 1.0 / (e12 * e12) - 1.0 / (e1 * e12) - 1.0 / (e1 * e2)
    double sum0 = 0.0;
    for(int i = 0; i < 4; i++)
    {
        sum0 += (a[i] / pow(e2, i));
    }
        
    cs += (sum0 * log(e1) / (n * pow(e2, 3)));
    cs *= (A * pow(n, 4) / (En * e1));
    return cs;
}


double Collisions1D_Recombination_TB::cross_section_each(double E1,double E2, int n)
{
    double Ry = 13.6 * const_e;
    double En = Ry / (n*n);
    double me = 9.109382616e-31;
    // a0: Bohr radius
    double a0 = 5.2917721067e-11;
    double E1_prime = E1 + E2 + En;
    double cs = 0.0;

    if(E1 >= E2)
    {
        cs = 0.5 * cross_section_sub(E1_prime, E2, n);
    }
    else
    {
        cs = 0.5 * cross_section_sub(E1_prime, E1, n);
    }   
    cs *= 2.0;

    double coefficient = 4.0 * pow(Ry, 1.5) * sqrt(2.0 * me) * n * n * const_pi * const_pi * pow(a0, 3) * E1_prime / (E1 *E2);
    cs *= coefficient;
    return cs;
}