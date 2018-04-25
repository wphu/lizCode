#ifndef BACKSCATTERING_H
#define BACKSCATTERING_H


#include <vector>
#include <math.h>

using namespace std;

class Backscattering
{
public:
    Backscattering(){};


    virtual ~Backscattering(){};


    // rnion and reion are number and energy backscattering coefficients, respectively
    // theta is angle of incidence in degree
    // ke is the incident kinetic energy in eV
    virtual void scatter(double &rnion, double &reion, double theta, double ke){};


private:

};



#endif
