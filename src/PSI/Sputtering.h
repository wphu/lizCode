#ifndef SPUTTERING_H
#define SPUTTERING_H

#include <vector>

class Sputtering
{

public:
    //! Constructor for Collisions between two species
    Sputtering(){};


    virtual ~Sputtering(){};

    virtual void init(){};
    virtual double phy_sput_yield(double ke, double theta){};

private:


};



#endif