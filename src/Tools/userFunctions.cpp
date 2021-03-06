#include <cstdlib>
#include <cmath>
#include "userFunctions.h"
#include <float.h>

//! inverse error function is taken from NIST
double userFunctions::erfinv (double x)
{
    if (x < -1 || x > 1)
        return NAN;

    if (x == 0)
        return 0;

    int  sign_x;
    if (x > 0) {
        sign_x = 1;
    } else {
        sign_x = -1;
        x = -x;
    }

    double r;
    if (x <= 0.686) {
        double x2 = x * x;
        r = x * (((-0.140543331 * x2 + 0.914624893) * x2 + -1.645349621) * x2 + 0.886226899);
        r /= (((0.012229801 * x2 + -0.329097515) * x2 + 1.442710462) * x2 + -2.118377725) * x2 + 1;
    } else {
        double y = sqrt (-log((1 - x) / 2));
        r = (((1.641345311 * y + 3.429567803) * y + -1.62490649) * y + -1.970840454);
        r /= ((1.637067800 * y + 3.543889200) * y + 1);
    }

    r *= sign_x;
    x *= sign_x;

    r -= (erf(r) - x) / (2 / sqrt (M_PI) * exp(-r*r));

    return r;
}

double userFunctions::erfinv2 (double x)
{
    double w, p;
    w = -log((1.0-x)*(1.0+x));

    if ( w < 5.000000 ) {
        w = w - 2.500000;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p*w;
        p = -3.5233877e-06 + p*w;
        p = -4.39150654e-06 + p*w;
        p = 0.00021858087 + p*w;
        p = -0.00125372503 + p*w;
        p = -0.00417768164 + p*w;
        p = 0.246640727 + p*w;
        p = 1.50140941 + p*w;
    } else {
        w = sqrt(w) - 3.000000;
        p = -0.000200214257;
        p = 0.000100950558 + p*w;
        p = 0.00134934322 + p*w;
        p = -0.00367342844 + p*w;
        p = 0.00573950773 + p*w;
        p = -0.0076224613 + p*w;
        p = 0.00943887047 + p*w;
        p = 1.00167406 + p*w;
        p = 2.83297682 + p*w;
    }
    return p*x;
}


int userFunctions::isFiniteNumber(double d)
{
    return (d<=DBL_MAX&&d>=-DBL_MAX);
}

//generate a gaussian random number R: mean=0, mean square=1.0
//Then you can get another gaussian random number R2 = R*V + E: mean=E, mean square=V;
double userFunctions::gaussrand()
{
    static double V1, V2, S;
    	    static int phase = 0;
    	    double X;

    	    if(phase == 0) {
        		do {
        		    double U1 = (double)rand() / RAND_MAX;
        		    double U2 = (double)rand() / RAND_MAX;

        		    V1 = 2 * U1 - 1;
        		    V2 = 2 * U2 - 1;
        		    S = V1 * V1 + V2 * V2;
        		} while(S >= 1 || S == 0);

        		X = V1 * sqrt(-2 * log(S) / S);
    	    }
            else
            {
                X = V2 * sqrt(-2 * log(S) / S);
            }

    	    phase = 1 - phase;

    	    return X;
}
