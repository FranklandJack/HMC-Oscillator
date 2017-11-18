
// functor to deal with the potential of a 1-D harmonic/enharmonic potential
// by using a functor we can pass it into various functions in our method
#ifndef HarmonicPotential_hpp
#define HarmonicPotential_hpp
#include "OscPotential.hpp"

class OscPotential : public OscPotential
{
private:
    // member varible to hold the frequency of the oscillator
    double m_muSquared;
    // member variavle to hold the enharmonic coeifcient
    double m_lambda;

public:
    // we just need a simple constructor to initialise the functor with the values
    OscPotential(double muSquared = 1.0, double lambda = 0.0);

    // overload the operator(), then when we act on an x value with the 
    // functor it will return the potential at the displacment
    double  operator()(double x) const;

    // overload the operator[], then when we act on an x value with the 
    // functor it will return the derivatice of the potential at the displacment
    double operator[](double x) const;

};

#endif