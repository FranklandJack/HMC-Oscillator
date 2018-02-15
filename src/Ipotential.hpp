
// functor to deal with the potential of a 1-D harmonic/enharmonic potential
// by using a functor we can pass it into various functions in our method
#ifndef Ipotential_hpp
#define Ipotential_hpp
#include <iostream>

class Ipotential
{
public:
    // overload the operator(), then when we act on an x value with the 
    // functor it will return the potential at the displacment
    virtual double  operator()(double x) const = 0;

    // overload the operator[], then when we act on an x value with the 
    // functor it will return the derivatice of the potential at the displacment
    virtual double operator[](double x) const = 0;

    virtual double groundStateEnergy(double meanXSquared, double meanXFourth) const = 0;

};

#endif