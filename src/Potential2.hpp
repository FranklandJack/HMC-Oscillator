#ifndef Potential2_hpp
#define Potential2_hpp
#include "OscPotential.hpp"

class Potential2 : OscPotential
{
private:
    // member varible to hold the frequency of the oscillator
    double m_fSquared;

public:
    // we just need a simple constructor to initialise the functor with the values
    Potential2(double lambda, double fSquared);

    // overload the operator(), then when we act on an x value with the 
    // functor it will return the potential at the displacment
    virtual double  operator()(double x) const override;

    // overload the operator[], then when we act on an x value with the 
    // functor it will return the derivatice of the potential at the displacment
    virtual double operator[](double x) const override;

};

#endif