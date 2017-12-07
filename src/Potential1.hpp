#ifndef Potential1_hpp
#define Potential1_hpp
#include "Ipotential.hpp"

class Potential1 : public Ipotential
{
private:
    // member varible to hold the frequency of the oscillator
    double m_lambda;
    double m_muSquared;

public:
    // we just need a simple constructor to initialise the functor with the values
    Potential1(double muSquared, double lambda);

    // overload the operator(), then when we act on an x value with the 
    // functor it will return the potential at the displacment
    virtual double  operator()(double x) const override;

    // overload the operator[], then when we act on an x value with the 
    // functor it will return the derivatice of the potential at the displacment
    virtual double operator[](double x) const override;

};

#endif