#ifndef HarmonicPotential_hpp
#define HarmonicPotential_hpp
#include "Ipotential.hpp"

class HarmonicPotential : public Ipotential
{
private:
    // member varible to hold the frequency of the oscillator
    double m_lambda;
    double m_muSquared;

public:
    // we just need a simple constructor to initialise the functor with the values
    HarmonicPotential(double muSquared, double lambda);

    // overload the operator(), then when we act on an x value with the 
    // functor it will return the potential at the displacment
    virtual double  operator()(double x) const override;

    // overload the operator[], then when we act on an x value with the 
    // functor it will return the derivatice of the potential at the displacment
    virtual double operator[](double x) const override;

};

#endif