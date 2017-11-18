
// functor to deal with the potential of a 1-D harmonic/enharmonic potential
// by using a functor we can pass it into various functions in our method
#ifndef OscPotential_hpp
#define OscPotential_hpp

class OscPotential
{
protected:
    // member variavle to hold the enharmonic coeifcient
    double m_lambda, m_muSquared;

public:
    // we just need a simple constructor to initialise the functor with the values
    OscPotential(double lambda, double muSquared);

    // overload the operator(), then when we act on an x value with the 
    // functor it will return the potential at the displacment
    virtual double  operator()(double x) const;

    // overload the operator[], then when we act on an x value with the 
    // functor it will return the derivatice of the potential at the displacment
    virtual double operator[](double x) const;

};

#endif