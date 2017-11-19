
#include "Potential2.hpp"

Potential2::Potential2(double lambda, double fSquared):m_fSquared(fSquared),OscPotential(lambda){}

double  Potential2::operator()(double x) const 
{
    return m_lambda*(x*x-m_fSquared)*(x*x-m_fSquared);
}

double Potential2::operator[](double x) const
{
    return 2.0 * m_lambda * (x*x-m_fSquared) * 2*x; 
}
