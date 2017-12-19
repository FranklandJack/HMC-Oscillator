#include "AnharmonicPotential.hpp"

AnharmonicPotential::AnharmonicPotential(double lambda, double fSquared):m_fSquared(fSquared),m_lambda(lambda){}

double  AnharmonicPotential::operator()(double x) const 
{
    return m_lambda*(x*x-m_fSquared)*(x*x-m_fSquared);
}

double AnharmonicPotential::operator[](double x) const
{
    return 2.0 * m_lambda * (x*x-m_fSquared) * 2*x; 
}
