#include "HarmonicPotential.hpp"

HarmonicPotential::HarmonicPotential(double muSquared, double lambda):m_muSquared(muSquared),m_lambda(lambda){}

double  HarmonicPotential::operator()(double x) const 
{
    return 0.5 * m_muSquared * x * x + m_lambda * x * x * x * x;
}

double HarmonicPotential::operator[](double x) const
{
    return m_muSquared * x + 4.0 * m_lambda * x * x * x; 
}
