#include "OscPotential.hpp"

OscPotential::OscPotential(double muSquared, double lambda):m_lambda(lambda), m_muSquared(muSquared){}

double  OscPotential::operator()(double x) const 
{
    return 0.5 * m_muSquared * x * x + m_lambda * x * x * x * x;
}

double OscPotential::operator[](double x) const
{
    return m_muSquared * x + 4.0 * m_lambda * x * x * x;
}