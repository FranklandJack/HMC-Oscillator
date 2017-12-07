#include "Potential1.hpp"

Potential1::Potential1(double muSquared, double lambda):m_muSquared(muSquared),m_lambda(lambda){}

double  Potential1::operator()(double x) const 
{
    return 0.5 * m_muSquared * x * x + m_lambda * x * x * x * x;
}

double Potential1::operator[](double x) const
{
    return m_muSquared * x + 4.0 * m_lambda * x * x * x; 
}
