#include "OscPotential.hpp"

OscPotential::OscPotential(double lambda):m_lambda(lambda){}

double  OscPotential::operator()(double x) const 
{
    return m_lambda * x * x * x * x;
}

double OscPotential::operator[](double x) const
{
    return  4.0 * m_lambda * x * x * x;
}