#include "Potential3.hpp"

Potential3::Potential3(double lambda, double fSquared):m_lambda(lambda), m_fSquared(fSquared){}

double Potential3::operator()(double x) const
{
	return m_lambda * (x*x - m_fSquared)*(x*x - m_fSquared)*(x*x - m_fSquared)*(x*x - m_fSquared);
}

double Potential3::operator[](double x) const
{
	return 4.0 * m_lambda * (x*x - m_fSquared) * (x*x - m_fSquared) * (x*x - m_fSquared) * 2.0 * x;
}