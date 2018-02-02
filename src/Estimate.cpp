#include "Estimate.hpp"

Estimate::Estimate(double mean, double error): m_mean(mean), m_error(error) {}

double Estimate::getMean() const
{
	return m_mean;
}

void Estimate::setMean(int mean)
{
	m_mean = mean;
}

double Estimate::getError() const
{
	return m_error;
}

void Estimate::setError(int error)
{
	m_error = error;
}


Estimate operator+(const Estimate& a, const Estimate& b)
{
	double mean = a.m_mean + b.m_mean;
	double error = sqrt(a.m_error * a.m_error + b.m_error * b.m_error);

	return Estimate(mean,error);
}
Estimate operator-(const Estimate& a, const Estimate& b)
{
	double mean = a.m_mean - b.m_mean;
	double error = sqrt(a.m_error * a.m_error + b.m_error * b.m_error);

	return Estimate(mean,error);

}
Estimate operator/(const Estimate& a, const Estimate& b)
{
	double mean = a.m_mean / b.m_mean;
	double error = a.m_mean * sqrt(a.m_error/a.m_mean * a.m_error/a.m_mean + b.m_error/b.m_mean * b.m_error/b.m_mean);

	return Estimate(mean,error);
}
Estimate operator*(const Estimate& a, const Estimate& b)
{
	double mean = a.m_mean * b.m_mean;
	double error = a.m_mean * sqrt(a.m_error/a.m_mean * a.m_error/a.m_mean + b.m_error/b.m_mean * b.m_error/b.m_mean);

	return Estimate(mean,error);
}

std::ostream& operator<<(std::ostream& out, const Estimate& a)
{
	out << a.m_mean << " +/- " << a.m_error;
	return out;
}