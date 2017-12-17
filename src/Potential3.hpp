#ifndef Potential3_hpp
#define Potential3_hpp
#include "Ipotential.hpp"

class Potential3 : public Ipotential
{

private:
	double m_lambda;
	double m_fSquared;

public:
	Potential3(double lambda, double fSquared);

	double operator()(double x) const;

	double operator[](double x) const;

};

#endif /* Potential3_hpp */