#ifndef OcticPotential_hpp
#define OcticPotential_hpp
#include "Ipotential.hpp"

class OcticPotential : public Ipotential
{

private:
	double m_lambda;
	double m_fSquared;

public:
	OcticPotential(double lambda, double fSquared);

	virtual double operator()(double x) const override;

	virtual double operator[](double x) const override;

};

#endif /* OcticPotential_hpp */