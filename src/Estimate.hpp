#ifndef Estimate_hpp
#define Estimate_hpp
#include <iostream>
#include <cmath>	
class Estimate
{
private:
	double m_mean;
	double m_error;

public:
	Estimate(double,double);
	double getMean() const;
	void setMean(int);
	double getError() const;
	void setError(int);
	friend Estimate operator+(const Estimate&, const Estimate&);
	friend Estimate operator-(const Estimate&, const Estimate&);
	friend Estimate operator/(const Estimate&, const Estimate&);
	friend Estimate operator*(const Estimate&, const Estimate&);
	friend std::ostream& operator<<(std::ostream&, const Estimate&);
};

#endif /* Estimate_hpp */