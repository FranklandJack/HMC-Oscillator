#ifndef Lattice1D_hpp
#define Lattice1D_hpp
#include <vector>
#include <random>
#include <iostream>
#include "Ipotential.hpp"

class Lattice1D
{
protected:
	int m_size;
	double m_spacing;
	double m_mass;
	std::vector<double> m_data;
	Ipotential *m_potential;

public:
	Lattice1D(int size, double spacing, double mass, Ipotential *potential);

	void   		   initialise(std::default_random_engine &generator, int min = -1, int max = 1);

	int    		   getSize() const;
	double 		   getSpacing() const;
	double 		   getMass() const;
	std::vector<double> getData() const;

	double 		   action() const;


	double& 	   operator[](int index);
	const double&  operator[](int index) const;
	friend std::ostream&  operator<<(std::ostream &out, const Lattice1D &lattice);

	double correlation(int t) const;
	std::vector<double> correlation(int t1, int t2) const;


	double meanX() const;
	double meanXSquared() const;

};

#endif /* Lattice1D_hpp */