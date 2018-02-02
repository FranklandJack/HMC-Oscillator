#ifndef HMCLattice1D_hpp
#define HMCLattice1D_hpp
#include <vector>
#include <random>
#include "Lattice1D.hpp"

class HMCLattice1D : public Lattice1D
{
private:
	std::vector<double> m_momenta;

public:
	HMCLattice1D(int size, double spacing, double mass, Ipotential *potential);

	void randomiseMomenta(std::default_random_engine& generator);

	double kineticEnergy() const;

	double hamiltonian() const;

	void leapFrog(int lfStepCount, double lfStepSize, double alpha = 1.0);

	bool leapFrogUpdate(std::default_random_engine& generator, int lfStepCount, double lfStepSize, double alpha = 1.0);
	
	
};

#endif /* HMCLattice1D_hpp */