#ifndef HMCLattice1D_hpp
#define HMCLattice1D_hpp
#include <vector>
#include <random>
#include "Lattice1D.hpp"

/**
 * \file
 * \brief Class to model a Hybrid Monte-Carlo lattice.
 *
 *  It inherits from the Lattice1D class and adds a momenta to go with each position variable. 
 * The whole system can then be evolved in time using Hamiltonian dynamics and a leapfrog method.
 */
class HMCLattice1D : public Lattice1D
{
private:
	/// Member variable to represent the conjugate momenta of the HMC algorithm.
	std::vector<double> m_momenta;

public:
	/**
	 * \brief Construct to create default lattice where position and momenta are zeros.
	 * \param size integer value representing the number of sites on the lattice.
	 * \param spacing floating point value representing the spacing between lattice sites.
	 * \param floating point value representing the mass of the particle.
	 * \param potential Ipotential pointer to the potential of the system.
	 *
	 * This constrictor will create a lattice where the particle position and momenta at each 
	 * time is initially zero.
	 *
	 */
	HMCLattice1D(int size, double spacing, double mass, Ipotential *potential);

	/**
	 * \brief Randomises momenta by drawing them from a Gaussian distribution.
	 * \param std::deafault_random_engine reference for generating random numbers.
	 * 
	 * As is standard practice in HMC momenta a drawn from a distribution of mean 0
	 * and variance 1.
	 */
	void randomiseMomenta(std::default_random_engine& generator);

	/**
	 * \brief Calculates the kinetic energy corresponding to the conjugate momenta.
	 * \return floating point value representing the kinetic energy.
	 *
	 * Kinetic energy is calculated according to the formula p*p/2, this is not the 
	 * kinetic energy of the particle and its velocity, but the kinetic energy for the 
	 * HMC algorithm.
	 */
	double kineticEnergy() const;

	/**
	 * \brief Calculates the Hamiltonian for the HMC system.
	 * \return floating point value representing the Hamiltonian.
	 *
	 * The HMC Hamiltonian is given by the sum of the Euclidean action of the system and
	 * the kinetic HMC term.
	 */
	double hamiltonian() const;


	void leapFrog(int lfStepCount, double lfStepSize, double alpha = 1.0);

	bool leapFrogUpdate(std::default_random_engine& generator, int lfStepCount, double lfStepSize, double alpha = 1.0);

	friend void leapFrog(const HMCLattice1D &currentLattice, HMCLattice1D &updatedLattice, int lfStepCount, double lfStepSize, double alpha);

};

#endif /* HMCLattice1D_hpp */