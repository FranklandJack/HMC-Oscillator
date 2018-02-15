#ifndef Lattice1D_hpp
#define Lattice1D_hpp

#include <vector>
#include <random>
#include <iostream>
#include "Ipotential.hpp"


/**
 * \file
 * \brief Class to model a 1D time lattice of a quantum mechanical particle.
 *
 * The class is intended for use in Monte-Carlo simulations of quantum mechanical systems.
 * It is a discrete time lattice in euclidean space that has a mass and a potential to determine 
 * the dynamics of the system. The class uses periodic boundary conditions, so on a lattice of N 
 * sites the (N+1)th site is identified with the site at 0.
 */
class Lattice1D
{
public:
	/// Member variable to hold the size of the lattice.
	int m_size;
	/// Member variable to hold the spacing between lattice sites.
	double m_spacing;
	/// Member variable to hold the mass of the particle.
	double m_mass;
	/// Member variable to hold the position of the particle at each point on the time lattice.
	std::vector<double> m_data;
	/**
	 * Member variable to point to the potential of the system, use if the interface class allows 
	 * user to derive new potentials and make use of virtual functions to have them behave differently.
	 */
	Ipotential *m_potential;

public:
	/**
	 * \brief Construct to create default lattice of zeros.
	 * \param size integer value representing the number of sites on the lattice.
	 * \param spacing floating point value representing the spacing between lattice sites.
	 * \param floating point value representing the mass of the particle.
	 * \param potential Ipotential pointer to the potential of the system.
	 *
	 * This constrictor will create a lattice where the particle position at each time is initially zero.
	 *
	 */
	Lattice1D(int size, double spacing, double mass, Ipotential *potential);

	/**
	 * \brief Construct to create an initialised lattice of values uniformly distributed in some range.
	 * \param size integer value representing the number of sites on the lattice.
	 * \param spacing floating point value representing the spacing between lattice sites.
	 * \param floating point value representing the mass of the particle.
	 * \param potential Ipotential pointer to the potential of the system.
	 * \param generator std::default_rand_engine reference for random number generation.
	 * \param min floating point value representing lower bound in the uniform distribution for the positions.
	 * \param max floating point value representing upper bound in the uniform distribution for the positions.
	 */
	Lattice1D(int size, 
				double spacing, 
				double mass, 
				Ipotential *potential, 
				std::default_random_engine &generator, 
				double min = -1, 
				double max = 1);

	/**
	 * \brief Initialises the lattice with uniform random values in a range set by the user.
	 * \param generator std::default_rand_engine reference for random number generation.
	 * \param min floating point value representing lower bound in the uniform distribution for the positions.
	 * \param max floating point value representing upper bound in the uniform distribution for the positions.
	 */
	void initialise(std::default_random_engine &generator, double min = -1, double max = 1);

	/**
	 * \brief getter method for the size of the lattice
	 * \return integer value representing size of the lattice.
	 */
	int getSize() const;

	/**
	 * \brief getter method for the lattice spacing.
	 * \return floating point value representing the spacing between lattice sites.
	 */
	double getSpacing() const;

	/**
	 * \brief getter method for the mass of the particle.
	 * \return floating point value representing the mass of the particle.
	 */
	double getMass() const;

	/**
	 * \brief getter method for the positions of the particle at each time step
	 * \return std::vector<double> instance representing the position of the particle at each position on the lattice.
	 */
	std::vector<double> getData() const;

	/**
	 * \brief function to calculate the Euclidean action on the lattice.
	 * \return floating point value representing the action of the system.
	 *
	 * Euclidean action is calculated according to the formula S_E = /int (T+V). Where the kinetic term is approximated 
	 * via a forward difference and the integral via a sum.
	 */
	double action() const;

	/**
	 * \brief operator overload for accessing position of particle at a specific site on the lattice.
	 * \param index integer value representing the index of site caller wants position of (can be negative).
	 * \return floating point reference to position of particle at specified index.
	 *
	 * This operator takes into account periodic boundary conditions so can be indexed for negative values or
	 * vales greater than m_size-1.
	 */ 
	double& operator[](int index);

	/**
	 * \brief Constant version of non-constant counter part for use with constant lattices. See non-constant 
	 * version for details.
	 */
	const double&  operator[](int index) const;

	/**
	 * \brief operator overload for streaming the lattice to an output stream.
	 * \param out std::ostream reference that lattice is being output to.
	 * \param lattice constant Lattice1D reference to be output.
	 * \return std::ostream reference so the operation can be chained.
	 */
	friend std::ostream&  operator<<(std::ostream &out, const Lattice1D &lattice);

	/** 
	 * \brief Calculates the correlation function at time t for the system.
	 * \param t integer value representing lattice time at which to calculate the correlation function.
	 * \return floating point value representing the value of the correlation at time t.
	 */
	double correlation(int t) const;

	/**
	 * \brief Calculates the correlation function for a range of values.
	 * \brief t1 integer value representing lower bound for time at which to calculate the correlation function.
	 * \brief t2 integer value representing upper bound for time at which to calculate the correaltion function.
	 * \return std::vector<double> representing the values of the correlation function at and between the lattice
	 * times specified by the caller.
	 */
	std::vector<double> correlation(int t1, int t2) const;

	/**
	 * \brief Calculates the expectation value of position on the lattice.
	 * \return floating point value representing the expectation value of position on the lattice.
	 */
	double meanX() const;

	/**
	 * \brief Calculates the expectation value of position squared on the lattice.
	 * \return floating point value representing the expectation value of position squared on the lattice.
	 */
	double meanXSquared() const;

	/**
	 * \brief Calculates the expectation value of position to the fourth on the lattice.
	 * \return floating point value representing the expectation value of position to the fourth on the lattice.
	 */
	double meanXFourth() const;

};

#endif /* Lattice1D_hpp */