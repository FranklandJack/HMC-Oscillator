#ifndef Latticefunctions_hpp
#define Latticefunctions_hpp

#include "OscPotential.hpp"
#include <vector>
#include <cmath>

//This is the total action for the lattice, not the normal PE energy
//this is inline with the definition of S being the potential in HMC 
double latticePotentialEnergy(const std::vector<double> &q, double latticeSpacing, double mass, const OscPotential &potentialEnergy);


// function to calculate the HMC kinetic energy term given by p*p/2m
double kineticEnergy(const std::vector<double> &momenta, double mass);


//calculates hmc hamiltonian  according to H(q,p) = p^2/2m + S(q)
double oscillatorHamiltonian(const std::vector<double>& p, const std::vector<double>& q, double latticeSpacing, double mass, const OscPotential &potentialEnergy);


// function to perform leapfrog only!, pass position and momentum vectors by reference and it will update them 

// function to perform leapfrog on whole configuration and return an updated configuration and update the momentum. This function will actually changes the momenta
// but will return an updated configuration, this allows us to keep the original configuration incase we need to reject the new one and since the momenta are discarded 
// anyway we will not need them 
void leapFrog(std::vector<double> &configuration, std::vector<double> &momentum, double latticeSpacing, int leapfrogSteps, double leapfrogStepSize, const OscPotential &potential, double mass = 1.0);

// function to calculate correlation

double correlationFunction(const std::vector<double> &configuration, int n);


#endif