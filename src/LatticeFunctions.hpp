#ifndef Latticefunctions_hpp
#define Latticefunctions_hpp

#include "Ipotential.hpp"
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

double latticeAction(const std::vector<double> &configuration, double latticeSpacing, double mass, const Ipotential *potential);

double kineticEnergy(const std::vector<double> &momenta);

double oscillatorHamiltonian(const std::vector<double>& p, const std::vector<double>& q, double latticeSpacing, double mass, const Ipotential *potential);

void leapFrog(std::vector<double> &configuration, std::vector<double> &momentum, double latticeSpacing, int leapfrogSteps, double leapfrogStepSize, const Ipotential *potential, double mass = 1.0);

void leapFrogTempering(std::vector<double> &configuration, std::vector<double> &momentum, double latticeSpacing, int leapfrogSteps, double leapfrogStepSize, const Ipotential *potential, double mass = 1.0, double alpha = 1.0);

double correlationFunction(const std::vector<double> &configuration, int t);


#endif