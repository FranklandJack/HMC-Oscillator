#include "metropolisUpdate.hpp"
bool metropolisUpdate(double energyBefore, double energyAfter, std::default_random_engine &generator)
{
	double deltaEnergy = energyAfter - energyBefore;
	static std::uniform_real_distribution<double> distribution(0.0,1.0);
	return distribution(generator) <= std::min(1.0,exp(-deltaEnergy));
}