#ifndef HMCOutput_hpp
#define HMCOutput_hpp
#include <iostream>
#include <iomanip>

class HMCOutput
{
public:
	double acceptanceRate;
	double acceptanceRateError;

	double action;
	double actionError;

	double kinetic;
	double kineticError;

	double dHamiltonian;
	double dHamiltonianError;

	double expdHamiltonian;
	double expdHamiltonianError;

	double position;
	double positionError;
	double positionIAC;

	double positionSquared;
	double positionSquaredError;
	double positionSquaredIAC;

	double positionFourth;
	double positionFourthError;
	double positionFourthIAC;

	double gsEnergy;
	double gsEnergyError;
	
	friend std::ostream& operator<<(std::ostream& out, const HMCOutput &output);
};

#endif /* HMCOutput_hpp */