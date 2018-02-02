#ifndef HMCInput_hpp
#define HMCInput_hpp
#include <iostream>
#include <iomanip>
class HMCInput
{
public:
	enum PotentialType
	{
	    Potential_Harmonic,
	    Potential_Anharmonic,
	    Potential_Octic,
	    Potential_MAX_POTENTIAL
	};

	// Lattice Parameters. 
    int    latticeSize;
    double latticeSpacing;

    // Oscillator Parameters.
    double mass;
    double muSquared;
    double lambda;
    double fSquared;

    // HMC Parameters.
    int    lfStepCount;
    double lfStepSize;

    // Other Parameters. 
    int configCount;
    int burnPeriod;
    int mInterval;

    // Choice of potential.
    PotentialType potentialChoice;

    // Histogram parameters.
    int    numBins;
    double histMaxValue;
    double histMinValue;

    // Tempering parameter sqrt(alpha)
    double temperingParameter;

    // Correlation range.
    int correlationRange;

    friend std::ostream& operator<<(std::ostream& out, const HMCInput &input);


};
#endif /* HMCInput_hpp */