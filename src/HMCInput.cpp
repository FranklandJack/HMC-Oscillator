#include "HMCInput.hpp"

std::ostream& operator<<(std::ostream& out, const HMCInput &input)
{
	// Tell user their input values to check they are correct.
    int outputColumnWidth = 30;
    out << "Input Parameters..." << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Lattice-Size: " << std::right << input.latticeSize << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Lattice-Spacing: " << std::right << input.latticeSpacing << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "LeapFrog-Step-Size: " << std::right << input.lfStepSize << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "LeapFrog-Steps " << std::right << input.lfStepCount<< '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Configurations: " << std::right << input.configCount << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Burn-Period: " << std::right << input.burnPeriod << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Measurement-Interval: " << std::right << input.mInterval << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Mass: " << std::right << input.mass << '\n';
    // Depending on which potential was used report correct parameters.
    switch(input.potentialChoice)
    {
        case HMCInput::Potential_Harmonic: 
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "V(x): " << std::right << "0.5.mu^2.x + lambda.x^4" << '\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "mu^2: " << std::right << input.muSquared << '\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Lambda: " << std::right << input.lambda << '\n';
            break;

        case HMCInput::Potential_Anharmonic: 
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "V(x): " << std::right << "lambda.(x^2-f^2)^2" << '\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "f^2: " << std::right << input.fSquared << '\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Lambda: " << std::right << input.lambda << '\n';
            break;

        case HMCInput::Potential_Octic: 
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "V(x): " << std::right << "lambda.(x^2-f^2)^4" << '\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "f^2: " << std::right << input.fSquared << '\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Lambda: " << std::right << input.lambda << '\n';
            break;
        default:
            out << "No potential selected, exiting program...";
            exit(1);
    }

    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Tempering-Parameter: " << std::right << input.temperingParameter << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Correlation-Range: " << std::right << input.correlationRange << '\n';
    return out;
}