#include "HMCOutput.hpp"

std::ostream& operator<<(std::ostream& out, const HMCOutput &output)
{
	int outputColumnWidth = 30;
	out << "Results..." << '\n';
	out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Acceptance-Rate: " << std::right << output.acceptanceRate << " +/- " << output.acceptanceRateError  << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "S:" << std::right << output.action << " +/- " << output.actionError << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "T:" << std::right << output.kinetic << " +/- " << output.kineticError << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "deltaH:" << std::right << output.dHamiltonian << " +/- " << output.dHamiltonianError << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "exp(-deltaH):" << std::right << output.expdHamiltonian << " +/- " << output.expdHamiltonianError << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "X: " << std::right << output.position << " +/- " << output.positionError << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "PositionIAC: " << std::right << output.positionIAC << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "X^2:" << std::right << output.positionSquared << " +/- " << output.positionSquaredError << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "PositionSquaredIAC: " << std::right << output.positionSquaredIAC << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "X^4:" << std::right << output.positionFourth << " +/- " << output.positionFourthError << '\n'; 
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "PositionFourthIAC: " << std::right << output.positionFourthIAC << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "E_0:" << std::right  << output.gsEnergy << " +/- " << output.gsEnergyError << '\n';

    return out;

}

 	