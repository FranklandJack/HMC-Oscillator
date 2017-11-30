#include "LatticeFunctions.hpp"




//This is the total action for the lattice, not the normal PE energy
//this is inline with the definition of S being the potential in HMC 
double latticePotentialEnergy(const std::vector<double> &configuration, double latticeSpacing, double mass, const OscPotential &potential)
{
    //create sum value which we add each term to
    double totalPotentialEnergy = 0.0;

    int latticeSize = configuration.size();

    for(int i = 0; i < latticeSize-1; ++i)
    {
        totalPotentialEnergy += (0.5 * mass * (configuration[i+1]-configuration[i]) * (configuration[i+1]-configuration[i])) / latticeSpacing + latticeSpacing * potential(configuration[i]);
    }

    totalPotentialEnergy += (0.5 * mass * (configuration[0]-configuration[latticeSize-1]) * (configuration[0]-configuration[latticeSize-1])) / latticeSpacing + latticeSpacing * potential(configuration[latticeSize-1]);

    
    return totalPotentialEnergy;
}

// function to calculate the HMC kinetic energy term given by p*p/2m
double kineticEnergy(const std::vector<double> &momentum, double mass)
{
    double totalKE = 0;

    for(const auto& p : momentum)
    {
        totalKE += p*p;
    }

    return totalKE/(2.0 * mass);

}

//calculates hmc hamiltonian  according to H(q,p) = p^2/2m + S(q)
double oscillatorHamiltonian(const std::vector<double>& p, const std::vector<double>& q, double latticeSpacing, double mass, const OscPotential &potential)
{
    return kineticEnergy(p,mass) + latticePotentialEnergy(q, latticeSpacing, mass, potential);
}


void leapFrog(std::vector<double> &configuration, std::vector<double> &momentum, double latticeSpacing, int lfStepCount, double lfStepSize, const OscPotential &potential, double mass)
{
    int latticeSize = configuration.size();

    // Intial half step in momentum.
    // First lattice site has special neighbour conditions. 
    momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * potential[configuration[0]]) * (lfStepSize/2.0);

    // All remaining sites except last have normal neighbours.
    for(int i = 1; i <= latticeSize - 2; ++i)
    {
        momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * potential[configuration[i]]) * (lfStepSize/2.0);
    }

    // Final site has special neighbour conditions.
    momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * potential[configuration[latticeSize-1]]) * (lfStepSize/2.0);


    // Full step in position. 
    for(int i = 0; i < latticeSize; ++i)
    {
        configuration[i] = configuration[i] + lfStepSize * momentum[i] / mass;
    }

    // N-1 Full steps in momentum and position.
    for(int n = 1; n < lfStepCount; ++n)
    {
        // First lattice site has special neighbour conditions. 
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize-1] ) + latticeSpacing * potential[configuration[0]]) * (lfStepSize);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * potential[configuration[i]]) * (lfStepSize);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * potential[configuration[latticeSize-1]]) * (lfStepSize);


        // Full step in position. 
        for(int i = 0; i < latticeSize; ++i)
        {
            configuration[i] = configuration[i] + lfStepSize * momentum[i] / mass;
        }

    }

    // Final half step in momentum.

    // First lattice site has special neighbour conditions. 
    momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize-1] ) + latticeSpacing * potential[configuration[0]]) * (lfStepSize/2.0);

    // All remaining sites except last have normal neighbours.
    for(int i = 1; i <= latticeSize - 2; ++i)
    {
        momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * potential[configuration[i]]) * (lfStepSize/2.0);
    }

    // Final site has special neighbour conditions.
    momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * potential[configuration[latticeSize-1]]) * (lfStepSize/2.0);

}

double correlationFunction(const std::vector<double> &configuration, int t)
{ 
    
    double sum = 0;
    double normalisation = 0;
    for(int i = 0; i < configuration.size(); ++i)
    {
        /*std::vector<double>::iterator ith = copyConfig.begin()+i;
        std::rotate(copyConfig.begin(),ith, copyConfig.end());*/
        sum += configuration[i]*configuration[(t+i)%(configuration.size())];
        normalisation += configuration[i]*configuration[i];
    }

    return sum/normalisation;
}

