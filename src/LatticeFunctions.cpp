#include "LatticeFunctions.hpp"




//This is the total action for the lattice, not the normal PE energy
//this is inline with the definition of S being the potential in HMC 
double latticePotentialEnergy(const std::vector<double> &configuration, double latticeSpacing, double mass, const Ipotential *potential)
{
    //create sum value which we add each term to
    double totalPotentialEnergy = 0.0;

    int latticeSize = configuration.size();

    for(int i = 0; i < latticeSize-1; ++i)
    {
        totalPotentialEnergy += (0.5 * mass * (configuration[i+1]-configuration[i]) * (configuration[i+1]-configuration[i])) / latticeSpacing + latticeSpacing * (*potential)(configuration[i]);
    }

    totalPotentialEnergy += (0.5 * mass * (configuration[0]-configuration[latticeSize-1]) * (configuration[0]-configuration[latticeSize-1])) / latticeSpacing + latticeSpacing * (*potential)(configuration[latticeSize-1]);

    
    return totalPotentialEnergy;
}

// function to calculate the HMC kinetic energy term given by p*p/2m
double kineticEnergy(const std::vector<double> &momentum)
{
    double totalKE = 0;

    for(const auto& p : momentum)
    {
        totalKE += p*p;
    }

    return totalKE/(2.0);

}

//calculates hmc hamiltonian  according to H(q,p) = p^2/2m + S(q)
double oscillatorHamiltonian(const std::vector<double>& p, const std::vector<double>& q, double latticeSpacing, double mass, const Ipotential *potential)
{
    return kineticEnergy(p) + latticePotentialEnergy(q, latticeSpacing, mass, potential);
}


void leapFrog(std::vector<double> &configuration, std::vector<double> &momentum, double latticeSpacing, int lfStepCount, double lfStepSize, const Ipotential *potential, double mass)
{
    int latticeSize = configuration.size();

    // Intial half step in momentum.
    // First lattice site has special neighbour conditions. 
    momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

    // All remaining sites except last have normal neighbours.
    for(int i = 1; i <= latticeSize - 2; ++i)
    {
        momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
    }

    // Final site has special neighbour conditions.
    momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);


    // Full step in position. 
    for(int i = 0; i < latticeSize; ++i)
    {
        configuration[i] = configuration[i] + lfStepSize * momentum[i];
    }

    // N-1 Full steps in momentum and position.
    for(int n = 1; n < lfStepCount; ++n)
    {
        // First lattice site has special neighbour conditions. 
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize-1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize);


        // Full step in position. 
        for(int i = 0; i < latticeSize; ++i)
        {
            configuration[i] = configuration[i] + lfStepSize * momentum[i];
        }

    }

    // Final half step in momentum.

    // First lattice site has special neighbour conditions. 
    momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize-1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

    // All remaining sites except last have normal neighbours.
    for(int i = 1; i <= latticeSize - 2; ++i)
    {
        momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
    }

    // Final site has special neighbour conditions.
    momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);

}

void leapFrogTempering(std::vector<double> &configuration, std::vector<double> &momentum, double latticeSpacing, int lfStepCount, double lfStepSize, const Ipotential *potential, double mass, double alpha)
{
    int latticeSize = configuration.size();

    for(int n = 0; n < lfStepCount/2; ++n)
    {

        // Intial half step in momentum.
        // First lattice site has special neighbour conditions. 

        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : momentum)
        {
            p *= alpha;
        }
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);


        // Full step in position. 
        for(int i = 0; i < latticeSize; ++i)
        {
            configuration[i] = configuration[i] + lfStepSize * momentum[i];
        }

        // half step in momentum.
        // First lattice site has special neighbour conditions. 
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);

        for(auto& p : momentum)
        {
            p *= alpha;
        }


    }

    if(0 != lfStepCount%2)
    {
        // Intial half step in momentum.
        // First lattice site has special neighbour conditions. 

        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : momentum)
        {
            p *= alpha;
        }
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);


        // Full step in position. 
        for(int i = 0; i < latticeSize; ++i)
        {
            configuration[i] = configuration[i] + lfStepSize * momentum[i];
        }

        // half step in momentum.
        // First lattice site has special neighbour conditions. 
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);

        for(auto& p : momentum)
        {
            p /= alpha;
        }

    }

    for(int n = 0; n < lfStepCount/2; ++n)
    {

        // Intial half step in momentum.
        // First lattice site has special neighbour conditions. 

        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : momentum)
        {
            p /= alpha;
        }
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);


        // Full step in position. 
        for(int i = 0; i < latticeSize; ++i)
        {
            configuration[i] = configuration[i] + lfStepSize * momentum[i];
        }

        // half step in momentum.
        // First lattice site has special neighbour conditions. 
        momentum[0] = momentum[0] - ((mass/latticeSpacing) * (2 * configuration[0] - configuration[1] - configuration[latticeSize - 1] ) + latticeSpacing * (*potential)[configuration[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= latticeSize - 2; ++i)
        {
            momentum[i] = momentum[i] - ((mass/latticeSpacing) * (2 * configuration[i] - configuration[i+1] - configuration[i-1] ) + latticeSpacing * (*potential)[configuration[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        momentum[latticeSize-1] = momentum[latticeSize-1] - ((mass/latticeSpacing) * (2 * configuration[latticeSize-1] - configuration[0] - configuration[latticeSize-2] ) + latticeSpacing * (*potential)[configuration[latticeSize-1]]) * (lfStepSize/2.0);

        for(auto& p : momentum)
        {
            p /= alpha;
        }


    }

}

double correlationFunction(const std::vector<double> &configuration, int t)
{ 
    
    double sum = 0;
    double normalisation = 0;
    for(int i = 0; i < configuration.size(); ++i)
    {
        sum           += configuration[i]*configuration[(t+i)%(configuration.size())];
        normalisation += configuration[i]*configuration[i];
    }

    return sum/normalisation;
}

double slope(const std::vector<double>& x, const std::vector<double>& y)
{
    const auto n    = x.size();
    const auto xSum = std::accumulate(x.begin(), x.end(), 0.0);
    const auto ySum = std::accumulate(y.begin(), y.end(), 0.0);
    const auto xxSum= std::inner_product(x.begin(),x.end(),x.begin(),0.0);
    const auto xySum= std::inner_product(x.begin(),x.end(),y.begin(),0.0);
    const auto slope= (n * xySum - xSum * ySum ) / (n * xxSum - xSum * xSum);
    return slope;
}

