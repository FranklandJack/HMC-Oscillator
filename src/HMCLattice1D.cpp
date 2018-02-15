#include "HMCLattice1D.hpp"
#include "metropolisUpdate.hpp"

HMCLattice1D::HMCLattice1D(int size, double spacing, double mass, Ipotential *potential) : Lattice1D(size, spacing, mass, potential), m_momenta(size,0.0) {}

void HMCLattice1D::randomiseMomenta(std::default_random_engine& generator)
{
	static std::normal_distribution<double> momentaDistribution(0.0, 1.0);

	for(auto& p : m_momenta)
	{
		p = momentaDistribution(generator);
	}
}

double HMCLattice1D::kineticEnergy() const
{
	double sum = 0;
	for(const auto& p : m_momenta)
	{
		sum += p*p/2.0;
	}

	return sum;
}

double HMCLattice1D::hamiltonian() const
{
	return kineticEnergy() + action();
}

void HMCLattice1D::leapFrog(int lfStepCount, double lfStepSize, double alpha)
{   
    for(int n = 0; n < lfStepCount/2; ++n)
    {

        // Intial half step in m_momenta.
        // First lattice site has special neighbour conditions. 

        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p *= alpha;
        }
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);


        // Full step in position. 
        for(int i = 0; i < m_size; ++i)
        {
            m_data[i] = m_data[i] + lfStepSize * m_momenta[i];
        }

        // half step in m_momenta.
        // First lattice site has special neighbour conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);

        for(auto& p : m_momenta)
        {
            p *= alpha;
        }


    }

    if(0 != lfStepCount%2)
    {
        // Intial half step in m_momenta.
        // First lattice site has special neighbour conditions. 

        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p *= alpha;
        }
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);


        // Full step in position. 
        for(int i = 0; i < m_size; ++i)
        {
            m_data[i] = m_data[i] + lfStepSize * m_momenta[i];
        }

        // half step in m_momenta.
        // First lattice site has special neighbour conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);

        for(auto& p : m_momenta)
        {
            p /= alpha;
        }

    }

    for(int n = 0; n < lfStepCount/2; ++n)
    {

        // Intial half step in m_momenta.
        // First lattice site has special neighbour conditions. 

        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p /= alpha;
        }
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);


        // Full step in position. 
        for(int i = 0; i < m_size; ++i)
        {
            m_data[i] = m_data[i] + lfStepSize * m_momenta[i];
        }

        // half step in m_momenta.
        // First lattice site has special neighbour conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);

        for(auto& p : m_momenta)
        {
            p /= alpha;
        }


    }

    for(auto& p : m_momenta)
    {
        p = -p;
    }
    
        
        
}

bool HMCLattice1D::leapFrogUpdate(std::default_random_engine& generator, 
								int lfStepCount, 
								double lfStepSize, 
								double alpha)
{
	std::vector<double> currentState = m_data;
	double currentHamiltonian = hamiltonian();

	leapFrog(lfStepCount, lfStepSize, alpha);

	double newHamiltonian = hamiltonian();

	if(metropolisUpdate(currentHamiltonian,newHamiltonian,generator))
	{
		return true;
	}

	else
	{
		m_data = currentState;
		return false;
	}


}


void leapFrog(const HMCLattice1D &currentLattice, HMCLattice1D &updatedLattice, int lfStepCount, double lfStepSize, double alpha = 1.0)
{
    // Define some local variables to make things cleaner.
    int latticeSize       = currentLattice.getSize();
    double latticeSpacing = currentLattice.getSpacing();

    
}