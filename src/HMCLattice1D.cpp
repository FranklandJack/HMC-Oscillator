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
    // Since there are possibly an odd number of time-steps we divide the trajectory into two halves.
    for(int n = 0; n < lfStepCount/2; ++n)
    {

        // For tempering multiply each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p *= sqrt(alpha);
        }
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */

        
        // First lattice site has special boundary conditions.
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal boundary conditions.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special boundary conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);
        


        // Full step in position. 
        for(int i = 0; i < m_size; ++i)
        {
            m_data[i] = m_data[i] + lfStepSize * m_momenta[i];
        }

        // Second half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */

        
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
        

        // For tempering we multiply the momenta by sqrt(alpha) after the second half step in momentum.
        for(auto& p : m_momenta)
        {
            p *= sqrt(alpha);
        }

        
    }

    // If there are an odd number of steps then to temper we multiply before the first momentum half step and 
    // divide after the second half step.
    if(0 != lfStepCount%2)
    {
        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p *= sqrt(alpha);
        }


        // First half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        // First half step has special boundary conditions.
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
        // Second half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        
        // Second half step in m_momenta.

        // First lattice site has special boundary conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special boundary conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);
        
        // After second half step in momentum we divide by the sqrt(alpha) to temper.
        for(auto& p : m_momenta)
        {
            p /= sqrt(alpha);
        }

    }

    for(int n = 0; n < lfStepCount/2; ++n)
    {
        // For tempering divide each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p /= sqrt(alpha);
        }

        // First half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        // First half step has special boundary conditions.
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

        // Final step in momenta.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        // First lattice site has special neighbour conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);
        

        // To temper we divide by sqrt(alpha) after second half step in momentum.
        for(auto& p : m_momenta)
        {
            p /= sqrt(alpha);
        }

    }   

    // Negate all momenta so the dynamics is reversible.
    for(auto& p : m_momenta)
    {
        p = -p;
    }      
        
}

void HMCLattice1D::leapFrog(int lfStepCount, double lfStepSize, double alpha, std::ostream& hamiltonianOutput)
{   

    // Print initial Hamiltonian.
    int stepCounter = 0;
    hamiltonianOutput << stepCounter++ << ' ' << hamiltonian() << '\n';

    // Since there are possibly an odd number of time-steps we divide the trajectory into two halves.
    for(int n = 0; n < lfStepCount/2; ++n)
    {

        // For tempering multiply each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p *= sqrt(alpha);
        }
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */

        
        // First lattice site has special boundary conditions.
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal boundary conditions.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special boundary conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);
        


        // Full step in position. 
        for(int i = 0; i < m_size; ++i)
        {
            m_data[i] = m_data[i] + lfStepSize * m_momenta[i];
        }

        // Second half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */

        
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
        

        // For tempering we multiply the momenta by sqrt(alpha) after the second half step in momentum.
        for(auto& p : m_momenta)
        {
            p *= sqrt(alpha);
        }

        hamiltonianOutput << stepCounter++ << ' ' << hamiltonian() << '\n';

        
    }

    // If there are an odd number of steps then to temper we multiply before the first momentum half step and 
    // divide after the second half step.
    if(0 != lfStepCount%2)
    {
        // For tempering multiple each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p *= sqrt(alpha);
        }


        // First half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        // First half step has special boundary conditions.
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
        // Second half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        
        // Second half step in m_momenta.

        // First lattice site has special boundary conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special boundary conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);
        
        // After second half step in momentum we divide by the sqrt(alpha) to temper.
        for(auto& p : m_momenta)
        {
            p /= sqrt(alpha);
        }

        hamiltonianOutput << stepCounter++  << ' ' << hamiltonian() << '\n';

    }

    for(int n = 0; n < lfStepCount/2; ++n)
    {
        // For tempering divide each momenta by sqrt(alpha) before first half step.
        for(auto& p : m_momenta)
        {
            p /= sqrt(alpha);
        }

        // First half step in momentum.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        // First half step has special boundary conditions.
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

        // Final step in momenta.
        /*
        // Can use overloaded access operator to insure periodic boundary conditions are applied.
        for(int i = 0; i < m_size; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * (*this)[i] - (*this)[i+1] - (*this)[i-1] ) + m_spacing * (*m_potential)[(*this)[i]]) * (lfStepSize/2.0);
        }
        */
        // First lattice site has special neighbour conditions. 
        m_momenta[0] = m_momenta[0] - ((m_mass/m_spacing) * (2 * m_data[0] - m_data[1] - m_data[m_size - 1] ) + m_spacing * (*m_potential)[m_data[0]]) * (lfStepSize/2.0);

        // All remaining sites except last have normal neighbours.
        for(int i = 1; i <= m_size - 2; ++i)
        {
            m_momenta[i] = m_momenta[i] - ((m_mass/m_spacing) * (2 * m_data[i] - m_data[i+1] - m_data[i-1] ) + m_spacing * (*m_potential)[m_data[i]]) * (lfStepSize/2.0);
        }

        // Final site has special neighbour conditions.
        m_momenta[m_size-1] = m_momenta[m_size-1] - ((m_mass/m_spacing) * (2 * m_data[m_size-1] - m_data[0] - m_data[m_size-2] ) + m_spacing * (*m_potential)[m_data[m_size-1]]) * (lfStepSize/2.0);
        

        // To temper we divide by sqrt(alpha) after second half step in momentum.
        for(auto& p : m_momenta)
        {
            p /= sqrt(alpha);
        }

        hamiltonianOutput << stepCounter++ << ' ' << hamiltonian() << '\n';

    }   

    // Negate all momenta so the dynamics is reversible.
    for(auto& p : m_momenta)
    {
        p = -p;
    }      
        
}