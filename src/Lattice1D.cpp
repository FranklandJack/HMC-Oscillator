#include "Lattice1D.hpp"

Lattice1D::Lattice1D(int size, double spacing, double mass, Ipotential *potential) : 
																		 m_size(size), 
																		 m_spacing(spacing), 
																		 m_mass(mass),
																		 m_data(size,0.0),
																		 m_potential(potential){}

Lattice1D::Lattice1D(int size, double spacing, double mass, Ipotential *potential, 
					std::default_random_engine &generator, double min, double max) :
																		 m_size(size), 
																		 m_spacing(spacing), 
																		 m_mass(mass),  
																		 m_potential(potential)
{
	std::uniform_real_distribution<double> initialDistribution(min, max);
	m_data.reserve(m_size);
	for(int i = 0; i < m_size; ++i)
	{
		m_data.push_back(initialDistribution(generator));
	}

}

void Lattice1D::initialise(std::default_random_engine &generator, double min, double max)
{
	std::uniform_real_distribution<double> initialDistribution(min, max);
	for(auto& x : m_data)
	{
		x = initialDistribution(generator);
	}
}

int Lattice1D::getSize() const
{
	return m_size;
}

double Lattice1D::getSpacing() const
{
	return m_spacing;
}

std::vector<double> Lattice1D::getData() const
{
	return m_data;
}

double Lattice1D::getMass() const
{
	return m_mass;
}

double& Lattice1D::operator[](int index)
{
	return m_data[(index+m_size)%m_size];
}

const double& Lattice1D::operator[](int index) const
{
	return m_data[(index+m_size)%m_size];
}

std::ostream& operator<<(std::ostream& out, const Lattice1D &lattice)
{
	for(const auto& x : lattice.m_data)
	{
		out << x << ' ';
	}
	return out;
}

double Lattice1D::action() const
{
	double sum = 0.0;

	for(int index = 0; index < m_size; ++index)
	{
		sum += (0.5 * m_mass * (m_data[(index+1)%m_size]-m_data[index]) * (m_data[(index+1)%m_size]-m_data[index])) / m_spacing + m_spacing * (*m_potential)(m_data[index]);
	}

    return sum;
}

double Lattice1D::correlation(int t) const
{
    double sum = 0;
    double normalisation = 0;
    for(int i = 0; i < m_size; ++i)
    {
        sum           += m_data[i]*m_data[(t+i)%(m_size)];
        normalisation += m_data[i]*m_data[i];
    }

    return sum/normalisation;
}

std::vector<double> Lattice1D::correlation(int t1, int t2) const
{
	std::vector<double> correlationData;
	correlationData.reserve(m_size);
	for(int t = t1; t <= t2; ++t1)
	{
		correlationData.push_back(correlation(t));
	}

	return correlationData;
}

double Lattice1D::meanX() const
{
	double sum = 0;
	for(const auto& x : m_data)
	{
		sum += x;
	}

	return sum/m_size;
}

double Lattice1D::meanXSquared() const
{
	double sum = 0;
	for(const auto& x : m_data)
	{
		sum += x*x;
	}

	return sum/m_size;

}

double Lattice1D::meanXFourth() const
{
	double sum = 0;
	for(const auto& x : m_data)
	{
		sum += x*x*x*x;
	}

	return sum/m_size;
}

