#include "Histogram.hpp"
Histogram::Histogram(double min, double max, int numberBins):m_min(min),m_max(max),m_binCount(numberBins)
{
    m_binWidth = (max - min) /m_binCount;
    m_counts   = new double [numberBins];
    for(int bin = 0; bin < numberBins; ++bin)
    {
        m_counts[bin] = 0;
    }
    m_totalCounts = 0;
}
Histogram::~Histogram()
{
    delete [] m_counts;
}

int Histogram::getBins() const
{
    return m_binCount;
}

int Histogram::count(int bin) const
{
    return m_counts[bin];
}

void Histogram::operator()(double x)
{
    int bin = static_cast<int>((x-m_min)/m_binWidth);
    if(bin>=0 && bin<m_binCount) 
    {
        ++m_counts[bin];
    }
    else
    {
        
    }
    ++m_totalCounts;
}

void Histogram::normalise()
{
    for(int i = 0; i < m_binCount; ++i)
        m_counts[i]/=(static_cast<double>(m_totalCounts)*m_binWidth);  
}

std::ostream& operator<<(std::ostream& out, const Histogram &histogram)
{
    out << histogram.m_min << ' ' << histogram.m_max << '\n';
    out << histogram.m_binWidth << '\n';
    for(int i = 0; i < histogram.m_binCount; ++i)
        out << histogram.m_counts[i] << '\n';
    return out;
}


