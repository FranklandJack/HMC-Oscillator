#include "Histogram.hpp"
Histogram::Histogram(double min, double max, int numberBins):m_min(min),m_max(max),m_binCount(numberBins)
{
	// Divide the range into bins to get the bin width.
    m_binWidth = (m_max - m_min) / m_binCount;

    // Allocate the array of bins.
    m_counts   = new double [numberBins];

    // Initialise each bin to have have zero samples in it.
    for(int bin = 0; bin < numberBins; ++bin)
    {
        m_counts[bin] = 0;
    }
    // Initially total number of samples is also zero.
    m_totalCounts = 0;
}
Histogram::~Histogram()
{
	// Dynamic memory allocation requires explicit destructor.
    delete [] m_counts;
}

int Histogram::getBins() const
{
    return m_binCount;
}

double Histogram::count(int bin) const
{
    return m_counts[bin];
}


void Histogram::operator()(double x)
{
	// Calculate relevant bin for this sample
    int bin = static_cast<int>((x-m_min)/m_binWidth);

    // Check that the bin actually exists for this histogram.
    if(bin>=0 && bin<m_binCount) 
    {
        ++m_counts[bin];
        ++m_totalCounts;
    }
    // All outliers are just ignored.
    
}

void Histogram::normalise()
{
	// Divide the frequency of each bin by the total number of counts and the bin width to normalise the histogram.
    for(int i = 0; i < m_binCount; ++i)
        m_counts[i]/=(static_cast<double>(m_totalCounts)*m_binWidth);  
}

std::ostream& operator<<(std::ostream& out, const Histogram &histogram)
{
	// Print minimum and maximum values on first line.
    out << histogram.m_min << ' ' << histogram.m_max << '\n';

    // Print the bin width on the second line.
    out << histogram.m_binWidth << '\n';

    // On remaining lines print the frequency of each bin per line.
    for(int i = 0; i < histogram.m_binCount; ++i)
        out << histogram.m_counts[i] << '\n';
    return out;
}


