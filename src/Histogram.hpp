#ifndef Histogram_hpp
#define Histogram_hpp
#include <iostream>
class Histogram
{
    private:
    	// Width of each bin in the histogram.
        double m_binWidth;
        // Min value of histogram.
        double m_min;
        // Max value of histogram.
        double m_max;
        // Array to hold each bin.
        double *m_counts;
        // Total number of samples in the histogram.
        int m_totalCounts;
        // Number of bins in the histogram.
		int m_binCount;
    public:
        Histogram(double min, double max, int numberBins);
        ~Histogram();

        // Returns total number of bins.
        int getBins() const;

        // Returns number of samples in given bin.
        double count(int bin) const;

        // Operator overload to add some sample to the histogram.
        void operator()(double x);

        // Prints the bins and their corresponding frequency.
        friend std::ostream& operator<<(std::ostream& out,const Histogram &histogram);

        // Normalises the histogram so that the area of all bins sums to one. 
        void normalise();

};

#endif /* Histogram_hpp */