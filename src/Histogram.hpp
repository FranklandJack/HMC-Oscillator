#ifndef Histogram_hpp
#define Histogram_hpp
#include <iostream>
class Histogram
{
    private:
        double m_binWidth;
        int m_binCount;
        int m_min;
        int m_max;
        double *m_counts;
        int m_totalCounts;

    public:
        Histogram(double min, double max, int numberBins);
        ~Histogram();

        int getBins() const;
        int count(int bin) const;

        void operator()(double x);

        friend std::ostream& operator<<(std::ostream& out,const Histogram &histogram);


        void normalise();

};

#endif /* Histogram_hpp */