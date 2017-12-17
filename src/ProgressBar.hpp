#ifndef ProgressBar_hpp
#define ProgressBar_hpp
#include <string>
#include <iostream>
#include <iomanip>
class ProgressBar
{
private:
    double      m_percentagePrecision;
    double      m_percentageComplete;
    double      m_barPrecision;
    double      m_increment;
    std::string m_bar;

public:
    ProgressBar(double finalValue, double barPrecision = 0.5, double percentagePrecision = 2, double percentageComplete = 0.0);

    void increment();

    friend std::ostream& operator<<(std::ostream &out, const ProgressBar& progressBar);


};
#endif  /*ProgressBar_hpp*/