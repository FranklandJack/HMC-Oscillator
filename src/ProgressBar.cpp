#include <ProgressBar.hpp>
#include <sstream>
#include <cmath>

ProgressBar::ProgressBar(double finalValue, double barPrecision, double precision, double percentageComplete): m_barPrecision{barPrecision}, m_percentagePrecision{precision}, m_percentageComplete{percentageComplete} 
{
	m_increment = 1.0/finalValue * 100;
	int initialNumberBars = static_cast<int>(percentageComplete*barPrecision);
	m_bar = std::string(initialNumberBars,'=');
}

void ProgressBar::increment()
{

	m_percentageComplete += m_increment;
	int currentNumberBars = static_cast<int>(m_percentageComplete*m_barPrecision);
	m_bar = std::string(currentNumberBars,'=');


}

std::ostream& operator<<(std::ostream &out, const ProgressBar& progressBar) 
{
		int    basePercentage     = static_cast<int>(progressBar.m_percentageComplete);
        int    decimalPercentage  = static_cast<int>(pow(10,progressBar.m_percentagePrecision) * (progressBar.m_percentageComplete-basePercentage));

		std::ostringstream barStringStream;
        barStringStream << std::setw(100*progressBar.m_barPrecision) << std::setfill(' ') << std::left << progressBar.m_bar;

        std::ostringstream percentageBaseStream;
        percentageBaseStream << std::setw(2) << std::setfill('0') << basePercentage; 

        std::ostringstream percentageDecimalStream;
        percentageDecimalStream << std::setw(progressBar.m_percentagePrecision) << std::setfill('0') << decimalPercentage;

        std::string percentage = percentageBaseStream.str() + '.' + percentageDecimalStream.str();
        out << '\r' << "Progress: " <<  percentage <<  "% " << '[' << barStringStream.str() << ']' << std::flush;
        return out;
}