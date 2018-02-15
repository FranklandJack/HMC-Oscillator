#include <iostream>
#include <fstream>
#include <cmath> 
#include <random> 
#include <chrono> 
#include <string> 
#include <vector>
#include <algorithm>
#include <boost/program_options.hpp>
#include "Ipotential.hpp"
#include "HarmonicPotential.hpp"
#include "AnharmonicPotential.hpp"
#include "LatticeFunctions.hpp"
#include "Histogram.hpp"
#include "ProgressBar.hpp"
#include "Timer.hpp"
#include "makeDirectory.hpp"
#include "getTimeStamp.hpp"
#include "HMCLattice1D.hpp"
#include "DataArray.hpp"
#include "HMCInput.hpp"
#include "HMCOutput.hpp"
#include "metropolisUpdate.hpp"



int main(int argc, const char * argv[]) 
{
    // Start Timing.
	Timer timer;
    
/*************************************************************************************************************************
****************************************************** Input Parameters **************************************************
**************************************************************************************************************************/

    // Lattice Parameters. 
    int    latticeSize;
    double latticeSpacing;

    // Oscillator Parameters.
    double mass;
    double muSquared;
    double lambda;
    double fSquared;

    // HMC Parameters.
    int    lfStepCount;
    double lfStepSize;

    // Other Parameters. 
    int configCount;
    int burnPeriod;
    int mInterval;

    // Choice of potential.
    HMCInput::PotentialType potentialChoice;

    // Histogram parameters.
    int    numBins;
    double histMaxValue;
    double histMinValue;

    // Tempering parameter sqrt(alpha)
    double temperingParameter;

    // Correlation range calculator.
    int correlationRange;

    // Set up optional command line argument.
    boost::program_options::options_description desc("Options for hmc oscillator program");

    // Add all optional command line arguments.
    desc.add_options()
        
        ("lattice-size,L", boost::program_options::value<int>(&latticeSize)->default_value(100), "The number of lattice sites")
        ("lattice-spacing,a", boost::program_options::value<double>(&latticeSpacing)->default_value(1.0), "The spacing between lattice sites")
        ("mass,m", boost::program_options::value<double>(&mass)->default_value(1.0), "The mass of the oscillator")
        ("mu-squared,u", boost::program_options::value<double>(&muSquared)->default_value(1), "The mu^2 of the oscillator")
        ("lambda,l", boost::program_options::value<double>(&lambda)->default_value(0.0), "The lambda value of the oscillator")
        ("f-squared,f", boost::program_options::value<double>(&fSquared)->default_value(0.0), "The f^2 value of the oscillator")
        ("lf-step-count,N", boost::program_options::value<int>(&lfStepCount)->default_value(5), "The number of leapfrog steps")
        ("lf-step-size,d", boost::program_options::value<double>(&lfStepSize)->default_value(0.2), "The leapfrog step size")
        ("configuration-count,c", boost::program_options::value<int>(&configCount)->default_value(100000), "The number of configurations")
        ("burn-period,b", boost::program_options::value<int>(&burnPeriod)->default_value(10000), "The burn period")
        ("measurement-interval,i", boost::program_options::value<int>(&mInterval)->default_value(4), "The number of steps between measurements")
        ("number-bins,B", boost::program_options::value<int>(&numBins)->default_value(100), "The number of bins in the position histogram")
        ("histogram-max-value,R", boost::program_options::value<double>(&histMaxValue)->default_value(4.0), "Maximum value in histogram range")
        ("histogram-min-value,r", boost::program_options::value<double>(&histMinValue)->default_value(-4.0), "Minimum value in histogram range")
        ("tempering-parameter,T", boost::program_options::value<double>(&temperingParameter)->default_value(1.0), "sqrt(alpha) that is the tempering parameter")
        ("correlation-range", boost::program_options::value<int>(&correlationRange)->default_value(10), "range to calculate the correlation function for.")
        ("anharmonic", "use alternative potential")
        ("help,h", "produce help message");

    // Make arguments available to program
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc), vm);
    boost::program_options::notify(vm);

    // If the user asks for help display it then exit.
    if(vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    // If the user specifies alternate potential need to let the program know.
    if(vm.count("anharmonic")) 
    {
        potentialChoice = HMCInput::Potential_Anharmonic;
    }


    else
    {
        potentialChoice = HMCInput::Potential_Harmonic;
    }

    // Construct an input object and print the values to the command line.
    HMCInput inputParameters
    {
		latticeSize,
     	latticeSpacing,
     	mass,
     	muSquared,
     	lambda,
     	fSquared,
        lfStepCount,
     	lfStepSize,
	    configCount,
	    burnPeriod,
	    mInterval,
		potentialChoice,
		numBins,
     	histMaxValue,
     	histMinValue,
     	temperingParameter,
     	correlationRange
    };

    std::cout << inputParameters << '\n';
    int outputColumnWidth = 10;

/*************************************************************************************************************************
****************************************************** Output Set Up *****************************************************
**************************************************************************************************************************/


	// Create string which holds unique time/date stamp.
	std::string outputName(makeDirectory(getTimeStamp()));

	// Create output file to hold the input parameters.
	std::ofstream inputParametersOutput(outputName + "/input.txt");

	// Create output file to hold numerical values calculated during the simulation.
	std::ofstream resultsOutput(outputName + "/results.txt");

	// Create output file to hold the wave function.
    std::ofstream wavefunctionOutput(outputName+"/wavefunction.dat");

    // Create output file to hold the correlation function data.
    std::ofstream correlationOutput(outputName+"/correlation.dat");

    // Create output file to hold the final configuration so it can be resused in future simulations.
    std::ofstream finalConfigOutput(outputName+"/finalConfiguration.dat");

    // Create output file for the mean position on each configuration.
    std::ofstream positionOutput(outputName+"/position.dat");

    
 
    // Create potentials for each type.
    HarmonicPotential       harmonicPotential(muSquared, lambda);
    AnharmonicPotential     anharmonicPotential(lambda, fSquared);
    Ipotential              *potential = nullptr;

    switch(potentialChoice)
    {
        case HMCInput::Potential_Harmonic: 
            potential = &harmonicPotential;
            break;

        case HMCInput::Potential_Anharmonic: 
            potential = &anharmonicPotential;
            break;

        default:
            std::cout << "No potential selected, exiting program";
            return 1;
    }

/*************************************************************************************************************************
*************************************************** Set up Measurements **************************************************
**************************************************************************************************************************/

    int    mCount                   = configCount/mInterval;    

    int    acceptance               = 0;
    
    DataArray positionData;
    positionData.reserve(mCount);

    DataArray positionSquaredData;
    positionSquaredData.reserve(mCount);

    DataArray positionFourthData;
    positionFourthData.reserve(mCount);

    DataArray actionData;
    actionData.reserve(mCount);

    DataArray keData;
    keData.reserve(mCount);

    DataArray dhData;
    dhData.reserve(mCount);

    DataArray expdhData;
    expdhData.reserve(mCount);

    DataArray gsEnergyData;
    gsEnergyData.reserve(mCount);
		

    
    
    std::vector<double> correlation(correlationRange,0);
    std::vector<double> correlationSquared(correlationRange,0);


    // Set up arrays to hold the wavefunction calculated on each measured configuration.
    std::vector<double> wavefunction(numBins,0.0);
    std::vector<double> wavefunctionSquared(numBins,0.0);
    std::vector<double> wavefunctionError(numBins,0.0);
    



/*************************************************************************************************************************
*************************************************** Prepare PRN generation ***********************************************
**************************************************************************************************************************/

    unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
    std::default_random_engine generator(seed);

	

/*************************************************************************************************************************
*************************************************** Prepare Lattice ******************************************************
**************************************************************************************************************************/

	HMCLattice1D lattice(latticeSize, latticeSpacing, mass, potential);
    lattice.initialise(generator);
    HMCLattice1D currentLattice = lattice;


/*************************************************************************************************************************
*************************************************** Do HMC **************************************************************
**************************************************************************************************************************/

    // Progress bar to inform use how far through simulation they are.
    
    ProgressBar progressBar(configCount+burnPeriod,0.5, 2, 0.0);
    for(int config = 0; config < configCount+burnPeriod; ++config)
    {
        std::cout << progressBar;
        progressBar.increment();
        

        // Randomise the Momenta.
        lattice.randomiseMomenta(generator);

        // Calculate the Hamiltonian.
        double hamiltonianBefore = lattice.hamiltonian();

        // Save original lattice.
        currentLattice = lattice;

        
        // Do a leapfrog update on the lattice.
        lattice.leapFrog(lfStepCount, lfStepSize, temperingParameter);

        // Calculate Hamiltonian after.
        double hamiltonianAfter = lattice.hamiltonian();

        // Do a Metropolis update and if it fails make sure we restore the original lattice.
        if(!(metropolisUpdate(hamiltonianBefore, hamiltonianAfter, generator)))
        {
            lattice = currentLattice;

        }
        // Otherwise if we are out of the burn period record the acceptance.
        else if(config >= burnPeriod)
        {
            ++acceptance;
        }
        // Measurements are made at the interval specified by the user once we have exceeded the burn period. 
        if(0 == config%mInterval && config >= burnPeriod)
        {   
            // Numerical values.
            positionData.push_back(lattice.meanX());
            positionSquaredData.push_back(lattice.meanXSquared());
            positionFourthData.push_back(lattice.meanXFourth());
            actionData.push_back(lattice.action()/latticeSize);
            keData.push_back(lattice.kineticEnergy()/latticeSize);
            dhData.push_back(hamiltonianAfter - hamiltonianBefore);
            expdhData.push_back(exp(hamiltonianBefore - hamiltonianAfter));

            // Wave Function.

            // Create histogram to record the wavefunction for this configuration.
            Histogram positionHistogram(histMinValue, histMaxValue, numBins);

            for(int index = 0; index < latticeSize; ++index)
            {
                positionHistogram(lattice[index]);
            }

            // Normalise the wavefunction and store it in the vector of wavefunctions.
            positionHistogram.normalise();

            for(int i = 0; i < wavefunction.size(); ++i)
            {
                double probabilty = positionHistogram.count(i);
                wavefunction[i] += probabilty/mCount;
                wavefunctionSquared[i] += probabilty * probabilty/mCount;
            }

            // Correlation function.

            // Store first n values of the correlation function on the lattice.
            for(int i = 0; i < correlation.size(); ++i)
            {   
                double correlationValue = lattice.correlation(i);
                correlation[i]         += correlationValue/mCount;
                correlationSquared[i]  += correlationValue*correlationValue/mCount;

            } 

            // Ground state energy.
            gsEnergyData.push_back((*potential).groundStateEnergy(lattice.meanXSquared(),lattice.meanXFourth()));

        }

    }
    progressBar.increment();
    std::cout << progressBar;
    std::cout << "\nDone!...\n" << std::endl;

/*************************************************************************************************************************
***************************************** Calculate Observables **********************************************************
**************************************************************************************************************************/

    double acceptanceRate = static_cast<double>(acceptance)/(configCount);
    
    // Calculate variance and standard error using the normal formulas.    
    double varianceAcceptance  = (acceptanceRate - acceptanceRate * acceptanceRate) * mCount/(mCount-1);
    double sdAcceptance        = sqrt(varianceAcceptance)/sqrt(mCount);
	
    for(int i = 0; i < wavefunction.size(); ++i)
    {
        double varianceWavefunction =   (wavefunctionSquared[i] - wavefunction[i] * wavefunction[i]) * mCount/(mCount-1);
        wavefunctionError[i]        = sqrt(varianceWavefunction)/sqrt(mCount);
    }
    
    std::vector<double> correlationError(correlation.size(),0);
    for(int i = 1; i < correlationError.size();++i)
    {
        double varianceCorrelation    = (correlationSquared[i] - correlation[i]*correlation[i]) * mCount/(mCount-1);
        correlationError[i]           = sqrt(varianceCorrelation)/sqrt(mCount);
    }

       
    double position 			= positionData.mean();
    double positionError 		= positionData.error();

    double positionSquared 		= positionSquaredData.mean();
    double positionSquaredError = positionSquaredData.error();

    double positionFourth 	    = positionFourthData.mean();
    double positionFourthError  = positionFourthData.error();

    double kineticEnergy	    = keData.mean();
    double kineticEnergyError   = keData.error();

    double action 			    = actionData.mean();
    double actionError          = actionData.error();

    double dh                   = dhData.mean();
    double dhError              = dhData.error();

    double expdh                = expdhData.mean();
    double expdhError           = expdhData.error();

    double gsEnergy             = gsEnergyData.mean();
    double gsEnergyError        = gsEnergyData.error();


/**********************************************************************************************************************
************************************** Output to command line ********************************************************
**********************************************************************************************************************/    
    // Construct an object to hold the results.
    HMCOutput results
    {
    	acceptanceRate*100,
    	sdAcceptance*100,
    	action,
    	actionError,
    	kineticEnergy,
    	kineticEnergyError,
    	dh,
    	dhError,
    	expdh,
    	expdhError,
    	position,
    	positionError,
    	positionSquared,
    	positionSquaredError,
    	positionFourth,
    	positionFourthError,
    	gsEnergy,
    	gsEnergyError
    };

    // Output the results to the command line 

    std::cout << results << '\n';

/**********************************************************************************************************************
************************************************* File Output *********************************************************
***********************************************************************************************************************/
       

    double histSpacing = (histMaxValue - histMinValue) / numBins;
    for(int i = 0;i < wavefunction.size();++i)
    {
        wavefunctionOutput << (histMinValue + histSpacing/2) + i * histSpacing << ' ' << wavefunction[i] << ' ' << wavefunctionError[i] << '\n';
    }
    
    for(int i = 0; i < correlation.size(); ++i)
    {
        correlationOutput << i << " " << correlation[i] << ' ' << correlationError[i] << '\n';
    }
    
    // Output the input parameters to their file.
    inputParametersOutput << inputParameters;

    // Output the numerical results to the file.
    resultsOutput << results;

    // Output final configuration.
    finalConfigOutput << lattice;

    // Output position data.
    positionOutput << positionData;


/**********************************************************************************************************************
************************************************* End Program *********************************************************
***********************************************************************************************************************/    
    std::cout << "Simulation Complete! Results have been outputed to the directory " << outputName << '\n'; 
    std::cout << "Time take to execute (s):   " << timer.elapsed() << std::endl << std::endl; 

    return 0;
}