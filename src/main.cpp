#include <iostream>//for any standard IO
#include <ctime> //for time()
#include <cstdlib> //for rand() and srand() and exit()
#include <fstream>//for file IO
#include <cmath> //will be used generally
#include <random> //used for random number generation
#include <chrono> //for system clock
#include <algorithm> //for max/min
#include <string> //for naming files
#include <vector> //for holding any 1-D configurations
#include <array>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "Ipotential.hpp"
#include "HarmonicPotential.hpp"
#include "AnharmonicPotential.hpp"
#include "OcticPotential.hpp"
#include "LatticeFunctions.hpp"
#include "Histogram.hpp"
#include "ProgressBar.hpp"
#include "Timer.hpp"
#include "makeDirectory.hpp"
#include "HMCLattice1D.hpp"
#include "DataArray.hpp"
#include "HMCInput.hpp"
#include "HMCOutput.hpp"
#include "metropolisUpdate.hpp"

// Lots of use of the standard library so use namespace here.
using namespace std;

// Simplifies the input options.
namespace po = boost::program_options;



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
    po::options_description desc("Options for hmc oscillator program");

    // Add all optional command line arguments.
    desc.add_options()
        
        ("lattice-size,L", po::value<int>(&latticeSize)->default_value(1000), "The number of lattice sites")
        ("lattice-spacing,a", po::value<double>(&latticeSpacing)->default_value(1.0), "The spacing between lattice sites")
        ("mass,m", po::value<double>(&mass)->default_value(1.0), "The mass of the oscillator")
        ("mu-squared,u", po::value<double>(&muSquared)->default_value(1), "The mu^2 of the oscillator")
        ("lambda,l", po::value<double>(&lambda)->default_value(0.0), "The lambda value of the oscillator")
        ("f-squared,f", po::value<double>(&fSquared)->default_value(0.0), "The f^2 value of the oscillator")
        ("lf-step-count,N", po::value<int>(&lfStepCount)->default_value(5), "The number of leapfrog steps")
        ("lf-step-size,d", po::value<double>(&lfStepSize)->default_value(0.2), "The leapfrog step size")
        ("configuration-count,c", po::value<int>(&configCount)->default_value(100000), "The number of configurations")
        ("burn-period,b", po::value<int>(&burnPeriod)->default_value(1000), "The burn period")
        ("measurement-interval,i", po::value<int>(&mInterval)->default_value(4), "The number of steps between measurements")
        ("number-bins,B", po::value<int>(&numBins)->default_value(100), "The number of bins in the position histogram")
        ("histogram-max-value,R", po::value<double>(&histMaxValue)->default_value(4.0), "Maximum value in histogram range")
        ("histogram-min-value,r", po::value<double>(&histMinValue)->default_value(-4.0), "Minimum value in histogram range")
        ("tempering-parameter,T", po::value<double>(&temperingParameter)->default_value(1.0), "sqrt(alpha) that is the tempering parameter")
        ("correlation-range", po::value<int>(&correlationRange)->default_value(10), "range to calculate the correlation function for.")
        ("anharmonic", "use alternative potential")
        ("octic", "use octic potential")
        ("help,h", "produce help message");

    // Make arguments available to program
    po::variables_map vm;
    po::store(po::parse_command_line(argc,argv,desc), vm);
    po::notify(vm);

    // If the user asks for help display it then exit.
    if(vm.count("help"))
    {
        cout << desc << "\n";
        return 1;
    }

    // If the user specifies alternate potential need to let the program know.
    if(vm.count("anharmonic")) 
    {
        potentialChoice = HMCInput::Potential_Anharmonic;
    }
    else if(vm.count("octic"))
    {
        potentialChoice = HMCInput::Potential_Octic;
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

    cout << inputParameters << '\n';
    int outputColumnWidth = 10;

    /*************************************************************************************************************************
    ****************************************************** Output Set Up *****************************************************
    **************************************************************************************************************************/

	// Create string which holds unique time/date stamp.
	string outputName(makeDirectory());

	// Create output file to hold the input parameters.
	ofstream inputParametersOutput(outputName + "/input.txt");

	// Create output file to hold numerical values calculated during the simulation.
	ofstream resultsOutput(outputName + "/results.txt");

	// Create output file to hold the wave function.
    ofstream wavefunctionOutput(outputName+"/wavefunction.dat");

    // Create output file to hold the correlation function data.
    ofstream correlationOutput(outputName+"/correlation.dat");

    // Create output file to hold the final configuration so it can be resused in future simulations.
    ofstream finalConfigOutput(outputName+"/finalConfiguration.dat");

    // Create output file for the mean position on each configuration.
    ofstream positionOutput(outputName+"/position.dat");

    
 
    // Create potentials for each type.
    HarmonicPotential       harmonicPotential(muSquared, lambda);
    AnharmonicPotential     anharmonicPotential(lambda, fSquared);
    OcticPotential          octicPotential(lambda, fSquared);
    Ipotential              *potential = nullptr;

    switch(potentialChoice)
    {
        case HMCInput::Potential_Harmonic: 
            potential = &harmonicPotential;
            break;

        case HMCInput::Potential_Anharmonic: 
            potential = &anharmonicPotential;
            break;

        case HMCInput::Potential_Octic: 
            potential = &octicPotential;
            break;

        default:
            cout << "No potential selected, exiting program";
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
		

    
    
    vector<double> correlation(correlationRange,0);
    vector<double> correlationSquared(correlationRange,0);


    // Set up arrays to hold the wavefunction calculated on each measured configuration.
    vector<double> wavefunction(numBins,0.0);
    vector<double> wavefunctionSquared(numBins,0.0);
    vector<double> wavefunctionError(numBins,0.0);
    



    /*************************************************************************************************************************
    *************************************************** Prepare PRN generation ***********************************************
    **************************************************************************************************************************/

    unsigned int seed = static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count());
    default_random_engine generator(seed);

    
    normal_distribution<double>       momentumDistribution(0.0, 1.0);
    uniform_real_distribution<double> positionDistribution(-1.0, 1.0);
    uniform_real_distribution<double> acceptanceDistribution(0.0, 1.0);
	

    /*************************************************************************************************************************
    *************************************************** Prepare Lattice ******************************************************
    **************************************************************************************************************************/

    // Create all vectors before iterations begin so they are not created on each iteration.
    
    vector<double> configuration;
    vector<double> originalConfiguration;
    vector<double> originalMomentum;
    vector<double> momentum;
    
    for (int i = 0; i < latticeSize; ++i)
    {
        configuration.push_back(positionDistribution(generator));
        momentum.push_back(momentumDistribution(generator));
    }
	

	HMCLattice1D lattice(latticeSize, latticeSpacing, mass, potential);


    /*************************************************************************************************************************
    *************************************************** Do HMC **************************************************************
    **************************************************************************************************************************/

    // Progress bar to inform use how far through simulation they are.
    
    ProgressBar progressBar(configCount+burnPeriod,0.5, 2, 0.0);
    for(int config = 0; config < configCount+burnPeriod; ++config)
    {
        cout << progressBar;
        progressBar.increment();
        
        // Randomize the momenta on the sites.
        for(auto& p : momentum)
        {
            p = momentumDistribution(generator);
        }

        // Save the original momenta and configurations since the leapFrog() function actually updates these vectors.
        originalConfiguration = configuration;
        originalMomentum      = momentum;
        
        // Do the leap frog update.
        leapFrogTempering(configuration, momentum, latticeSpacing, lfStepCount, lfStepSize, potential, mass, temperingParameter);
        
        // Calculate the Hamiltonian for the whole configuration before and after the update.
        double hamiltonianBefore = oscillatorHamiltonian(originalMomentum, originalConfiguration, latticeSpacing, mass, potential);
        double hamiltonianAfter = oscillatorHamiltonian(momentum, configuration, latticeSpacing, mass, potential);

        // Calculate the change in the Hamiltonian over the update.
        double deltaH = hamiltonianAfter - hamiltonianBefore;
        
        // Do the metropolis update - if its sucessful and we are out of the burn period record it.
        if(!metropolisUpdate(hamiltonianBefore, hamiltonianAfter, generator))
        {
            configuration = originalConfiguration;
            momentum      = originalMomentum;
        }

        // Otherwise record that the state has been accepted, if we have got past the burn period.
        else if(config >= burnPeriod)
        {
            ++acceptance;
        }
		
        /*
		if(lattice.leapFrogUpdate(generator, lfStepCount, lfStepSize, temperingParameter));
		{
			++acceptance;
		}
        */
        // Measurements are made at the interval specified by the user once we have exceeded the burn period. 
        if(0 == config%mInterval && config >= burnPeriod)
        {   
        	//positionData.push_back(lattice.meanX());
        	//positionSquaredData.push_back(lattice.meanX());
        	
            // Quantities of interest.
            double meanX         = 0.0;
            double meanXSquared  = 0.0;
            double meanXFourth   = 0.0;

            double meanKE        = 0.0;
            double meanPE        = 0.0;

            double meanDeltaH    = 0.0;
            double meanExpdeltaH = 0.0;

            double meanGSEnergy  = 0.0;

            // Create histogram to record the wavefunction for this configuration.
            Histogram positionHistogram(histMinValue, histMaxValue, numBins);


            for(const auto& x : configuration)
            {
                meanX        += x;
                meanXSquared += x*x;
                meanXFourth  += x*x*x*x;
                positionHistogram(x);

            }

            // Normalise the wavefunction and store it in the vector of wavefunctions.
            positionHistogram.normalise();

            for(int i = 0; i < wavefunction.size(); ++i)
            {
                double probabilty = positionHistogram.count(i);
                wavefunction[i] += probabilty/mCount;
                wavefunctionSquared[i] += probabilty * probabilty/mCount;
            }

            meanX        /= latticeSize;
            positionData.push_back(meanX);

            meanXSquared /= latticeSize;
            positionSquaredData.push_back(meanXSquared);

            meanXFourth  /= latticeSize;
            positionFourthData.push_back(meanXFourth);

            // Quicker to output position values to file on each iteration rather than storing them and outputting them later.s
            positionOutput << (config-burnPeriod)/mInterval << ' ' << meanX << '\n';

            meanPE                    = latticeAction(configuration, latticeSpacing, mass, potential)/latticeSize;            
			actionData.push_back(meanPE);

            meanKE                    = kineticEnergy(momentum)/latticeSize;
            keData.push_back(meanKE);

            meanDeltaH                = deltaH;
            dhData.push_back(deltaH);

            meanExpdeltaH             = exp(-deltaH);
            expdhData.push_back(meanExpdeltaH);

            switch(potentialChoice)
            {
            case HMCInput::Potential_Harmonic: 
                meanGSEnergy = muSquared * meanXSquared + 3.0 * lambda * meanXFourth;
                break;

            case HMCInput::Potential_Anharmonic: 
                meanGSEnergy = -4.0 * fSquared * lambda * meanXSquared + 3.0 * lambda * meanXFourth + lambda * fSquared * fSquared;
                break;

            case HMCInput::Potential_Octic: 
                meanGSEnergy = muSquared * meanXSquared + 3.0 * lambda * meanXFourth;
                break;
            }

            gsEnergyData.push_back(meanGSEnergy);


            // Store first n values of the correlation function on the lattice.
            for(int i = 0; i < correlation.size(); ++i)
            {   
                double correlationValue = correlationFunction(configuration,i);
                correlation[i]         += correlationValue/mCount;
                correlationSquared[i]  += correlationValue*correlationValue/mCount;

            } 

            
        }

        
        

    }
    progressBar.increment();
    cout << progressBar;
    cout << "\nDone!...\n" << endl;
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
    
    vector<double> correlationError(correlation.size(),0);
    for(int i = 1; i < correlationError.size();++i)
    {
        double varianceCorrelation    = (correlationSquared[i] - correlation[i]*correlation[i]) * mCount/(mCount-1);
        correlationError[i]           = sqrt(varianceCorrelation)/sqrt(mCount);
    }

       
    double position = positionData.mean();
    double positionError = positionData.error();

    double positionSquared = positionSquaredData.mean();
    double positionSquaredError = positionSquaredData.error();

    double positionFourth = positionFourthData.mean();
    double positionFourthError = positionFourthData.error();

    double kineticEnergy	   = keData.mean();
    double kineticEnergyError  = keData.error();

    double action 			   = actionData.mean();
    double actionError         = actionData.error();

    double dh                  = dhData.mean();
    double dhError             = dhData.error();

    double expdh               = expdhData.mean();
    double expdhError          = expdhData.error();

    double gsEnergy            = gsEnergyData.mean();
    double gsEnergyError        = gsEnergyData.error();
    /*

    double varianceX           = (averageX_Squared - averageX * averageX) * mCount/(mCount-1);;
    double sdX                 = sqrt(varianceX)/sqrt(mCount);

    double varianceXSquared    = (averageXSquared_Squared - averageXSquared * averageXSquared) * mCount/(mCount-1);
    double sdXSquared          = sqrt(varianceXSquared)/sqrt(mCount);

    double varianceXFourth     = (averageXFourth_Squared - averageXFourth * averageXFourth) * mCount/(mCount-1);
    double sdXFourth           = sqrt(averageXFourth)/sqrt(mCount);

    double varianceKE          = (averageKE_Squared - averageKE * averageKE) * mCount/(mCount-1);
    double sdKE                = sqrt(varianceKE)/sqrt(mCount);

    double variancePE          = (averagePE_Squared - averagePE * averagePE) * mCount/(mCount-1);
    double sdPE                = sqrt(variancePE)/sqrt(mCount);

    double varianceDeltaH      = (averageDeltaH_Squared - averageDeltaH*averageDeltaH) * mCount/(mCount-1);
    double sdDeltaH            = sqrt(varianceDeltaH)/sqrt(mCount);

    double varainceExpDeltaH   = (averageExpDeltaH_Squared - averageExpDeltaH * averageExpDeltaH) * mCount/(mCount-1);
    double sdExpDeltaH         = sqrt(varainceExpDeltaH)/sqrt(mCount);

    double varianceGSEnergy    = (averageGSEnergy_Squared - averageGSEnergy * averageGSEnergy) * mCount/(mCount-1);
    double sdGSEnergy          = sqrt(varianceGSEnergy)/sqrt(mCount);
	*/


    /*
    
    double averageDeltaE_Squared = 0;

    for(int i = 1; i < correlation.size(); ++i)
    {
        if(correlation[i]>=0)
        {
            double deltaE = -1/static_cast<double>(i) * log(correlation[i]);
            averageDeltaE += deltaE / (correlation.size()-1);
            averageDeltaE_Squared += deltaE*deltaE/(correlation.size()-1);
        }

    }    

    double varianceDeltaE = (averageDeltaE_Squared - averageDeltaE*averageDeltaE) * (correlation.size()/(correlation.size()-1));
    double sdDeltaE       = sqrt(varianceDeltaE)/sqrt(correlation.size());
    */

    /**********************************************************************************************************************
    ************************************** Output to command line ********************************************************
    **********************************************************************************************************************/    
    // Construct an object to hold the results.
    /*
    HMCOutput results
    {
    	acceptanceRate*100,
    	sdAcceptance*100,
    	averagePE,
    	sdPE,
    	averageKE,
    	sdKE,
    	averageDeltaH,
    	sdDeltaH,
    	averageExpDeltaH,
    	sdExpDeltaH,
    	averageX,
    	sdX,
    	averageXSquared,
    	sdXSquared,
    	averageXFourth,
    	sdXFourth,
    	averageGSEnergy,
    	sdGSEnergy
    };
    */
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

    // Output the results to the commandline 

    cout << results << '\n';

    /**********************************************************************************************************************
    ************************************************* File Output *********************************************************
    ***********************************************************************************************************************/
       
    //wavefunctionOutput << histMinValue << ' ' << histMaxValue << '\n' <<  (histMaxValue - histMinValue) / numBins << '\n';
    double histSpacing = (histMaxValue - histMinValue) / numBins;
    for(int i = 0;i < wavefunction.size();++i)
    {
        wavefunctionOutput << (histMinValue + histSpacing/2) + i * histSpacing << ' ' << wavefunction[i] << ' ' << wavefunctionError[i] << '\n';
    }
    
    for(int i = 0; i < correlation.size(); ++i)
    {
        correlationOutput << i << " " << correlation[i] << ' ' << correlationError[i] << '\n';
    }

    for(int i = 0; i <configuration.size();++i)
    {
        finalConfigOutput << i << ' ' << configuration[i] <<'\n';
    }

    // Output the input parameters to their file.
    inputParametersOutput << inputParameters;

    // Output the numerical results to the file.
    resultsOutput << results;



    /**********************************************************************************************************************
    ************************************************* End Program *********************************************************
    ***********************************************************************************************************************/    
   cout << "Simulation Complete! Results have been outputed to the directory " << outputName << '\n'; 
   cout << "Time take to execute (s):   " << timer.elapsed() << endl << endl; 
   return 0;
}