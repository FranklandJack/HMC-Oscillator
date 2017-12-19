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
// Lots of use of the standard library so use namespace here.
using namespace std;

// Simplifies the input options.
namespace po = boost::program_options;

enum PotentialType
{
    Potential_Harmonic,
    Potential_Anharmonic,
    Potential_Octic,
    Potential_MAX_POTENTIAL

};

int main(int argc, const char * argv[]) 
{
    // Begin timing.
    auto start = chrono::system_clock::now();

    // Tell user that simulation is starting.
    cout << "Setting up Simulation...\n";
    
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
    PotentialType potentialChoice;

    // Histogram parameters.
    int    numBins;
    double histMaxValue;
    double histMinValue;

    // Tempering parameter sqrt(alpha)
    double temperingParameter;

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
        potentialChoice = Potential_Anharmonic;
    }
    else if(vm.count("octic"))
    {
        potentialChoice = Potential_Octic;
    }

    else
    {
        potentialChoice = Potential_Harmonic;
    }

    // Tell user their input values to check they are correct.
    int outputColumnWidth = 30;
    cout << "Input Parameters..." << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "Lattice Size: " << right << latticeSize << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "Lattice Spacing: " << right << latticeSpacing << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "LeapFrog Step Size: " << right << lfStepSize << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "#LeapFrog Steps " << right << lfStepCount<< '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "#Configurations: " << right << configCount << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "Burn Period: " << right << burnPeriod << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "Measurement Interval: " << right << mInterval << '\n';
    cout << setw(outputColumnWidth) << setfill(' ') << left << "Mass: " << right << mass << '\n';
    // Depending on which potential was used report correct parameters.
    switch(potentialChoice)
    {
        case Potential_Harmonic: 
            cout << setw(outputColumnWidth) << setfill(' ') << left << "V(x): " << right << "0.5.mu^2.x + lambda.x^4" << '\n';
            cout << setw(outputColumnWidth) << setfill(' ') << left << "mu^2: " << right << muSquared << '\n';
            cout << setw(outputColumnWidth) << setfill(' ') << left << "Lambda: " << right << lambda << '\n' << '\n';
            break;

        case Potential_Anharmonic: 
            cout << setw(outputColumnWidth) << setfill(' ') << left << "V(x): " << right << "lambda.(x^2-f^2)^2" << '\n';
            cout << setw(outputColumnWidth) << setfill(' ') << left << "f^2: " << right << fSquared << '\n';
            cout << setw(outputColumnWidth) << setfill(' ') << left << "Lambda: " << right << lambda << '\n' << '\n';
            break;

        case Potential_Octic: 
            cout << setw(outputColumnWidth) << setfill(' ') << left << "V(x): " << right << "lambda.(x^2-f^2)^4" << '\n';
            cout << setw(outputColumnWidth) << setfill(' ') << left << "f^2: " << right << fSquared << '\n';
            cout << setw(outputColumnWidth) << setfill(' ') << left << "Lambda: " << right << lambda << '\n' << '\n';
            break;
    }

    /*************************************************************************************************************************
    ****************************************************** Create Output Directory ************************************************
    **************************************************************************************************************************/

    // Create string from the time the program started.
    time_t startTime = chrono::system_clock::to_time_t(start);
    string outputName = ctime(&startTime);

    // Strip out and replace difficult characters.
    std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ' ' ? '_' : ch;});
    std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ':' ? '-' : ch;});
    outputName.erase(std::remove(outputName.begin(), outputName.end(), '\n'), outputName.end());

    // Create directory path from the string.
    boost::filesystem::path outPath = outputName;
    
    // If user calls program more than once a second so that directories will be overwritten append an index.
    for(int i = 2; boost::filesystem::exists(outPath) && i < 100; ++i)
    {
        stringstream ss;
        ss << outPath << "(" << i << ")";
        outPath = ss.str();
    }

    // Create the directory for output.
    boost::filesystem::create_directory(outPath);
    
    
    // Create potentials for each type.
    HarmonicPotential       harmonicPotential(muSquared, lambda);
    AnharmonicPotential     anharmonicPotential(lambda, fSquared);
    OcticPotential          octicPotential(lambda, fSquared);
    Ipotential              *potential = nullptr;

    switch(potentialChoice)
    {
        case Potential_Harmonic: 
            potential = &harmonicPotential;
            break;

        case Potential_Anharmonic: 
            potential = &anharmonicPotential;
            break;

        case Potential_Octic: 
            potential = &octicPotential;
            break;
    }


    //Create output file for position.
    ofstream positionOutput(outputName+"/position.dat");

    /*************************************************************************************************************************
    *************************************************** Set up Measurements **************************************************
    **************************************************************************************************************************/
    int    mCount                   = configCount/mInterval;         
    int    acceptance               = 0;

    double tunnelRate               = 0;

    double averageX                 = 0.0;
    double averageX_Squared         = 0.0;

    double averageXSquared          = 0.0;
    double averageXSquared_Squared  = 0.0;

    double averageXFourth           = 0.0;
    double averageXFourth_Squared   = 0.0;

    double averagePE                = 0.0;
    double averagePE_Squared        = 0.0;

    double averageKE                = 0.0;
    double averageKE_Squared        = 0.0;


    double averageDeltaH            = 0.0;
    double averageDeltaH_Squared    = 0.0;

    double averageExpDeltaH         = 0.0;
    double averageExpDeltaH_Squared = 0.0;

    double averageGSEnergy          = 0.0;
    double averageGSEnergy_Squared  = 0.0;

    int    correlationRange         = 10;
    
    
    vector<double> correlation(correlationRange,0);
    vector<double> correlationSquared(correlationRange,0);

    //Histogram positionHistogram2(histMinValue, histMaxValue, numBins);

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
        
        // Do the metropolis update - if the update fails we need to restore the original configuration.
        if(exp(-deltaH) < acceptanceDistribution(generator))
        {
            configuration = originalConfiguration;
            momentum      = originalMomentum;
        }

        // Otherwise record that the state has been accepted, if we have got past the burn period.
        else if(config >= burnPeriod)
        {
            ++acceptance;
        }

        
        // Measurements are made at the interval specified by the user once we have exceeded the burn period. 
        if(0 == config%mInterval && config >= burnPeriod)
        {   
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
            meanXSquared /= latticeSize;
            meanXFourth  /= latticeSize;

            averageX                 += meanX;
            averageX_Squared         += meanX * meanX;

            // Quicker to output position values to file on each iteration rather than storing them and outputting them later.s
            positionOutput << meanX << '\n';

            averageXSquared          += meanXSquared;
            averageXSquared_Squared  += meanXSquared * meanXSquared;

            averageXFourth           += meanXFourth;
            averageXFourth_Squared   += meanXFourth * meanXFourth;

            meanPE                    = latticePotentialEnergy(configuration, latticeSpacing, mass, potential)/latticeSize;
            averagePE                += meanPE;
            averagePE_Squared        += meanPE*meanPE;

            meanKE                    = kineticEnergy(momentum, mass)/latticeSize;
            averageKE                += meanKE;
            averageKE_Squared        += meanKE * meanKE;

            meanDeltaH                = deltaH;
            averageDeltaH            += meanDeltaH;
            averageDeltaH_Squared    += deltaH*deltaH;

            meanExpdeltaH             = exp(-deltaH);
            averageExpDeltaH         += meanExpdeltaH;
            averageExpDeltaH_Squared += meanExpdeltaH*meanExpdeltaH;

            switch(potentialChoice)
            {
            case Potential_Harmonic: 
                meanGSEnergy = muSquared * meanXSquared + 3.0 * lambda * meanXFourth;
                break;

            case Potential_Anharmonic: 
                meanGSEnergy = -4.0 * fSquared * lambda * meanXSquared + 3.0 * lambda * meanXFourth + lambda * fSquared * fSquared;
                break;

            case Potential_Octic: 
                meanGSEnergy = muSquared * meanXSquared + 3.0 * lambda * meanXFourth;
                break;
            }

            averageGSEnergy         += meanGSEnergy;
            averageGSEnergy_Squared += meanGSEnergy*meanGSEnergy;


            // Store first n values of the correlation function on the lattice.
            for(int i = 0; i < correlation.size(); ++i)
            {   
                double correlationValue = correlationFunction(configuration,i);
                correlation[i]         += correlationValue/mCount;
                correlationSquared[i]  += correlationValue*correlationValue/mCount;

            } 

            // If their is a sign switch between two adjacent lattice sites record this as a tunnel.
            for(int i = 0; i < correlation.size(); ++i)
            {
                if(configuration[i]/configuration[(i+1)%configuration.size()] < 0) tunnelRate+= 1.0/configuration.size();
            }
            
        }

        
        

    }
    progressBar.increment();
    cout << progressBar;
    cout << "\nDone!...\n" << endl;
    /*************************************************************************************************************************
    ***************************************** Calculate Observables **********************************************************
    **************************************************************************************************************************/

    
    // Average over all measurements.
    averageX                 /= mCount;
    averageX_Squared         /= mCount;
    averageXSquared          /= mCount;
    averageXSquared_Squared  /= mCount;
    averageXFourth           /= mCount;
    averageXFourth_Squared   /= mCount;
    averagePE                /= mCount;
    averagePE_Squared        /= mCount;
    averageKE                /= mCount;
    averageKE_Squared        /= mCount;
    averageDeltaH            /= mCount;
    averageDeltaH_Squared    /= mCount;
    averageExpDeltaH         /= mCount;
    averageExpDeltaH_Squared /= mCount;
    averageGSEnergy          /= mCount;
    averageGSEnergy_Squared  /= mCount;
    tunnelRate               /= mCount;

    double acceptanceRate = static_cast<double>(acceptance)/(configCount);

    // Calculate variance and standard error using the normal formulas.    
    double varianceAcceptance  = (acceptanceRate - acceptanceRate * acceptanceRate) * mCount/(mCount-1);
    double sdAcceptance        = sqrt(varianceAcceptance)/sqrt(mCount);

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

    double averageDeltaE   = 0;
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

    /**********************************************************************************************************************
    ************************** OUTPUT TO TERMINAL *************************************************************************
    **********************************************************************************************************************/
    cout << "Output...\n";

    cout << setw(outputColumnWidth) << setfill(' ') << left << "Acceptance Rate: " << right << acceptanceRate*100 << " +/- " << sdAcceptance*100  << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<S>:" << right << averagePE << " +/- " << sdPE << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<T>:" << right << averageKE << " +/- " << sdKE << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<deltaH>:" << right << averageDeltaH << " +/- " << sdDeltaH << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<exp(-deltaH)>:" << right << averageExpDeltaH << " +/- " << sdExpDeltaH << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<X>: " << right << averageX << " +/- " << sdX << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<X^2>:" << right << averageXSquared << " +/- " << sdXSquared << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "<x^4>:" << right << averageXFourth << " +/- " << sdXFourth << endl; 
    cout << setw(outputColumnWidth) << setfill(' ') << left << "E_0:" << right  << averageGSEnergy << " +/- " << sdGSEnergy << endl;
    cout << setw(outputColumnWidth) << setfill(' ') << left << "E_1 - E_0:" << right << averageDeltaE << " +/- " << sdDeltaE << endl;
    if(potentialChoice == Potential_Anharmonic || potentialChoice == Potential_Octic || lambda != 0)
    {
    cout << setw(outputColumnWidth) << setfill(' ') << left << "Tunnel Rate:" << right << tunnelRate << endl;
    }

    /**********************************************************************************************************************
    ************************************************* File Output *********************************************************
    ***********************************************************************************************************************/

    ofstream wavefunctionOutput(outputName+"/wavefunction.dat");
    
    //wavefunctionOutput << histMinValue << ' ' << histMaxValue << '\n' <<  (histMaxValue - histMinValue) / numBins << '\n';
    double histSpacing = (histMaxValue - histMinValue) / numBins;
    for(int i = 0;i < wavefunction.size();++i)
    {
        wavefunctionOutput << (histMinValue + histSpacing/2) + i * histSpacing << ' ' << wavefunction[i] << ' ' << wavefunctionError[i] << '\n';
    }
    
    
    ofstream correlationOutput(outputName+"/correlation.dat");
    for(int i = 0; i < correlation.size(); ++i)
    {
        correlationOutput << i << " " << correlation[i] << ' ' << correlationError[i] << '\n';
    }

    ofstream finalConfigOutput(outputName+"/finalConfiguration.dat");
    for(int i = 0; i <configuration.size();++i)
    {
        finalConfigOutput << i << ' ' << configuration[i] <<'\n';
    }

    /**********************************************************************************************************************
    ************************************************* End Program *********************************************************
    ***********************************************************************************************************************/    
   cout << "Simulation Complete! Results have been outputed to the directory " << outputName << '\n'; 
   auto end = chrono::system_clock::now();
   auto elapsed = chrono::duration_cast<chrono::seconds>(end - start);
   cout << "Time take to execute (s):   " << elapsed.count() << endl << endl; 
   return 0;
}