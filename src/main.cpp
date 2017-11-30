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
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "OscPotential.hpp"
#include "Potential1.hpp"
#include "Potential2.hpp"
#include "LatticeFunctions.hpp"
#include "Histogram.hpp"

// Lots of use of the standard library so use namespace here.
using namespace std;

// Simplifies the input options.
namespace po = boost::program_options;

int main(int argc, const char * argv[]) 
{
    // Begin timing.
    auto start = chrono::system_clock::now();

    
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

    // Set up optional command line argument.
    po::options_description desc("Options for hmc oscillator program");

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
        ("potential,p", "use alternative potential")
        ("help,h", "produce help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc,argv,desc), vm);
    po::notify(vm);

    // If the user asks for help display it then exit.
    if(vm.count("help"))
    {
        cout << desc << "\n";
        return 1;
    }

    /*************************************************************************************************************************
    ****************************************************** Create Output File ************************************************
    **************************************************************************************************************************/
    time_t startTime = chrono::system_clock::to_time_t(start);
    string outputName = ctime(&startTime);
    std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ' ' ? '_' : ch;});
    std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ':' ? '-' : ch;});
    outputName.erase(std::remove(outputName.begin(), outputName.end(), '\n'), outputName.end());
    boost::filesystem::path outPath = outputName;
    

    for(int i = 2; boost::filesystem::exists(outPath) && i < 10; ++i)
    {
        stringstream ss;
        ss << outPath << "(" << i << ")";
        outPath = ss.str();
    }
    boost::filesystem::create_directory(outPath);
    
    
    // Create potentials for each type.
    OscPotential oscPotential(lambda);
    Potential1   potential1(muSquared, lambda);
    Potential2   potential2(lambda, fSquared);
    OscPotential& potential = potential1;

    // If user specified alternative potential use that.
    if(vm.count("potential"))
    {
        potential = potential2;
    }
    


    /*************************************************************************************************************************
    *************************************************** Set up Measurements **************************************************
    **************************************************************************************************************************/
    int    mCount                   = configCount/mInterval;         
    int    acceptance               = 0;

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

    double averageExpDeltaH         = 0.0;
    double averageExpDeltaH_Squared = 0.0;

    int correlationRange = 20;
    int coOutputInterval = 1000;
    
    vector<double> correlation(correlationRange,0);
    vector<double> correlationSquared(correlationRange,0);

    int numBins = 40;
    Histogram positionHistogram(-2.0, 2.0, numBins);
    



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

   
    vector<double> configuration;
    vector<double> momentum;
    for (int i = 0; i < latticeSize; ++i)
    {
        configuration.push_back(positionDistribution(generator));
        momentum.push_back(momentumDistribution(generator));
    }

    /*************************************************************************************************************************
    *************************************************** Do HMC **************************************************************
    **************************************************************************************************************************/

    for(int config = 0; config < configCount+burnPeriod; ++config)
    {
        if(config>burnPeriod && 0==config%1000)
        {
            cout << '\r' << config << " configurations completed..." << flush;
        }

        for(auto& p : momentum)
        {
            p = momentumDistribution(generator);
        }

        vector<double> originalConfiguration = configuration;
        vector<double> originalMomentum      = momentum;
    
        leapFrog(configuration, momentum, latticeSpacing, lfStepCount, lfStepSize, potential, mass);
    
        double hamiltonianBefore = oscillatorHamiltonian(originalMomentum, originalConfiguration, latticeSpacing, mass, potential);
        double hamiltonianAfter = oscillatorHamiltonian(momentum, configuration, latticeSpacing, mass, potential);

    
        double deltaH = hamiltonianAfter - hamiltonianBefore;
    
        if(exp(-deltaH) < acceptanceDistribution(generator))
        {
            configuration = originalConfiguration;
            momentum      = originalMomentum;
        }

        else if(config >= burnPeriod)
        {
            ++acceptance;
        }

    
        if(0 == config%mInterval && config >= burnPeriod)
        {   
            double meanX         = 0.0;
            double meanXSquared  = 0.0;
            double meanXFourth   = 0.0;

            double meanKE        = 0.0;
            double meanPE        = 0.0;

            double meanExpdeltaH = 0.0;


            for(const auto& x : configuration)
            {
                meanX        += x;
                meanXSquared += x*x;
                meanXFourth  += x*x*x*x;
                positionHistogram(x);

            }

            meanX        /= latticeSize;
            meanXSquared /= latticeSize;
            meanXFourth  /= latticeSize;

            averageX                 += meanX;
            averageX_Squared         += meanX * meanX;

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

            meanExpdeltaH             = exp(-deltaH);
            averageExpDeltaH         += meanExpdeltaH;
            averageExpDeltaH_Squared += meanExpdeltaH*meanExpdeltaH;

            for(int i = 0; i < correlation.size(); ++i)
            {   
                double correlationValue = correlationFunction(configuration,i);
                correlation[i] += correlationValue;
                correlationSquared[i] += correlationValue*correlationValue;
            } 
            
        }
        /*
        if(config>=burnPeriod && 0==config%coOutputInterval)
        {
            for(int i = 0; i < correlation.size(); ++i)
            {   
                double correlationValue = correlationFunction(configuration,i)/100;
                correlation[i] += correlationValue;
                correlationSquared[i] += correlationValue*correlationValue;
            } 
        }
        */
        
        


    }

    cout << "\n Done!..." << endl;

    /*************************************************************************************************************************
    ***************************************** Calculate Observables **********************************************************
    **************************************************************************************************************************/
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
    averageExpDeltaH         /= mCount;
    averageExpDeltaH_Squared /= mCount;

    for(auto &x : correlation) x /= mCount;
    
    for(auto &x : correlationSquared) x /= mCount;
     
    

    double acceptanceRate = static_cast<double>(acceptance)/(configCount);

    double varianceAcceptance  = acceptanceRate - acceptanceRate * acceptanceRate;
    double sdAcceptance        = sqrt(varianceAcceptance)/sqrt(mCount);

    double varianceX           = averageX_Squared - averageX * averageX;
    double sdX                 = sqrt(varianceX)/sqrt(mCount);

    double varianceXSquared    = averageXSquared_Squared - averageXSquared * averageXSquared;
    double sdXSquared          = sqrt(varianceXSquared)/sqrt(mCount);

    double varianceXFourth     = averageXFourth_Squared - averageXFourth * averageXFourth;
    double sdXFourth           = sqrt(averageXFourth)/sqrt(mCount);

    double varianceKE          = averageKE_Squared - averageKE * averageKE;
    double sdKE                = sqrt(varianceKE)/sqrt(mCount);

    double variancePE          = averagePE_Squared - averagePE * averagePE;
    double sdPE                = sqrt(variancePE)/sqrt(mCount);

    double varainceExpDeltaH   = averageExpDeltaH_Squared - averageExpDeltaH * averageExpDeltaH;
    double sdExpDeltaH         = sqrt(varainceExpDeltaH)/sqrt(mCount);

    vector<double> sdCorrelation(correlation.size(),0);
    for(int i = 0; i < sdCorrelation.size();++i)
    {
        double varianceCorrelation = correlationSquared[i] - correlation[i]*correlation[i];
        sdCorrelation[i] = sqrt(varianceCorrelation)/sqrt(mCount);
    }

    double groundStateEnergy;
    double sdGroundStateEnergy;

    if(vm.count("potential"))
    {
        groundStateEnergy   = -2.0*lambda*fSquared*averageXSquared + 3.0 * lambda * averageXFourth;
        sdGroundStateEnergy = sqrt(4.0*lambda*lambda*fSquared*fSquared*averageXSquared + 9.0*lambda*lambda*varianceXFourth)/sqrt(mCount);  
    }

    else
    {
        groundStateEnergy   = muSquared*averageXSquared + 3.0 * lambda * averageXFourth;
        sdGroundStateEnergy = sqrt(muSquared*muSquared*varianceXSquared + 9.0*lambda*lambda*varianceXFourth)/sqrt(mCount);
    }  

    double averageDeltaE   = 0;
    double averageDeltaE_Squared = 0;

    for(int i = 0; i < correlation.size(); ++i)
    {
        double deltaE = -1/static_cast<double>(i) * log(correlation[i]);
        cout << deltaE << endl;
        averageDeltaE += deltaE /correlation.size();
        averageDeltaE_Squared += deltaE*deltaE/correlation.size();

    }    

    double sdDeltaE = averageDeltaE_Squared - averageDeltaE*averageDeltaE;




    

    /**********************************************************************************************************************
    ************************** OUTPUT TO TERMINAL *************************************************************************
    **********************************************************************************************************************/
    cout << endl;

    cout << "Input__"               << endl;

    cout << "Lattice Size              : " << latticeSize    << endl;
    cout << "Lattice Spacing           : " << latticeSpacing << endl;
    cout << "LeapFrog Stepsize         : " << lfStepSize     << endl;
    cout << "#LeapFrog steps           : " << lfStepCount    << endl;
    cout << "#Configurations           : " << configCount    << endl;
    cout << "Burn period               : " << burnPeriod     << endl;
    cout << "Measurement interval      : " << mInterval      << endl;

    cout << endl;

    cout << "mass                      : " << mass << endl;
    if(vm.count("potential"))
    {
    cout << "f^2                       : " << fSquared << endl;             
    }
    else
    {
    cout << "mu^2                      : " << muSquared << endl;
    }
    
    cout << "lambda                    : " << lambda << endl;

    cout << endl;

    cout << endl;

    cout << "Output__"              << endl;

    cout << "Acceptance Rate           : " << acceptanceRate         << " +/- " << sdAcceptance  << endl;
    cout << "<S>                       : " << averagePE              << " +/- " << sdPE          << endl;
    cout << "<T>                       : " << averageKE              << " +/- " << sdKE          << endl;
    cout << "<exp(-deltaH)>            : " << averageExpDeltaH       << " +/- " << sdExpDeltaH   << endl;
    cout << "<X>                       : " << averageX               << " +/- " << sdX           << endl;
    cout << "<X^2>                     : " << averageXSquared        << " +/- " << sdXSquared    << endl;
    cout << "E_0                       : " << groundStateEnergy      << " +/- " << sdGroundStateEnergy << endl;
    cout << "E_1 - E_0                 : " << averageDeltaE                 << " +/- " << sdDeltaE      << endl;

    /**********************************************************************************************************************
    ************************************************* File Output *********************************************************
    ***********************************************************************************************************************/

    ofstream wavefunction(outputName+"/wavefunction.dat");
    positionHistogram.normalise();
    wavefunction << positionHistogram;

    ofstream correlationOutput(outputName+"/correlation.dat");
    for(int i = 0; i < correlation.size(); ++i)
    {
        correlationOutput << i << " " << correlation[i] << ' ' << sdCorrelation[i] << '\n';
    }

    /**********************************************************************************************************************
    ************************************************* End Program *********************************************************
    ***********************************************************************************************************************/

   auto end = chrono::system_clock::now();
   auto elapsed = chrono::duration_cast<chrono::seconds>(end - start);
   cout << "Time take to execute (s):   " << elapsed.count() << endl << endl; 
   return 0;
}