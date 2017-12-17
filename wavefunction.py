#!/usr/bin/env python3

# import matplotlib for plotting input
import matplotlib.pyplot as plt
import numpy as np
# import sys so that we may access command line arguments
import sys 

import statistics as stat


# define main method
if __name__ == "__main__":


    # prepare list to hold probabilities
    wavefunctions = []

    # get number of samples
    samples = len(sys.argv) - 1

    # get min max values and position values from first file
    filein = open(sys.argv[1],"r")

    # read lines into buffer
    inputData = filein.readlines()

    # first line of file must be min and max values. 
    minmax      = inputData[0]
    minmaxtoken = minmax.split(" ")
    minimum     = float(minmaxtoken[0])
    maximum     = float(minmaxtoken[1])

    # second line is the bin width
    binwidth    = float(inputData[1])

    # create an array of positions consisting of centre of each bin 
    position    = np.arange(minimum+binwidth/2, maximum, binwidth)
    

    for i in range(0, samples):

        # get ith input file
        filein = open(sys.argv[i+1],"r")

        # read the lines from ith input file
        inputData = filein.read().splitlines()
        inputData = inputData[2:]
        inputData = [float(i) for i in inputData]

        wavefunctions.append(inputData)


        # close file for safety
        filein.close()

    
    meanWavefunction    = list(map(stat.mean,zip(*wavefunctions)))
    wavefunctionVariance = [0] * len(wavefunctions[0])
    if samples != 1:
        wavefunctionVariance = list(map(stat.variance,zip(*wavefunctions)))

    wavefunctionError    = [np.sqrt(i/samples) for i in wavefunctionVariance]

    
    # plot histogram with various customizations
    plt.plot(position,meanWavefunction,'x',color='red', markersize='4',label="Measured")
    plt.errorbar(position,meanWavefunction,yerr=wavefunctionError,linestyle='none')
    plt.title("Oscillator Wavefunction")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|\psi(x)|^2$')
    plt.legend()
   
    
    # show the plot
    plt.show()


		
