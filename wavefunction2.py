#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys 
import statistics as stat

if __name__ == "__main__":


    # Get file.
    filein = open(sys.argv[1],"r")

    # Read lines into buffer.
    inputData = filein.readlines()

    # First line of file must be min and max values. 
    minmax      = inputData[0]
    minmaxtoken = minmax.split(" ")
    minimum     = float(minmaxtoken[0])
    maximum     = float(minmaxtoken[1])

    # Second line is the bin width.
    binwidth    = float(inputData[1])

    # Create an array of positions consisting of centre of each bin. 
    position    = np.arange(minimum+binwidth/2, maximum, binwidth)

    wavefunction = []
    wavefunctionError = []
    
    for i in range(2,len(inputData)):
        token       = inputData[i].split(" ")
        probability = float(token[0])
        error       = float(token[1].strip('\n'))

        wavefunction.append(probability)
        wavefunctionError.append(error)

    
    # plot histogram with various customizations
    plt.plot(position,wavefunction,'x',color='red', markersize='3',label="Measured")
    plt.errorbar(position,wavefunction,yerr=wavefunctionError,linestyle='none')
    plt.title("Oscillator Wavefunction")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|\psi(x)|^2$')
    plt.legend()
   
    
    # show the plot
    plt.show()


		
