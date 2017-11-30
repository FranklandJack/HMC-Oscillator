#!/usr/bin/env python2.7

# import matplotlib for plotting input
import matplotlib.pyplot as plt
import numpy as np
# import sys so that we may access command line arguments
import sys 

def wavefunction(x):
    return 0.59 * np.exp(-1.1*x*x)

# define main method
if __name__ == "__main__":


    # prepare list to hold probabilites
    probabilty = []

    # prepare list to hold probabilty squared
    probabiltySquared = []

    # prepare list to hold error in probabilty
    error = []

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
    position    = np.arange(minimum+binwidth/2, maximum, binwidth)

    for i in range (0,len(position)):
        probabilty.append(0)
        probabiltySquared.append(0)
        error.append(0)

    
    # rest of data will be specific to each file
    for i in range (0,samples):
         # get file from the first command line argument
        filein = open(sys.argv[i+1],"r")

        # read lines into buffer
        inputData = filein.readlines()

        # prepare list to hold data
        yData = []
        xData = []
    
        # loop through remainder of file getting data from each line
        for line in range(2,len(inputData)):

            # add data to list defined above
            value = float(inputData[line].rstrip('\n'))
            probabilty[line-2] += value/samples
            probabiltySquared[line-2] += value*value/samples


        # close file for saftey
        filein.close()

    for i in range (0, len(probabilty)):
        error[i] = np.sqrt((probabiltySquared[i] - probabilty[i]*probabilty[i])/samples)




    
       
    
    
    print(error)
    trueData = np.arange(minimum,maximum,0.01)
    
    # plot histogram with various customizations
    plt.plot(position,probabilty,'x',color='red', markersize='4',label="Measured")
    plt.errorbar(position,probabilty,yerr=error,linestyle='none')
    plt.plot(trueData,wavefunction(trueData),label="Theory",color='blue',linewidth='0.3')
    plt.title("Harmonic Oscillator Wavefunction")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|\psi(x)|^2$')
    plt.legend()
   
    # save the plot with file name the same as its title
    #plt.savefig(plotTitle)
    
    # show the plot
    plt.show()


		
