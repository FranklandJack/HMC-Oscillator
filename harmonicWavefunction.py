#!/usr/bin/env python2.7

# import matplotlib for plotting input
import matplotlib.pyplot as plt
import numpy as np
# import sys so that we may access command line arguments
import sys 

def wavefunction(x):
    return 1/np.sqrt(np.pi) * np.exp(-x*x)

# define main method
if __name__ == "__main__":

    # get file from the first command line argument
    filein = open(sys.argv[1],"r")

    # read lines into buffer
    inputData = filein.readlines()

    # prepare list to hold data
    yData = []
    xData = []
    

    # first line of file must be min and max values. 
    minmax      = inputData[0]
    minmaxtoken = minmax.split(" ")
    minimum     = float(minmaxtoken[0])
    maximum     = float(minmaxtoken[1])


    # second line is the bin width
    binwidth    = float(inputData[1])
    
    xData = np.arange(minimum+binwidth/2, maximum, binwidth)

   
    # loop through remainder of file getting data from each line
    for line in range(2,len(inputData)):

        # add data to list defined above
        value = float(inputData[line].rstrip('\n'))
        yData.append(value)
      
       
    
    # close file for saftey
    filein.close()

    trueData = np.arange(minimum,maximum,0.01)
    
    # plot histogram with various customizations
    plt.plot(xData,yData,'x',color='red', markersize='4',label="Measured")
    plt.plot(trueData,wavefunction(trueData),label="Theory",color='blue',linewidth='0.3')
    plt.title("Harmonic Oscillator Wavefunction")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|\psi(x)|^2$')
    plt.legend()
   
    # save the plot with file name the same as its title
    #plt.savefig(plotTitle)
    
    # show the plot
    plt.show()


		
