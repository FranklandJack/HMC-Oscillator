#!/usr/bin/env python2.7

# import matplotlib for plotting input
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
# import sys so that we may access command line arguments
import sys 

# define main method
if __name__ == "__main__":



	# get file from the first command line argument
    filein = open(sys.argv[1],"r")

    # read lines into buffer
    inputData = filein.readlines()

    # prepare list to hold data
    correlationData = []
    latticeTimeData = []
    errorData       = []
   
    # loop through remainder of file getting data from each line
    for line in inputData:

    	token = line.split(" ")
    	latticeTime = int(token[0])
    	latticeTimeData.append(latticeTime)
    	correlationValue = float(token[1])
    	correlationData.append(correlationValue)
        errorValue = float(token[2])
        errorData.append(errorValue)
  
    
    # close file for saftey
    filein.close()

    
    plt.figure(1)
    plt.subplot(211)
 
    plt.plot(latticeTimeData,correlationData,'-')
    plt.errorbar(latticeTimeData,correlationData,yerr=errorData,linestyle='none')
  
    plt.title("Corelation Function vs. latticeTime")
    plt.xlabel("latticeTime")
    plt.ylabel("correlation Function")
    plt.yscale("linear")
    
   
    plt.subplot(212)
    plt.plot(latticeTimeData,correlationData,'-')
    plt.errorbar(latticeTimeData,correlationData,yerr=errorData,linestyle='none')
  
    plt.title("Corelation Function vs. latticeTime")
    plt.xlabel("latticeTime")
    plt.ylabel("correlation Function")
    plt.yscale('log')


    
    # show the plot
    plt.show()