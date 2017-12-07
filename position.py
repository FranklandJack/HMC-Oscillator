#!/usr/bin/env python2.7

# import matplotlib for plotting input
import matplotlib.pyplot as plt
import numpy as np
# import sys so that we may access command line arguments
import sys 

# define main method
if __name__ == "__main__":


    # prepare list to hold probabilites
    position = []

    # prepare list to hold probabilty squared
    positionSquared = []

    # prepare list to hold error in probabilty
    error = []

    # get number of samples
    samples = len(sys.argv) - 1

    #read in number of configurations from first file.
    filein = open(sys.argv[1],"r")

    # read lines into buffer
    inputData = filein.readlines()

    #count lines
    for j in range(0, len(inputData)):
        position.append(0)
        positionSquared.append(0)
        error.append(0)


    
    # rest of data will be specific to each file
    for i in range(0,samples):
         # get file from the first command line argument
        filein = open(sys.argv[i+1],"r")

        # read lines into buffer
        inputData = filein.readlines()

    
        # loop through remainder of file getting data from each line
        for line in range(0,len(inputData)):

            # add data to list defined above
            value = float(inputData[line].rstrip('\n'))
            position[line] += value/samples
            positionSquared[line] += value*value/samples


        # close file for saftey
        filein.close()

    for i in range(0, len(position)):
        error[i] = np.sqrt((positionSquared[i] - position[i]*position[i])/samples)



    configurations = np.arange(0,len(position),1)
    

    
    # plot histogram with various customizations
    plt.plot(configurations,position,'-',color='red', linewidth='0.5',label="Measured")
    plt.errorbar(configurations,position,yerr=error,linestyle='none')
    plt.legend()
   
    # save the plot with file name the same as its title
    #plt.savefig(plotTitle)
    
    # show the plot
    plt.show()


        
