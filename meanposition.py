#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys 

if __name__ == "__main__":

    latticeTime  = []
    positionData = []



    filein = open(sys.argv[1],"r")
    inputData = filein.readlines()
    
    for line in inputData:
        token       = line.split(" ")
        
        time = float(token[0])
        position    = float(token[1])
        latticeTime.append(time)
        positionData.append(position)


    
    # plot histogram with various customizations
    plt.plot(positionData,latticeTime,'-x',linewidth=0.3)
    plt.title(r'$<x>$')
    plt.xlabel(r'$<x>$')
    plt.ylabel('configuration')
    plt.legend()
   
    
    # show the plot
    plt.show()