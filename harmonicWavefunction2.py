#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys 
import statistics as stat

def discreteWavefunction(x):
    return np.sqrt(np.sqrt(5)/(2*np.pi)) * np.exp(-0.5*np.sqrt(5)*x*x)

def continuumWavefunction(x):
    return 1/np.sqrt(np.pi)*np.exp(-x*x)

if __name__ == "__main__":

    positionData          = []
    wavefunctionData      = []
    wavefunctionErrorData = []

    filein = open(sys.argv[1],"r")
    inputData = filein.readlines()
    
    for line in inputData:
        token       = line.split(" ")
        position    = float(token[0])
        probability = float(token[1])
        error       = float(token[2].strip('\n'))

        positionData.append(position)
        wavefunctionData.append(probability)
        wavefunctionErrorData.append(error)


    

    trueData = np.arange(positionData[0],positionData[len(positionData)-1],0.01)
    
    # plot histogram with various customizations
    plt.plot(positionData,wavefunctionData,'x',color='red', markersize='4',label="Measured")
    plt.errorbar(positionData,wavefunctionData,yerr=wavefunctionErrorData,linestyle='none')
    plt.plot(trueData,continuumWavefunction(trueData),label="Continuum Theory",color='blue',linewidth='0.3')
    plt.plot(trueData,discreteWavefunction(trueData),'--',label="Discrete Theory",color='green')
    plt.title("Oscillator Wavefunction")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|\psi(x)|^2$')
    plt.legend()
   
    
    # show the plot
    plt.show()


		
