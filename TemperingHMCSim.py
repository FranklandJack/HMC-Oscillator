#!/usr/bin/env python3

import matplotlib.pyplot as plt

def grad_U(l,f,x):
    return 4 * l * x*x*x - 4 * l * f * x

def grad_V(u,l,x):
    return u * x + 4 * l * x * x * x

def hamiltonian(x,p,l,f,m):
    return p*p/2*m + l * (x*x-f)*(x*x-f)

if __name__ == '__main__':

    positionData = []
    momentumData = []
    hamiltonians = []
    
    iPosition = 1
    iMomentum = 1

    lfSteps    = 200
    lfStepSize = 0.1
    alpha      = 1.01

    u          = 1
    l          = 1
    f          = 1
    mass       = 1

    position = iPosition
    momentum  = iMomentum

    for x in range(0,int(lfSteps/2)):

        momentum = momentum * alpha

        momentum_halfStep = momentum - lfStepSize/2.0 * grad_V(l, f, position)

        position_fullStep = position + lfStepSize * momentum_halfStep/mass

        momenum_fullStep  = momentum_halfStep - lfStepSize/2.0 * grad_V(l, f, position_fullStep)

        momentum = momenum_fullStep

        momentum = momentum * alpha

        position = position_fullStep

        positionData.append(position)
        momentumData.append(momentum)
        hamiltonians.append(hamiltonian(position,momentum,l,f,mass))

    if 0 != lfSteps % 2:

        momentum = momentum * alpha

        momentum_halfStep = momentum - lfStepSize/2.0 * grad_V(l, f, position)

        position_fullStep = position + lfStepSize * momentum_halfStep/mass

        momenum_fullStep  = momentum_halfStep - lfStepSize/2.0 * grad_V(l, f, position_fullStep)

        momentum = momenum_fullStep

        momentum = momentum / alpha

        position = position_fullStep

        positionData.append(position)
        momentumData.append(momentum)
        hamiltonians.append(hamiltonian(position,momentum,l,f,mass))


    for x in range(0,int(lfSteps/2)):

        momentum = momentum/alpha

        momentum_halfStep = momentum - lfStepSize/2.0 * grad_V(l, f, position)

        position_fullStep = position + lfStepSize * momentum_halfStep/mass

        momenum_fullStep  = momentum_halfStep - lfStepSize/2.0 * grad_V(l, f, position_fullStep)

        momentum = momenum_fullStep

        momentum = momentum/alpha

        position = position_fullStep

        positionData.append(position)
        momentumData.append(momentum)
        hamiltonians.append(hamiltonian(position,momentum,l,f,mass))


    

    plt.figure(1)
    plt.subplot(211)

    stepNumbers = range(0,lfSteps)
    plt.plot(stepNumbers, hamiltonians, 'x', color = 'black', markersize='4')
    plt.xlabel('Leapfrog Step Number')
    plt.ylabel('Value of Hamiltonian')


    plt.subplot(212)

    plt.plot(positionData, momentumData, 'x', color='black', markersize='4')
    plt.xlabel('Position')
    plt.ylabel('Momentum')


    plt.show()

