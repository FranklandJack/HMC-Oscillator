#!/usr/bin/env python2
from argparse import ArgumentParser
import numpy

# Quick program to caluclate energy gs energy of QHM on discrete time lattice.

# Parse command line arguments.
parser = ArgumentParser(description=__doc__)
parser.add_argument('--spacing', '-a', default=1, type=float, help='Lattice Spacing.')
parser.add_argument('--size', '-N', default=1000, type=int, help='Lattice Size.')
parser.add_argument('--mass', '-m', default=1, type=float,help='Mass of particle.')
parser.add_argument('--muSquared', '-u', default=1, type=float,help='mu parameter of harmonic potential.')
args = parser.parse_args()


# create parameters.
a 	= args.spacing
N    	= args.size
m    	= args.mass
muSquared   = args.muSquared

# Use virial theorem to calculate the gs energy. E_0 = mu^2 * <x^2>.

R = 1 + 0.5*a*a*muSquared/m - a*numpy.sqrt(muSquared)*numpy.sqrt(1/m + 0.25*a*a*muSquared/(m*m))

xSquared = 1/(2*numpy.sqrt(muSquared)*numpy.sqrt(m+0.25*a*a*muSquared))*((1+numpy.power(R,N))/(1-numpy.power(R,N)))

energy = muSquared*xSquared

print(energy)



