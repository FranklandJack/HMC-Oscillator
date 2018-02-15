#!/usr/bin/env python2
from argparse import ArgumentParser
import numpy

# Quick program to calculate energy gap of QHM on discrete time lattice.

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

# Use the fact discrete lattice is an altered qho.
omega = numpy.sqrt(muSquared/m)*numpy.sqrt(1+0.25*a*a*muSquared/m)

# Then normal formula for energy gap. E_n = (n+1/2)* omega = > deltaE = omega.

print(omega)



