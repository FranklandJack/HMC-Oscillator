#!/usr/bin/env python3
import numpy as np
import sys

def expectationXSquared(mu,a,N):

	R = 1 + a*a*mu*mu/2.0 - a*muSquared*np.sqrt(1+a*a*mu*mu/4)
	return 1/(2*mu*np.sqrt(1+a*a*mu*mu/4)) * (1 + np.power(R,N))/(1 - np.power(R,N))


if __name__ == '__main__':
	a = float(sys.argv[1])
	muSquared = float(sys.argv[2])
	N = float(sys.argv[3])
	print(expectationXSquared(np.sqrt(muSquared),a,N))
