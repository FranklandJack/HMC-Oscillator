#!/usr/bin/env python3
import numpy as np
import sys

def expectationXSquared(mu,m,a,N):

	R = 1 + a*a*mu*mu/(2.0*m) - a*mu*np.sqrt(1+a*a*mu*mu/(4*m))/(np.sqrt(m))
	return 1/(2*mu*np.sqrt(m+a*a*mu*mu/(4))) * (1 + np.power(R,N))/(1 - np.power(R,N))


if __name__ == "__main__":
	muSquared = float(sys.argv[1])
	m  = float(sys.argv[2])
	a = float(sys.argv[3])
	N = float(sys.argv[4])
	print(expectationXSquared(np.sqrt(muSquared),m,a,N))
