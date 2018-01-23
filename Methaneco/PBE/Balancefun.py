#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code is under MIT license. See the License.txt file.
Module coding for PBE useful functions

Antonin Affholder
antonin.affholder@ens.fr
"""


import numpy as np 
from scipy.stats import truncnorm
from math import pi
import scipy.interpolate as intrp


def p(x1,x2,sig):
	"""
	Defines the partition function
	something like the probability of a cell with internal coordinate x2 dividing into 2 cells in state x1
	"""
	mu = x2/2
	clipa = 0
	clipb = x2
	a,b = (clipa-mu)/sig,(clipb-mu)/sig
	return(truncnorm.pdf(x1,a,b,loc=mu,scale=sig))


def Pmap(X,sig):
	"""
	Returns the map of the p function, it is a long operation so better do it only once
	Optimization is not so easy
	"""
	M = np.array([p(X,x2,sig) for x2 in X]).transpose() #In that way, M[a][:] is the profile of where can the particle come from?
	return(M)


def Dotmap(Profile,Qanamap):
	"""
	Returns the X derivative term using the discretization
	Can be optimized, we should find adequate initial conitions to our problem
	"""
	L = Profile*Qanamap
	dL = [-L[0]]+[L[i-1]-L[i] for i in range(1,len(L)-1)]+[L[-2]]
	return(dL)

