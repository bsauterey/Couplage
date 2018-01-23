#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code is under MIT license. See the License.txt file.
Module containing functions for rates computation

Boris Sauterey
boris.sauterey@ens.fr
"""

from Metabolisme.Energy import *
import numpy as np
import scipy.interpolate as intrp

def Yl(lam):
	"""
	Returns the effective stoichiometry taking into account
	times catabolic reaction has to run to fuel anabolic reaction
	order is H C N G X
	"""
	global Catabolism
	global Anabolism
	return(np.add(lam*Catabolism,Anabolism))

def LambdaCat(dgana,dgdiss,dgcat):
	"""
	returns lambdacat the number of times catabolic reaction has to run
	in order to fuel anabolism
	SHOULD RETURN SMTH WHEN DGCAT > 0
	"""
	return(-((dgana+dgdiss)/dgcat))

def Slim(vec,S): #TO OPTIMIZE TAKES TOO MUCH TIME
	"""
	Returns the concentration of the limiting substrate
	C is the vector containing the concentrations
	S contains the stoiciometric coefficients
	"""
	S2 = np.array(S)
	C2 = np.array(vec)
	stbalanced = np.abs(np.array(C2)/np.array(S2))
	return(vec[list(stbalanced).index(np.min(stbalanced))])

def Monod(qmax,k,x):
	"""
	Monod function of x with
	qmax the maximum rate (trait)
	k the half saturation constant (trait)
	"""
	return(qmax*(x/(k+x)))

def Mreq(mg,deltaGcat):
	"""
	Computes the requested maintenance rate from
	mg the maintenance free energy trait and the catabolic energy
	**REF**
	"""
	return(-mg/deltaGcat)

def QCat(deltaGcat,H,C,qmax,ks):
	"""
	Computes the rate at which catabolic reaction occurs **REF**
	"""
	if deltaGcat < 0:
		return(Monod(qmax,ks,Slim([H,C],[-1,-0.25])))
	else:
		return(0)

def QMet(deltaGcat,qmax,ks,Smetlim):
	"""
	Returns the rate at which the metabolism runs
	**REF**
	"""
	if deltaGcat < 0:
		return(Monod(qmax,ks,Smetlim))
	else:
		return(0)

def QAna(deltaGcat,deltaGana,lam,qcat,qmet,mreq,qmax,ks,Sanalim):
	"""
	Computes the rate at which anabolic reaction runs
	lam is lambdacat
	qcat is catabolic reaction rate
	qmet is overall metabolism rate
	mreq is the required rate for maintenance, obtained with Mreq function
	qmax is the trait for max rate
	ks is the trait for half saturation constant
	Sanalim is the limiting substrate concentration for anabolism
	**REF**
	"""

	global dgdiss

	if deltaGcat < 0 and lam > 0 and qcat > mreq:
#		return((1/lam)*(qmet-mreq))
		return((1/lam)*(qcat-mreq)*Monod(qmax,ks,Sanalim))/qmax
#	elif deltaGana + dgdiss < 0:
#		return(Monod(qmax,ks,Sanalim))
	else:
		return(0)

def Decay(mreq,qcat,deltaGcat,kd):
	"""
	Computes the decay rate of the bacterial computation (maybe a problem, can go over 1)
	mreq is the maintenance energy rate requirement
	qcat is catabolic reaction rate
	deltaGcat is catabolic reaction energy yield
	kd is decay rate trait
	**REF**
	"""
	if mreq > qcat and deltaGcat <= 0:
		return(kd*(mreq-qcat)/mreq)
	elif deltaGcat > 0:
		return(kd)
	else:
		return(0)

def Gamma(thresh,slope,gmax,x):
	"""
	Division intensity function
	Returns the probability for a cell to divide
	**REF**
	"""
	#c = (4*slope)/gmax
	#a = np.exp(-c*thresh) #We can bound the two so that the slope is defined adimensionnaly 
	return(gmax*(1/(1+np.exp(-slope*(np.log10(x)-np.log10(thresh))))))



