#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code is under MIT license. See the License.txt file.
This package contains the energetical constants **Including constants specific to H2 autotrophs metabolism : NoC and Gamma***

Boris Sauterey
boris.sauterey@ens.fr
"""

import numpy as np
from math import *

Catabolism = np.array([-1,-0.25,0,0.25,0])
Anabolism = np.array([-2.1,-1,-0.2,0,1])

R = 8.314 #Gas Constant J/K/mol
TS = 298

## Metabolic characteristics of H2 autotrophs metabolism
NoC = 1 #Carbon source chain length
gamma = 4 #Carbon oxidation state in carbon source

## Standard energy for H2 autotrophs metabolism **REFERENCE**
deltaG0Cat = -36800 #%J/molH
deltaG0Ana = -29200 #J/molC

deltaH0Cat = -58200 #J/molH
deltaH0Ana = -79900 #J/mol

def Dgdiss(NoC,gamma):
	"""
	Returns the dissipaed energy during metabolism **REFERENCE**
	"""
	return((200 + 18*(6 - NoC)**1.8 + exp((-0.2 - gamma)**2)**0.16*(3.6 + 0.4*NoC))*1000)

dgdiss = Dgdiss(NoC,gamma)

def DeltaG0(T,DeltaG0s,DeltaH0s):
    """
    Compute the DeltaG0 value when the temperature is different from the standard temperature TS = 298 K
    """
    global TS
    return(DeltaG0s*(T/TS)+DeltaH0s*((TS-T)/TS))

def DeltaG(T,DeltaG0s,C,S,DeltaH0s):
    """
    Computes the Gibbs energy of the reaction
    C is the list of reactants and products
    S is their stoichiometric coefficients
    (Must be in the same order : H,C,N,G,X!!)
    DeltaG0 and DeltaH0 are for T=298K
    """
    global R
    global TS

    DeltaG0T = DeltaG0(T,DeltaG0s,DeltaH0s)
    DeltaG = DeltaG0T + R*T*np.log(np.product(np.array(C)**np.array(S)))
    return(DeltaG)



def DeltaGcat(T,H,C,G):
	"""
	Computes the catabolic gibbs energy for the H2 autotrophs
	"""
	global deltaG0Cat
	global deltaH0Cat

	return(DeltaG(T,deltaG0Cat,[H,C,G],[-1,-0.25,0.25],deltaH0Cat))

def DeltaGana(T,H,C,N,x):
	"""
	Computes the anabolic gibbs energy for the H2 autotrophs
	"""
	global deltaG0Ana
	global deltaH0Ana

	return(DeltaG(T,deltaG0Ana,[H,C,N,x],[-2.1,-1,-0.2,1],deltaH0Ana))