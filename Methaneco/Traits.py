#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code is under MIT license. See the License.txt file.
This file contains the traits and traits functions and traits tradeoffs

Boris Sauterey
boris.sauterey@ens.fr
"""

from math import pi
import numpy as np

rc         = 1e-6                     # µm
Vc         = (4/3)*pi*rc**3           # µm3
Qc         = (18E-12*Vc**(0.94))/10   # molX.Cell-1       Menden-Deuer and Lessard 2000
ks         = 1e-12                    # molX.L-1          arbitrary
qmax       = 1e-1                     # (d.(molX.L-1))-1  Gonzalez Cabaleiro 2015 PLOS
qmax       = qmax*Qc/Vc               # (d.Cell)-1
mg         = 4500                     # J.(molX.h-1)      Gonzalez Cabaleiro 2015 PLOS
mg         = 4500*24*Qc               # J.(Cell.d-1) 
kd         = 1                        # d-1               Batstone et al 2002 in GC 2015 ISME
mort       = 0.1                      # d-1               arbitrary 
thresh     = 10*Qc                    # molX.Cell-1       arbitrary
slope      = 10                       #                   arbitrary 
gmax       = 1                        # d-1               arbitrary

starters = [rc,Vc,Qc,ks,qmax,mg,kd,thresh,slope,gmax] # a good start for traits
