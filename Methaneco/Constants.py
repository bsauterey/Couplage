#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
This code is under MIT license. See the License.txt file.
Environmental constants

Boris Sauterey
boris.sauterey@ens.fr
"""
from math import pi
import numpy as np

#Environmental constants
T          = 298
Av         = 6e23
p          = 1e5      # Pa     
R          = 8.31     # J.(K.mol)-1
Toc        = 283      # K
S          = 5.1e18   # cm2
g          = 9.81     #
Mair       = 28.16e-3 # kg.mol-1
ntot       = (p*S)/(g*Mair)
QH         = 1.1e-1   * 1e1          # m(x100).d-1       Kharecha
QC         = 4.1e-2   * 1e1          # m(x100).d-1       Kharecha
QN         = 1E-6     * 1e1          # m(x100).d-1       arbitrary
QG         = 3.9e-2   * 1e1          # m(x100).d-1       Kharecha

#Photolysis of CH4
a          = 1.438e15 #yr-1
b          = 2.0291

#CO conversion
Vco        = 0.030 #cm.yr-1 

#Volcanism
#Volc       = 2.6e-7 #pH2 = 200 ppm
#Volc       = 6.5e-7 #pH2 = 500 ppm
Volc       = 1.05e-6 #pH2 = 800 ppm

## Trait
rc         = 1e-6                     # µm
Vc         = (4/3)*pi*rc**3           # µm3
Qc         = (18E-12*Vc**(0.94))/10   # molX.Cell-1       Menden-Deuer and Lessard 2000
ks         = 1e-12                    # molX.L-1          arbitrary
qmax       = 1e-1                     # (d.(molX.L-1))-1  Gonzalez Cabaleiro 2015 PLOS; Kral et al. 1998
qmax       = qmax*Qc/Vc               # (d.Cell)-1
mg         = 4500                     # J.(molX.h-1)      Gonzalez Cabaleiro 2015 PLOS
mg         = 4500*24*Qc               # J.(Cell.d-1) 
kd         = 1                        # d-1               Batstone et al 2002 in GC 2015 ISME
mort       = 0.1                      # d-1               arbitrary 
thresh     = 10*Qc                    # molX.Cell-1       arbitrary
slope      = 10                       #                   arbitrary 
#gmax       = 1                        # d-1               arbitrary
gmax       = 1                        # d-1               arbitrary