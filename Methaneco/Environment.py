#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
This code is under MIT license. See the License.txt file.
Environmental constants

Boris Sauterey
boris.sauterey@ens.fr
"""

# The QX are obtained by multiplying the piston velocity given in Kharecha's paper (cm.s-1) by the depth of the mixing 
# column (assumed of ~100 meters, similar to modern oceans), by converting it in d-1, the time unit of the model, and in 
# dm, to check the concentration in mol.L-1.    
QC   = 4.8E-3 * 86400 * 1e-3 * 1e2
QH   = 1.3e-2 * 86400 * 1e-3 * 1e2
QG   = 4.5E-3 * 86400 * 1e-3 * 1e2 
QN   = 1E-6   * 1e3

# The Xinf for H2, CO2, and CH4 are obtained by multiplying the partial atmospheric pressure in X (assumed fixed for now, 
# and for an atmospheric pressure of 1 bar) by the the solubility of the chemical element (the mixing ratios of H2, CO2, 
# CH4 are assumed to be equal to 1000, 2500, and 10 ppm respectively). It is fixed arbitrarily for NH4: 
Hinf = 1000*1e-6 * 1 * 7.8e-4
Cinf = 2500*1e-6 * 1 * 1e-3
Ginf = 10*1e-6   * 1 * 1.4e-3
Ninf = 1e-8

Sinf = [Hinf,Cinf,Ninf,Ginf]