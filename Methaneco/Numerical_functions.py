"""
This code is under MIT license. See the License.txt file.
This file contains a set of functions aiming at peforming some basic numerical manipulations

Boris Sauterey
boris.sauterey@ens.fr
"""

import numpy as np

def diff(X):

    DeltaX = [x - X[i - 1] for i, x in enumerate(X)][1:]
    
    return(DeltaX)

def find(lst, a):
    vec = []
    for i, x in enumerate(lst):
        if x==a:
            vec.append(i)
    return vec