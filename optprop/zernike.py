"""
This module contains functions required to compute the Zernike polynomials
"""

import pandas as pd
import pathlib
import os
import numpy as np
from scipy.special import gamma

package_dir = pathlib.Path(__file__).parent.resolve()
zernike_path = os.path.join(package_dir,"data","zernike_index.csv")

zernike_index = pd.read_csv(zernike_path, index_col=False, header=None)
zernike_index = zernike_index.to_numpy(dtype=int)

# def zaf(m,i,theta):
#     """
#     This function computes the azimuthal factor within the Zernike Polynomial

#     Args:
#         m: Zernike index m
#         i: Zernike index i
#         theta: Zernike polar angle
#     Returns:
#         G: Zernike azimuthal factor 
#     """

#     if i % 2 != 0:  # When i is odd
#         if isinstance(i,int) == False:
#             print("i is not integer")
#         G = np.sin(m*theta)
#     else:
#         G = np.cos(m*theta)
#     return G

def zrf(n,m,r):
    """
    This function computes the radial factor within the Zernike Polynomial

    Args:
        n: Zernike index n
        m: Zernike index m
        r: Radial coordinates [m]
    Returns:
        R: Zernike radial factor
    """

    R = 0
    for s in range(0,(n-m)//2+1):
        R_num = (-1)**s * gamma(n-s + 1)
        R_den = gamma(s+1) * gamma(((n+m)/2)-s + 1) * gamma(((n-m)/2)-s + 1)
        R += (R_num/R_den) * (r**(n-2*s))
    return R

def zernike(i,r,theta):
    """
    This function computes the Zernike circular polynomial of mode number 'i'
    
    Args:
        i: Mode number (1-based indexing)
        r: Radial coordinates (radius) [m]
        theta: Polar angle
    Returns:
        Z: Zernike circular polynomial at mode 'i'
    """

    if i <= 0:
        raise ValueError("Mode number i must be greater than 0")

    n = zernike_index[i-1,0]
    m = zernike_index[i-1,1]
    i = zernike_index[i-1,2]    # Converts i to 1-based indexing

    if i % 2 == 0:  # When i is even
        Z = np.sqrt(2*(n+1)) * zrf(n,m,r) * np.cos(m*theta)
    else:           # When i is odd
        Z = np.sqrt(2*(n+1)) * zrf(n,m,r) * np.sin(m*theta)
    return Z