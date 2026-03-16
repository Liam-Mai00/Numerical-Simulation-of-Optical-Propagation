"""
aperture_functions.py is a module containing common (aperture) functions found in the appendix of J. Schmidt's book.
Reference:
 - Numerical Simulation of Optical Propagation - J. Schmidt (Appendix B MATLAB Code Listings)
"""
import numpy as np
from scipy.special import jv

def rect(x,D=1):
    """
    Creates a rectangular function
    
    Args:
        x: Input array
    Returns:
        y: Rectangular function with width D
    """
    x = np.abs(x)
    y = (x < D/2).astype(float)
    y[x == D/2] = 0.5
    return y

def tri(t):
    """
    Creates a triangular function

    Args:
        t: Input array
    Returns:
        y: Triangular function output array
    """
    t = np.abs(t)
    y = np.zeros_like(t)
    y[t < 1.0] = 1.0 - t[t < 1.0]
    return y

def circ(x,y,D):
    """
    Creates a circular aperture.

    Args:
        x: Meshgrid of x values
        y: Meshgrid of y values
        D: Circle diameter
    Returns:
        z: Array of aperture with diameter D
    """
    r = np.sqrt(x**2 + y**2)
    z = (r < D/2).astype(float)
    z[r == D/2] = 0.5
    return z

def jinc(x):
    y = np.ones_like(x)
    idx = x != 0
    y[idx] = 2.0*jv(1,(np.pi*x[idx]))/(np.pi*x[idx])
    return y