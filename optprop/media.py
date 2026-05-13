"""
Module containing functions related to the media of propagation
"""
import numpy as np
from .grid import square_meshgrid
from .transforms import cart2pol
from numpy.random import randn
from .ft_functions import ift2

def ft_phase_screen(r0,N,delta,L0,l0):
    """
    Phase screen using the Fourier method

    Args:
        r0: Atmospheric coherence diameter [m]
        N: Number of sample points
        delta: Grid spacing [m]
        L0: Outer scale size [m]
        l0: Inner scale size [m]
    Returns:
        phz: Phase screen (real values)
    """

    del_f = 1/(N*delta)
    fx,fy = square_meshgrid(N,del_f)
    th,f = cart2pol(fx,fy)
    fm = 5.92/l0/(2*np.pi)
    f0 = 1/L0
    PSD_phi = 0.023*r0**(-5/3)*np.exp(-(f/fm)**2)\
        /(f**2+f0**2)**(11/6)
    PSD_phi[N//2,N//2] = 0
    cn = (randn(N,N) + 1j*randn(N,N)) * np.sqrt(PSD_phi)*del_f
    phz = np.real(ift2(cn,1))
    return phz
    