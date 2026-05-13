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

def ft_sh_phase_screen(r0,N,delta,L0,l0):
    """
    Function generating phase screens with random draws
    """
    D = N*delta
    phz_hi = ft_phase_screen(r0,N,delta,L0,l0)
    x,y = square_meshgrid(N,delta)
    phz_lo = np.zeros_like(phz_hi)
    for p in range(1,4):
        del_f = 1/(3**p*D)
        fx = np.arange(-1,2) * del_f
        fx,fy = np.meshgrid(fx,fx)
        th,f = cart2pol(fx,fy)
        fm = 5.92/l0/(2*np.pi)
        f0 = 1/L0
        PSD_phi = 0.023*r0**(-5/3)*np.exp(-(f/fm)**2)\
        /(f**2+f0**2)**(11/6)
        PSD_phi[1,1] = 0
        cn = (randn(3,3) + 1j*randn(3,3)) * np.sqrt(PSD_phi)*del_f
        SH = np.zeros((N,N))
        for ii in range(0,9):
            SH = SH + cn[ii] * \
                np.exp(1j*2*np.pi*(fx[ii]*x+fy[ii]*y))
        phz_lo = phz_lo + SH
    phz_lo = np.real(phz_lo) - np.mean(np.real(phz_lo))
    return phz_lo, phz_hi
    