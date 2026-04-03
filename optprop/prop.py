"""
This module contains functions for various propagation techniques.
"""
import numpy as np
from .ft_functions import ft2

def fraunhofer_prop(Uin,wvl,d1,Dz):
    """
    This function performs the Frounhofer propagation using the diffraction integral
    
    Args:
        Uin:    The field at the source plane (x1,y1)
        wvl:    Wavelength [m]
        d1:     Sampling interval at the source plane [m]
        Dz:     Propagation distance [m]
    Returns:
        Uout:   Complex field at the observation plane (x2,y2)
    """
    N = np.shape(Uin)[0]                    # Assume 2D square matrix
    k = (2*np.pi)/wvl                       # Wavenumber
    fX = np.arange(-N/2,N/2,1) / (N*d1)     # Frequency samples
    x2,y2 = np.meshgrid(wvl*Dz*fX,wvl*Dz*fX)
    Uout = (np.exp(1j*(k/(2*Dz))*(x2**2+y2**2))/(1j*wvl*Dz)) * ft2(Uin,d1)
    return Uout,x2,y2   # May need to reorder this.

def one_step_prop(Uin,wvl,d1,Dz):
    """
    This function performs the one-step propagation using the Fresnel integral

    Args:
        Uin: Input field for one-step propagation (x1,y1)
        wvl: Wavelength [m]
        d1: Sampling interval at source plane [m]
        Dz: Propagation distance [m]
    Returns:
        Uout: Complex field at the observation plane (x2,y2)
    """
    N = np.shape(Uin)[0]
    k = (2*np.pi)/wvl
    x1 = np.arange(-N/2,N/2,1) * d1
    x1,y1 = np.meshgrid(x1,x1) # Assume square grid
    d2 = (wvl*Dz)/(N*d1)
    x2 = np.arange(-N/2,N/2,1) * d2
    x2,y2 = np.meshgrid(x2,x2)
    Uout = (1/(1j*wvl*Dz)) * np.exp(1j*k/(2*Dz)*(x2**2 + y2**2)) * ft2(Uin*np.exp(1j*k/(2*Dz)*(x1**2+y2**2)), d1)
    return Uout,x2,y2