"""
This module contains functions involving propagation from pupil plane to focal plane. Chapter 4.
"""
import numpy as np
from .ft_functions import ft2

def lens_against_ft(Uin, wvl, d1, f):
    """
    Computes the complex field at the observation plane given that a transparent object is placed against (before or after) the lens. (Equation 4.10)

    Args:
        Uin: Input field
        wvl: Wavelength [m]
        d1: Spacial sampling interval [m]
        f: Focal length [m]
    Returns:
        Uout: Complex field at observation plane
    """
    N = np.shape(Uin)[0]
    k = (2*np.pi)/wvl
    fX = np.arange(-N/2,N/2,1) / (N*d1) # Assuming x2 = y2
    x2,y2 = np.meshgrid(wvl * f * fX, wvl * f * fX)
    Uout = (np.exp(1j * (k/(2*f)) * (x2**2 + y2**2)) / (1j*wvl*f)) * ft2(Uin,d1)
    return Uout

def lens_in_front_ft(Uin, wvl, d1, f, d):
    """
    Computes the complex field at the observation plane given that a transparent object is placed at distance d before the lens. (Equation 4.12)

    Args:
        Uin: Input field
        wvl: Wavelength [m]
        d1: Spacial sampling interval [m]
        f: Focal length [m]
        d: Distance between the object placed before the lens and the lens [m]
    Returns:
        Uout: Complex field at observation plane
    """
    N = np.shape(Uin)[0]
    k = (2*np.pi)/wvl
    fX = np.arange(-N/2,N/2,1) / (N*d1)
    x2 = fX * wvl * f
    x2,y2 = np.meshgrid(x2, x2) # Assuming x2 = y2
    Uout = (1/(1j*wvl*f)) * np.exp(1j * (k/(2*f)) * (1-(d/f)) * (x2**2 + y2**2)) * ft2(Uin,d1)
    return Uout

def lens_behind_ft(Uin,wvl,d1,f,d):
    """
    Computes the complex field at the observation plane given that a transparent object is placed at distance d after the lens. (Equation 4.16)

    Args:
        Uin: Input field
        wvl: Wavelength [m]
        d1: Spacial sampling interval [m]
        f: Focal length [m]
        d: Distance between the object placed after the lens and the lens [m]
    Returns:
        Uout: Complex field at observation plane
    """
    N = np.shape(Uin)[0]
    k = (2*np.pi)/wvl
    fX = np.arange(-N/2,N/2,1) / (N*d1)
    x2 = fX * wvl * f
    x2,y2 = np.meshgrid(x2,x2)
    Uout = (f/d) * (1/(1j*wvl*d)) * np.exp(1j * (k/(2*d)) * (x2**2+y2**2)) * ft2(Uin,d1)
    return Uout