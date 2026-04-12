"""
This module contains functions for various propagation techniques.
"""
import numpy as np
from .ft_functions import ft2, ift2
from .grid import square_meshgrid

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
        x2: 2D meshgrid of x coordinates at observation
        y2: 2D meshgrid of y coordinates at observation
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
        x2: 2D meshgrid of x coordinates at observation
        y2: 2D meshgrid of y coordinates at observation
    """
    N = np.shape(Uin)[0]                # Number of samples
    k = (2*np.pi)/wvl                   # Wavenumber
    x1 = np.arange(-N/2,N/2,1) * d1     # 1D x1 coordinates
    x1,y1 = np.meshgrid(x1,x1)          # Assume square grid
    d2 = (wvl*Dz)/(N*d1)                # Define x2 spacing
    x2 = np.arange(-N/2,N/2,1) * d2     # 1D x2 coordinates
    x2,y2 = np.meshgrid(x2,x2)          # Assume square grid
    Uout = (1/(1j*wvl*Dz)) * np.exp(1j*k/(2*Dz)*(x2**2+y2**2)) * ft2(Uin*np.exp(1j*k/(2*Dz)*(x1**2+y1**2)), d1)
    return Uout,x2,y2

def two_step_prop(Uin,wvl,d1,d2,Dz):
    """
    This function performs the two-step propagation using the Fresnel integral

    Args:
        Uin: Input field for two-step propagation (x1,y1)
        wvl: Wavelength [m]
        d1: Sampling interval at source plane [m]
        d2: Sampling interval at observation plane [m]
        Dz: Propagation distance [m]
    Returns:
        Uout: Complex field at the observation plane (x2,y2)
        x2: 2D meshgrid of x coordinates at observation
        y2: 2D meshgrid of y coordinates at observation
    """
    N = np.shape(Uin)[0]
    k = (2*np.pi)/wvl
    x1 = np.arange(-N/2,N/2,1) * d1
    x1,y1 = np.meshgrid(x1,x1)
    m = d2/d1
    Dz1 = Dz/(1-m)                  # Equation 6.25
    d1a = wvl*np.abs(Dz1)/(N*d1)    # Equation 6.19
    x1a = np.arange(-N/2,N/2,1) * d1a
    x1a,y1a = np.meshgrid(x1a,x1a)
    
    Uitm = (1/(1j*wvl*Dz1)) * np.exp(1j*(k/(2*Dz1))*(x1a**2+y1a**2)) * ft2(Uin*np.exp(1j*(k/(2*Dz1))*(x1**2+y1**2)),d1)
    Dz2 = Dz - Dz1
    x2 = np.arange(-N/2,N/2,1) * d2
    x2,y2 = np.meshgrid(x2,x2)
    Uout = (1/(1j*wvl*Dz2)) * np.exp(1j*(k/(2*Dz2))*(x2**2+y2**2)) * ft2(Uitm*np.exp(1j*(k/(2*Dz2))*(x1a**2+y1a**2)),d1a)
    return Uout,x2,y2

def ang_spec_prop(Uin,wvl,d1,d2,Dz):
    """
    This function performs the angular spectrum propagation

    Args:
        Uin: Input field for angular spectrum propagation (x1,y1)
        wvl: Wavelength [m]
        d1: Sampling interval at source plane [m]
        d2: Sampling interval at observation plane [m]
        Dz: Propagation distance [m]
    Returns:
        Uout: Complex field at the observation plane (x2,y2)
        x2: 2D meshgrid of x coordinates at observation
        y2: 2D meshgrid of y coordinates at observation
    """
    N = np.shape(Uin)[0]
    k = (2*np.pi)/wvl
    m = d2/d1
    df1 = 1/(N*d1)
    x1,y1 = square_meshgrid(N,d1)
    fX,fY = square_meshgrid(N,df1)
    x2,y2 = square_meshgrid(N,d2)
    Q3 = np.exp(1j*(k/2)*((m-1)/(m*Dz))*(x2**2+y2**2))
    Q2 = np.exp(-1j*np.pi**2*((2*Dz)/(k*m))*(fX**2+fY**2))
    Q1 = np.exp(1j*(k/2)*((1-m)/(Dz))*(x1**2+y1**2))
    Uout = Q3*ift2(Q2*ft2(Q1*Uin, d1), df1)
    return Uout,x2,y2