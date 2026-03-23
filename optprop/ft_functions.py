"""
Discrete Fourier Transform functions from Chapter 2
When numerically calculating the DFT (through numpy or MATLAB libraries), factoring in the sampling interval is omitted.
In this case, it is factored in. For more information, check the dft_scaling.ipynb file under Notes.
"""
import numpy as np
from numpy.fft import fft, fft2, fftshift, ifft, ifft2, ifftshift

def ft(x,delta):
    """
    Performs a 1D DFT on array x and scales with the sampling interval

    Args:
        x: Input array in spacial domain
    Returns:
        X: Scaled DFT output of x
    """
    X = fftshift(fft(fftshift(x))) * delta
    return X

def ft2(x,delta):
    """
    Performs a 2D DFT on array x and scales with the sampling interval

    Args:
        x: 2D Input array in spacial domain
    Returns:
        X: Scaled 2D DFT output of x
    """
    X = fftshift(fft2(fftshift(x))) * (delta**2)
    return X

def ift(X,delta_f):
    """
    Performs a 1D IDFT on array x and scales with the sampling interval

    Args:
        x: Input array in frequency domain
    Returns:
        X: Scaled IDFT output of x
    """
    N = len(X)
    x = ifftshift(ifft(ifftshift(X))) * (N*delta_f) # Uses the ifftshift to cancel effects of fftshift of ft when N is odd
    return x

def ift2(X,delta_f):
    """
    Performs a 2D IDFT on array x and scales with the sampling interval

    Args:
        x: 2D Input array in frequency domain
    Returns:
        X: Scaled 2D IDFT output of x
    """
    N = len(X)
    x = ifftshift(ifft2(ifftshift(X))) * (N*delta_f)**2 # Uses the ifftshift to cancel effects of fftshift of ft2 when N is odd
    return x

def myconv(A,B,delta):
    """
    Performs a 1D convolution of A and B using FT technique

    Args:
        A: Input array A
        B: Input array B
    Returns:
        C: Convolution of A and B using FT method
    """
    N = len(A)
    C = ift(ft(A,delta) * ft(B,delta), 1/(N*delta))
    return C

def myconv2(A,B,delta):
    """
    Performs a 2D convolution of A and B using FT technique

    Args:
        A: Input array A
        B: Input array B
    Returns:
        C: 2D Convolution of A and B using FT method
    """
    N = np.shape(A)[0]
    var1 = ft2(A,delta)
    var2 = ft2(B,delta)
    delta_f = 1/(N*delta)
    C = ift2(ft2(var1 * var2, delta),delta_f)
    return C