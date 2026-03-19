"""
Discrete Fourier Transform functions from Chapter 2
When numerically calculating the DFT (through numpy or MATLAB libraries), factoring in the sampling interval is omitted.
In this case, it is factored in. For more information, check the dft_scaling.ipynb file under Notes.
"""

from numpy.fft import fft, fft2, fftshift, ifft, ifft2

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
    x = fftshift(ifft(fftshift(X))) * len(X) * delta_f
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
    x = fftshift(ifft2(fftshift(X))) * (N*delta_f)**2