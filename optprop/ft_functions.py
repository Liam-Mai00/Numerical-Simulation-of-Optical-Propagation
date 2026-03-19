"""
Discrete Fourier Transform functions from Chapter 2
When numerically calculating the DFT (through numpy or MATLAB libraries), factoring in the sampling interval is omitted.
In this case, it is factored in. For more information, check the dft_scaling.ipynb file under Notes.
"""

from numpy.fft import fft, fft2, fftshift, ifft, ifft2

def ft(x,delta):
    X = fftshift(fft(fftshift(x))) * (delta)
    return X

def ft2(x,delta):
    X = fftshift(fft2(fftshift(x))) * (delta**2)
    return X

def ift(X,delta_f):
    N = len(X)
    x = fftshift(ifft(fftshift(X))) * (N*delta_f)
    return x

def ift2(X,delta_f):
    N = len(X)
    x = fftshift(ifft2(fftshift(X))) * (N*delta_f)**2