from numpy.fft import fft, fft2, fftshift, ifft, ifft2

def ft(x,delta):
    X = fftshift(fft(fftshift(x))) * delta
    return X

def ft2(x,delta):
    X = fftshift(fft2(fftshift(x))) * (delta**2)
    return X

def ift(X,delta_f):
    x = fftshift(ifft(fftshift(X))) * delta_f
    return x