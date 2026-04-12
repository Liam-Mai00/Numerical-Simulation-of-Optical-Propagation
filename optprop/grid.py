import numpy as np

def square_meshgrid(N,delta):
    """
    This function addresses the square meshgrid repitition problem
    
    Args:
        N: Number of samples
        delta: Spacing in spacial or frequency domain
    Returns:
        X: 2D mesh of x coordinates
        Y: 2D mesh of y coordinates
    """
    x = np.arange(-N/2,N/2,1) * delta
    X,Y = np.meshgrid(x,x)
    return X,Y