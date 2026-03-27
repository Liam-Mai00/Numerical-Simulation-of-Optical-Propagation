"""
This module contains useful transformations
"""

import numpy as np

def cart2pol(x, y):
    """
    Transforms from Cartesian to Polar geometry

    Args:
        x: x coordinates
        y: y coordinates
    Returns:
        rho: radius
        phi: angle
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)

def pol2cart(rho, phi):
    """
    Transforms from Polar to Cartesian geometry

    Args:
        rho: radius
        phi: angle
    Returns:
        x: x coordinates
        y: y coordinates
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)