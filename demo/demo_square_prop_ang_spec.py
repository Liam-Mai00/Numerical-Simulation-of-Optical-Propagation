"""
Evaluating Fresnel diffraction integral using angular-spectrum method
Listing 7.2
"""
import numpy as np
import matplotlib.pyplot as plt
from optprop.grid import square_meshgrid
from optprop.aperture_functions import rect
from optprop.prop import ang_spec_prop

#### Define Propagation Parameteres ####
D1 = 2e-3   # Diameter of the source aperture [m]
D2 = 4e-3   # Diameter of the observation aperture [m]
wvl = 1e-6  # Optical wavelength [m]
k = (2*np.pi)/wvl   # Wavenumber
Dz = 0.1            # Propagation distance [m]
delta1 = 9.4848e-6  # Source grid spacing [m]
delta2 = 28.1212e-6 # Observation grid spacing [m]
Nmin = D1/(2*delta1) + D2/(2*delta2) + \
    (wvl*Dz)/(2*delta1*delta2)  # Minimum number of samples
N = 2**np.ceil(np.log2(Nmin))   # Obtain next power of 2 for fft

#### Propagate ####
x1,y1 = square_meshgrid(N,delta1)
ap = rect(x1/D1) * rect(y1/D2)
Uout,x2,y2 = ang_spec_prop(ap,wvl,delta1,delta2,Dz)

#### Plot Figs ####
fig,ax = plt.subplots()
ax.plot(x2[0,:],np.abs(Uout[int(N)//2,:])**2)

plt.show()
a = 0