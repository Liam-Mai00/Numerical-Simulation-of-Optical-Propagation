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
N = int(2**np.ceil(np.log2(Nmin)))   # Obtain next power of 2 for fft

#### Propagate ####
x1,y1 = square_meshgrid(N,delta1)
ap = rect(x1/D1) * rect(y1/D1)
Uout,x2,y2 = ang_spec_prop(ap,wvl,delta1,delta2,Dz)

#### Plot Figs ####
fig,ax = plt.subplots(1,2)
ax[0].plot(x2[0,:],np.abs(Uout[N//2,:])**2,marker='x')
ax[0].set_xlabel(r"$x_2[m]$",math_fontfamily='cm')
ax[0].set_ylabel("Irradiance")
ax[0].set_xlim(-2e-3,2e-3)
ax[0].set_ylim(0,1.5)
ax[0].ticklabel_format(axis='x',style='sci',scilimits=(-3,-3))

ax[1].plot(x2[0,:],np.angle(Uout[N//2,:]),marker='x')
ax[1].set_xlabel(r"$x_2[m]$",math_fontfamily='cm')
ax[1].set_ylabel("Phase[rad]")
ax[1].set_xlim(-2e-3,2e-3)
ax[1].ticklabel_format(axis='x',style='sci',scilimits=(-3,-3))

plt.show()