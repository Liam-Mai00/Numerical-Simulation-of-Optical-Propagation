import numpy as np
import matplotlib.pyplot as plt
from optprop.grid import square_meshgrid
from optprop.aperture_functions import rect
from optprop.prop import one_step_prop

"""
Demonstrates the one-step propagation of a square aperture 
"""
D1 = 2e-3               # Diameter of source aperture [m]
D2 = 3e-3               # Diameter of observation region of interest [m]
delta1 = D1/50          # Source grid spacing [m]
wvl = 1e-6              # Optical wavelength [m]
k = (2*np.pi)           # Wavenumber
Dz = 0.5                # Propagation distance [m]
Nmin = D1*wvl*Dz/(delta1*(wvl*Dz - D2*delta1))    # Minimum number of sample points required (Equation 7.31)
N = (2**np.ceil(np.log2(Nmin))).astype(int)         # Makes N the next available power of 2
x1,y1 = square_meshgrid(N,delta1)
ap = rect(x1/D1) * rect(y1/D1)                      # Create square distributed field
Uout,x2,y2 = one_step_prop(ap,wvl,delta1,Dz)

fig,ax = plt.subplots(1,2)
ax[0].plot(x2[0,:],np.abs(Uout[N//2,:])**2,marker='x')
ax[0].set_xlabel(r"$x_2[m]$",math_fontfamily='cm')
ax[0].set_ylabel("Irradiance")
ax[0].set_xlim(-1.5e-3,1.5e-3)
ax[0].set_ylim(0,1.2)
ax[0].ticklabel_format(axis='x',style='sci',scilimits=(-3,-3))

ax[1].plot(x2[0,:],np.angle(Uout[N//2,:]),marker='x')
ax[1].set_xlabel(r"$x_2[m]$",math_fontfamily='cm')
ax[1].set_ylabel("Irradiance")
ax[1].set_xlim(-1.5e-3,1.5e-3)
ax[1].set_ylim(-0.5,2.5)
ax[1].ticklabel_format(axis='x',style='sci',scilimits=(-3,-3))
plt.show()