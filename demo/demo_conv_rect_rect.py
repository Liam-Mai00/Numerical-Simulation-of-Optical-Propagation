"""
2D convolution between two rectangular functions
"""
import numpy as np
import matplotlib.pyplot as plt
from optprop.ft_functions import myconv
from optprop.aperture_functions import rect,tri

#### Configure Optical Parameters
N = 64         # Number of samples
L = 8          # Spacial extent of grid [m]
delta = L/N     # Spacial sample spacing [m]
F = 1/(N*delta) # Frequency sample spacing [1/m]
w = 2           # Width of rectangle [m]

#### Perform 2D Convolution
x = np.arange(-N/2,N/2,1) * delta
A = rect(x/w)
print(f"A={A}\n")
B = A
print(f"B={B}\n")
C = myconv(A,B,delta)
print(f"C={C}")
C_cont = w * tri(x/w)

#### Display Figs
fig,ax = plt.subplots(1,2)
ax[0].plot(A)
ax[0].set_xlabel("$x[m]$",math_fontfamily='cm')
ax[0].set_ylabel("$y[0]$",math_fontfamily='cm')
ax[1].plot(B)

fig,ax = plt.subplots()
ax.scatter(x,C,marker='x',label="Real C")
ax.scatter(x,np.abs(C),edgecolors='red',marker='s',facecolors='none',label="Magnitude C")
ax.scatter(x,C_cont,edgecolors='green',marker='o',facecolors='none',label="C_cont")
ax.set_ylim(0,2.25)
ax.set_xlim(-4,4)
ax.legend()
plt.show()