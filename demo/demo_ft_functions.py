import numpy as np
import matplotlib.pyplot as plt
from optical_functions.ft_functions import ft

# 1D Case
ampl = 1                    # Signal Amplitude
freq = 10                   # Signal frequency [Hz]
samp_freq = 100             # Sampling frequency [Hz]
samp_period = 1/samp_freq   # Sampling preiod [s]
phase = 0
duration = 5                # Total signal duration [s]

t = np.arange(0,duration,samp_period)
signal = ampl*np.cos(2*np.pi*freq*t + phase)
f_signal = ft(signal,samp_period)

fig1,ax1 = plt.subplots()
ax1.plot(t,signal)

fig2,ax2 = plt.subplots()
ax2.plot()

plt.show()