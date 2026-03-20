"""
This is a python demo of a 1D Fourier Transform.
"""
import numpy as np
import matplotlib.pyplot as plt

DISPLAY_FIGS = True

# Signal parameters
ampl = 1                    # Signal Amplitude
freq = 2                    # Signal frequency [Hz]
samp_freq = 100             # Sampling frequency [Hz]
samp_period = 1/samp_freq   # Sampling preiod [s]
phase = 0
duration = 3                # Total signal duration [s]

# Create cosine wave
t = np.arange(0,duration,samp_period)
signal = ampl*np.cos(2*np.pi*freq*t + phase)
f_signal = np.fft.fft(signal)
f = np.fft.fftfreq(len(f_signal),samp_period)

# Plot fft
if DISPLAY_FIGS:
    fig1,ax1 = plt.subplots()
    ax1.plot(t,signal,marker='x')
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel(r"$f(t)$",math_fontfamily='cm')
    ax1.set_title(f"Discrete time signal at {freq}Hz")

    fig2,ax2 = plt.subplots()
    ax2.plot(f,np.abs(f_signal))
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel(r"$\scrF\{f(t)\}$",math_fontfamily='cm')
    ax2.set_title(r"Frequency spectrum of $f(t)$",math_fontfamily='cm')

    plt.show()