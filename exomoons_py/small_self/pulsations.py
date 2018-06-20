"""
Test script for including stellar pulsations in the beta Pictoris exoring light curve simulations. Creates a sine wave when run.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.io import ascii
import matplotlib as mpl
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt


def stellar_puls(time, flux, amplitude, frequency):
    """Generates the sine signal to overlay with the light curve."""
    t = np.linspace(np.min(time), np.max(time),
                    flux.size)  # create array equivalent in size and min and max on x-axis like the input flux
    phi = np.pi / 2.  # define a phase (can be anything, it just shifts the sine wave on the x-axis)
    sine = amplitude * np.sin(2 * np.pi * frequency * t + phi)  # create the sine wave
    flux_with_pulse = flux + sine  # add the sine wave to the input flux

    return flux_with_pulse


def cd_to_uHz(freq_cd):
    """Transforms the input frequency from cycles per day to microHerz."""
    freq_uHz = freq_cd * 11.574074  # [microHz]

    return freq_uHz


def cd_to_periodDay(freq_cd):
    """Transforms the input frequency from cycles per day to a period in units of days."""
    period_day = 1 / freq_cd  # [days]

    return period_day


### beta Pic pulsation frequencies and amplitudes (Koen at al. 2003a)

freq1 = 47.055  # [c/d] - cycles per day
freq2 = 38.081  # [c/d]
freq3 = 52.724  # [c/d]

amp1 = 1.63  # [mmag] (B filter)
amp2 = 1.50  # [mmag] (B filter)
amp3 = 1.07  # [mmag] (B filter)

fs = 100000  # sampling rate
t = np.linspace(65., -65., fs)
phi = 0
x = amp1 * np.sin(2 * np.pi * freq1 * t + phi)
y = amp2 * np.sin(2 * np.pi * freq2 * t + phi)
z = amp3 * np.sin(2 * np.pi * freq3 * t + phi)

u = x + y + z  # sum of all of the different pulsation frequencies

# plt.plot(t,x)
# plt.plot(t,y)
# plt.plot(t,z)
plt.plot(t, u)
plt.xlabel('time')
plt.ylabel('amplitude')
plt.show()
