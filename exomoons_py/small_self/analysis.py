"""
Data analysis for the light curves. Meaning, I can extract single features and run the test if it is a
significant detection.
"""

import sys, getopt, os

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column

import exomoons_py.modules.bring as bring

print('##################  STARTING  #################################\n')

######################################################################################
data = "b17perc_m-gal_a-special_v13300_halfHill_3_gapsplit"  # which data folder should be used for the plots
modelnr = 4  # which model number to use from that folder (range from 0 to 8; careful with numbering of models (data files in folder, 0-8) and numbering of light curves and ring systems in plots (1-9))
######################################################################################
print('Working in folder: ', data, '\n')

# the following parameters should be copied from "grid_lightcurve_plots.py"
# from here --
days_side = 30  # how many days left and right from ecliipse midpoint I want to go # 30 for scaled, 5 for M_Earth
y_range = 20  # half of height of ring system plots # 20 for scaled, 35 for M_Earth at 20% R_Hill
fig_title = "b = 3% Hill sphere radius, m = M_Earth, a = 20% Hill sphere radius"
# to here --

sigma = 5  # confidence interval to check detection
controlpt = 8  # how many light curve points to check back with the last 24h

### analyze the light curves

# for i,axis in enumerate(ax_light_list):   #temp: repalce modelnr with i when going for loop
print("Analyzing light curve #", modelnr + 1, ' (model', modelnr, ")")
print('')
# Importing the light curve data
time, flux, flux_errxx = bring.lightcurve_import(data + '/model' + str(modelnr) + '.dat')
# interpolation for 5min-sampling
new_time, new_flux = bring.interpol_5min(time, flux)

# introduce gaussian scatter
flux_noise = bring.gaussian_scatter(new_flux)

# bin to 1h
bin_time_1h, bin_means_1h = bring.binned_1h(new_time, flux_noise)

### new stuff

true_arr = np.zeros(len(bin_means_1h) - 30)

for n in range(len(bin_means_1h) - 30):
    mean24 = np.mean(bin_means_1h[n:(n + 24)])
    rms24 = np.std(bin_means_1h[n:(n + 24)])

    mean5 = np.mean(bin_means_1h[(n + 24):(n + 24 + controlpt)])
    totval = mean24 + sigma * rms24
    if (mean5 > totval):
        test = True
    else:
        test = False

    true_arr[n] = test

woisndas = np.where(true_arr == True)  # find where it's True
overlap = (true_arr * 1) + 0.94  # convert to 1 and 0 instead of True and False and make the 0s to 0.94
overlap[woisndas] = 1.  # ones should stay ones

plt.scatter(bin_time_1h, bin_means_1h, s=5)  # the full light curve
# plt.plot(bin_time_1h[:-30],overlap,color='r',linewidth=2.0)   # detection indication overlaid
plt.tick_params(width=2, length=10)
plt.xlim(50., 60.)
plt.ylim(0.92, 1.02)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlabel('Days', fontsize=30)
plt.ylabel('Transmission', fontsize=30)
plt.title("feature_4")
fig1 = plt.gcf()  # Save the figure before you show() by calling plt.gcf() for "get current figure", then you can call savefig() on this Figure object at any time. Because after plt.show() is called, a new figure is created.
plt.show()
# fig1.savefig(data + "/detection_" + str(modelnr) + '_sig' + str(sigma) + '-pt' + str(controlpt) + ".pdf", format='pdf', bbox_inches='tight')
fig1.savefig(data + "/feature4_model4.pdf", format='pdf', bbox_inches='tight')
