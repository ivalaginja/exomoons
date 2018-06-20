###########################
#    Iva Laginja, 2016    #
# last change, 03-04-2017 #
###########################

"""
This code combines previous codes to produce one single rings system and light curve plot.
i copied the important parts from 'bring_disk_sim_data.py' and 'grid_lightcurve_plots.py' into here.
"""

import sys, getopt, os

sys.path.append('/Users/kenworthy/Dropbox/python_workbooks/lib')
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from matplotlib.gridspec import GridSpec

import exomoons_py.modules.bring as bring
import exomoons_py.main.draw_rings as draw_rings
import exomoons_py.main.bring_disk_sim_data as disksim

# constants (taken from disk_sim.py)
G = 6.6738480e-11  # m3 kg-1 s-2
yr = 365.242189669 * 86400  # sec
msol = 1.98855e30  # kg
rsol = 6.5500e8  # m
mjup = 1.8986e27  # kg
rjup = 6.9911e7  # m
mearth = 5.97219e24  # kg
rearth = 6.371e3  # kg (Wikipedia)
mmoon = 7.3476e22  # kg
au = 1.49597870700e11  # m
pc = 3.0856e16  # m

# -------------------------------------------- #
# from data_creation.py

#################
i_deg = 72.  # physically: i = 90-i_deg
phi_deg = 32.
#################

targetdir = "thin_Lecav/"  # don't forget the / at the end!!!
data = 'thin_Lecav'
impact = 0.1722  # impact parameter in terms of parts of the Hill radius => b = impact * R_Hill

# If targetdir doesn't exist yet, create it
if os.path.isdir(targetdir) == False:
    os.mkdir(targetdir)

# Vorlage: bring_disk_data(impact,i_in,phi_in,output,targetdir,modelnumber,paramfile=False)

# create data
print('')
print("MAKING MODEL 0 NOW:")
norm_val = disksim.bring_disk_data(impact, i_deg, phi_deg, "model0", targetdir, str(0), paramfile=True)

print('\n---------------------------')
print('impact = ', impact)
print('DATA SAVED TO ', targetdir)
print('---------------------------')

# ----------------------------------------------- #
# from grid_lightcurve_plots.py

days_side = 100  # how many days left and right from eclipse midpoint I want to go; # 100 for b=17%; 30 for scaled, 50 for M_Earth, 100 for b=30%
y_range = 42  # half of height of ring system plots; # 42 for b=17%; 20 for scaled, 35 for M_Earth at 20% R_Hill, 60 for b=30%
fig_title = "both m and a from 3 *inner* moons / b = 17.2% Hill sphere radius, m = galilean *5, a = supergalilean *40, v = 13.3 km/s, Hill sphere half filled"

fig_grid = plt.figure()  # create the figure
# fig_grid.tight_layout()
# plt.suptitle(fig_title)   # include title of the whole thing

gs1 = GridSpec(3, 1)  # first row of grids (i = 60)
ax_rings_1 = plt.subplot(gs1[0:2, 0])
ax_light_1 = plt.subplot(gs1[2, 0])

#####################
### Make light curve
axis = ax_light_1

# Importing the light curve data
time, flux, flux_errxx = bring.lightcurve_import(data + '/model' + str(0) + '.dat', 'time', 'flux', 'flux_rms')
time_sc, flux_sc, flux_errxx_sc = bring.lightcurve_import(data + '/model' + str(0) + '_scat' + '.dat', 'sc_deg',
                                                          'sc_light', 'pseudo_error')

# testing
# plt.plot(time,flux)
# plt.plot(time_sc,flux_sc)
# plt.show()

# interpolation for 5min-sampling
new_time, new_flux = bring.interpol_5min(time, flux)
new_time_sc, new_flux_sc = bring.interpol_5min(time_sc, flux_sc)

# introduce gaussian scatter
flux_noise = bring.gaussian_scatter(new_flux)
flux_noise_sc = bring.gaussian_scatter(new_flux_sc)

# bin to 1h
bin_time_1h, bin_means_1h = bring.binned_1h(new_time, flux_noise)
bin_time_1h_sc, bin_means_1h_sc = bring.binned_1h(new_time_sc, flux_noise_sc)

# add scattered and occulted light curves
sum_light = bin_means_1h + bin_means_1h_sc

# plot the light curve
# bring.plot_1h_sampling(axis,bin_time_1h,bin_means_1h,new_time,new_flux,days_side)   # this fills the light curve grids

# add scattered light on extra axis but same plot (yes that works)
# ax_add = axis.twinx()
# ax_add.plot(new_time,new_flux_sc,color='blue')
# ax_add.set_ylim(0.,2.)
# ax_add.set_xlim(-days_side,days_side)
# ax_add.set_xlim(-days_side+150,days_side-45)

# add added light on extra axis but same plot (yes that works)
# ax_add = axis.twinx()
# ax_add.plot(new_time,new_flux+new_flux_sc,color='green')
# ax_add.set_ylim(0.,2.)
# ax_add.set_xlim(-days_side,days_side)
# ax_add.set_xlim(-days_side+150,days_side-45)

flux_tot = bin_means_1h + bin_means_1h_sc
curve_tot = new_flux + new_flux_sc

# plot the added curve
bring.plot_1h_sampling(axis, bin_time_1h_sc, flux_tot, new_time, curve_tot,
                       days_side)  # this fills the light curve grids

######################
### Make ring system
plot = ax_rings_1

in_fits = data + '/model' + str(0) + '.fits'
draw_rings.plotting_rings(plot, in_fits, days_side, -y_range, y_range, data)  # this fills in the ring system grids

print('-------------------------------')
print('PLOT SAVING TO FOLDER: ', data)
print('-------------------------------')

fig_grid.savefig(data + "/" + data + ".pdf", bbox_inches='tight')
plt.show(fig_grid)
