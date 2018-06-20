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
from astropy.table import Table,Column
from astropy.io import ascii

import pyfits

import exorings
import j1407
import bring
import draw_rings
import bring_disk_sim_data


from scipy.optimize import fmin
#from scipy.ndimage import convolve
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from matplotlib.patches import PathPatch
from scipy import stats

from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.gridspec import GridSpec

# constants (taken from disk_sim.py)
G = 6.6738480e-11 # m3 kg-1 s-2
yr = 365.242189669 * 86400  # sec
msol = 1.98855e30 # kg
rsol = 6.5500e8 # m
mjup = 1.8986e27 # kg
rjup = 6.9911e7 # m
mearth = 5.97219e24 # kg
rearth = 6.371e3 # kg (Wikipedia)
mmoon = 7.3476e22 # kg
au = 1.49597870700e11 # m
pc = 3.0856e16 # m

# -------------------------------------------- #
# from data_creation.py

i_deg = np.arange(1.,90.,5.)
phi_deg = np.arange(0.,90.,5.)

#################
#i_deg = 45.         # physically: i = 90-i_deg
#phi_deg = 40.
#################

targetdir = "test_single/many/"   # don't forget the / at the end!!!
data = 'test_single/many'
impact = 0.1722       # impact parameter in terms of parts of the Hill radius => b = impact * R_Hill

# If targetdir doesn't exist yet, create it
if os.path.isdir(targetdir) == False:
    os.mkdir(targetdir)

for i,n in enumerate(i_deg):
    for k,m in enumerate(phi_deg):

        # Vorlage: bring_disk_data(impact,i_in,phi_in,output,targetdir,modelnumber,paramfile=False)

        # create data
        print ''
        print "MAKING MODEL 0 NOW:"
        bring_disk_sim_data.bring_disk_data(impact,n,m,"model_i"+str(n)+"_phi"+str(m),targetdir,str(0),paramfile=False)


        print ''
        print '---------------------------'
        print 'impact = ',impact
        print 'DATA SAVED TO ',targetdir
        print '---------------------------'

        # ----------------------------------------------- #
        # from grid_lightcurve_plots.py

        days_side = 100   # how many days left and right from eclipse midpoint I want to go; # 100 for b=17%; 30 for scaled, 50 for M_Earth, 100 for b=30%
        y_range = 42     # half of height of ring system plots; # 42 for b=17%; 20 for scaled, 35 for M_Earth at 20% R_Hill, 60 for b=30%
        fig_title = "both m and a from 3 *inner* moons / b = 17.2% Hill sphere radius, m = galilean *5, a = supergalilean *40, v = 13.3 km/s, Hill sphere half filled"

        fig_grid = plt.figure()   # create the figure
        #fig_grid.tight_layout()
        #plt.suptitle(fig_title)   # include title of the whole thing

        gs1 = GridSpec(3, 1)   # first row of grids (i = 60)
        ax_rings_1 = plt.subplot(gs1[0:2,0])
        ax_light_1 = plt.subplot(gs1[2,0])

        #####################
        ### Make light curve
        axis = ax_light_1

        # Importing the light curve data
        time, flux, flux_errxx = bring.lightcurve_import(data + "/model_i"+str(n)+"_phi"+str(m) + '.dat','time','flux','flux_rms')

        # interpolation for 5min-sampling
        new_time, new_flux = bring.interpol_5min(time, flux)

        # introduce gaussian scatter
        flux_noise = bring.gaussian_scatter(new_flux)

        # bin to 1h
        bin_time_1h, bin_means_1h = bring.binned_1h(new_time,flux_noise)

        # plot the light curve
        bring.plot_1h_sampling(axis,bin_time_1h,bin_means_1h,new_time,new_flux,days_side)   # this fills the light curve grids

        ######################
        ### Make ring system
        plot = ax_rings_1

        in_fits = data + "/model_i"+str(n)+"_phi"+str(m) + '.fits'
        draw_rings.plotting_rings(plot,in_fits,days_side,-y_range,y_range,data)   # this fills in the ring system grids

        print '-------------------------------'
        print 'PLOT SAVING TO FOLDER: ', data
        print '-------------------------------'

        fig_grid.savefig(data + "/" + "model_i"+str(n)+"_phi"+str(m) + ".pdf", bbox_inches='tight')
        #plt.show(fig_grid)






