###########################
#    Iva Laginja, 2016    #
###########################
"""
Making gridded light curves and ring system images from the data in the folder 'data', created by 'data_creation.py'.

Set variable 'data' to the folder you want to use the light curves and ring systems from.
Uses 'draw_rings.py'.
"""

import sys
import getopt
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import exomoons_py.modules.bring as bring
import exomoons_py.main.draw_rings as draw_rings
from exomoons_py.config import CONFIG_INI


if __name__ == '__main__':

    print('##################  STARTING  #################################\n')

    ######################################################################################
    local_path = CONFIG_INI.get('data_paths', 'local_data_path')
    current_exp = CONFIG_INI.get('data_paths', 'curr_data_path')

    data = os.path.join(local_path, current_exp)  # which data folder should be used for the plots
    ######################################################################################
    print('Working in folder: ', data, '\n')

    # Keeping track of code runtime
    start_time = time.time()

    days_side = 100  # how many days left and right from eclipse midpoint I want to go; # 100 for b=17%; 30 for scaled, 50 for M_Earth, 100 for b=30%
    y_range = 42  # half of height of ring system plots; # 42 for b=17%; 20 for scaled, 35 for M_Earth at 20% R_Hill, 60 for b=30%
    # fig_title = "both m and a from 3 *inner* moons / b = 17.2% Hill sphere radius, m = galilean *5, a = supergalilean *40, v = 13.3 km/s, Hill sphere half filled"
    fig_title = "setting_up"

    #################
    ### Define grid
    #################
    fig_grid = plt.figure()  # create the figure

    plt.suptitle(fig_title)  # include title of the whole thing

    gs1 = GridSpec(3, 3)  # first row of grids (i = 60)
    gs1.update(top=0.98, bottom=0.67, hspace=0.01, wspace=0.07)
    ax_rings_1 = plt.subplot(gs1[0:2, 0])
    ax_light_1 = plt.subplot(gs1[2, 0])
    ax_rings_2 = plt.subplot(gs1[0:2, 1])
    ax_light_2 = plt.subplot(gs1[2, 1])
    ax_rings_3 = plt.subplot(gs1[0:2, 2])
    ax_light_3 = plt.subplot(gs1[2, 2])

    gs2 = GridSpec(3, 3)  # second row of grids (i = 45)
    gs2.update(top=0.65, bottom=0.34, hspace=0.01, wspace=0.07)
    ax_rings_4 = plt.subplot(gs2[0:2, 0])
    ax_light_4 = plt.subplot(gs2[2, 0])
    ax_rings_5 = plt.subplot(gs2[0:2, 1])
    ax_light_5 = plt.subplot(gs2[2, 1])
    ax_rings_6 = plt.subplot(gs2[0:2, 2])
    ax_light_6 = plt.subplot(gs2[2, 2])

    gs3 = GridSpec(3, 3)  # third row of grids (i = 30)
    gs3.update(top=0.32, bottom=0.02, hspace=0.01, wspace=0.07)
    ax_rings_7 = plt.subplot(gs3[0:2, 0])
    ax_light_7 = plt.subplot(gs3[2, 0])
    ax_rings_8 = plt.subplot(gs3[0:2, 1])
    ax_light_8 = plt.subplot(gs3[2, 1])
    ax_rings_9 = plt.subplot(gs3[0:2, 2])
    ax_light_9 = plt.subplot(gs3[2, 2])

    # creating lists for help in axis manipulation (labels, ticks, axis limits)
    ax_light_list = [ax_light_1, ax_light_2, ax_light_3, ax_light_4, ax_light_5, ax_light_6, ax_light_7, ax_light_8,
                     ax_light_9]
    ax_rings_list = [ax_rings_1, ax_rings_2, ax_rings_3, ax_rings_4, ax_rings_5, ax_rings_6, ax_rings_7, ax_rings_8,
                     ax_rings_9]
    no_x_list = [ax_rings_1, ax_rings_2, ax_rings_3, ax_rings_4, ax_rings_5, ax_rings_6, ax_rings_7, ax_rings_8, ax_rings_9,
                 ax_light_1, ax_light_2, ax_light_3, ax_light_4, ax_light_5, ax_light_6]
    no_y_list = [ax_rings_1, ax_rings_2, ax_rings_3, ax_rings_4, ax_rings_5, ax_rings_6, ax_rings_7, ax_rings_8, ax_rings_9,
                 ax_light_2, ax_light_3, ax_light_5, ax_light_6, ax_light_8, ax_light_9]

    # name_axes(plt.gcf())
    bring.make_x_ticklabels_invisible(plt.gcf(), no_x_list)
    bring.make_y_ticklabels_invisible(plt.gcf(), no_y_list)
    bring.make_x_labels(plt.gcf(), (ax_light_7, ax_light_8, ax_light_9))
    bring.make_y_labels(plt.gcf(), (ax_light_1, ax_light_4, ax_light_7))
    bring.x_tick_dist(plt.gcf(), ax_light_list, 20.)  # the number in here defines the interval for the x ticks
    bring.set_axes_limits(plt.gcf(), -days_side, days_side, 0.91, 1.02, ax_light_list)

    print('---------------------------')
    print('MAKING LIGHT CURVES:\n')
    ########################
    ### the light curves
    ########################

    for i, axis in enumerate(ax_light_list):
        print("Making light curve #", i + 1)
        print('')
        # Importing the light curve data
        time, flux, flux_errxx = bring.lightcurve_import(os.path.join(data, 'model' + str(i) + '.dat'), 'time', 'flux', 'flux_rms')
        time_sc, flux_sc, flux_errxx_sc = bring.lightcurve_import(os.path.join(data, 'model' + str(i) + '_scat' + '.dat'), 'sc_deg',
                                                                  'sc_light', 'pseudo_error')

        # interpolation for 5min-sampling
        new_time, new_flux = bring.interpol_5min(time, flux)
        new_time_sc, new_flux_sc = bring.interpol_5min(time_sc, flux_sc)

        # introduce gaussian scatter
        flux_noise = bring.gaussian_scatter(new_flux)
        flux_noise_sc = bring.gaussian_scatter(new_flux_sc)

        # bin to 1h
        bin_time_1h, bin_means_1h = bring.binned_1h(new_time, flux_noise)
        bin_time_1h_sc, bin_means_1h_sc = bring.binned_1h(new_time_sc, flux_noise_sc)

        # Introduce stellar pulsations to light curve (Koen at al. 2003a)
        '''
        freq1 = 47.055 # [c/d] - cycles per day
        freq2 = 38.081 # [c/d]
        freq3 = 52.724 # [c/d]
    
        amp1 = 1.63 # [mmag] (B filter)
        amp2 = 1.50 # [mmag] (B filter)
        amp3 = 1.07 # [mmag] (B filter)
    
        amp_array = np.array([amp1*0.001,amp2*0.001,amp3*0.001]) # need to make it compatible with the transmission scale on the ligth curves (thus the factor of 0.001)
        freq_array = np.array([freq1,freq2,freq3])
        
        bin_means_1h,sine_tot = bring.stellar_puls(bin_time_1h,bin_means_1h,amp_array,freq_array)
        '''

        # Introduce forward scattering
        '''
        timexx, scat_flux, scat_flux_errxx = bring.lightcurve_import(os.path.join(data, 'model' + str(i) + '_scat.dat'), 'sc_deg', 'sc_light', 'pseudo_error')
        
        diff = scat_flux.size - bin_means_1h.size       # Original light curves don't have same dimensions,
        side = (diff-1)/2                               # so here I make sure I cut them to the same size
        scat_flux = scat_flux[side+1:-side]-0.3         # while still keeping them centered
    
        bin_means_1h = bin_means_1h + scat_flux
        print bin_means_1h[5000:5100]
        '''

        # Plot the light curve
        bring.plot_1h_sampling(axis, bin_time_1h, bin_means_1h, new_time, new_flux,
                               days_side)  # this fills the light curve grids

        # add scattered curve
        ax_add = axis.twinx()
        ax_add.plot(new_time, new_flux + new_flux_sc)
        ax_add.set_ylim(60., 100.)
        ax_add.set_xlim(-days_side, days_side)

    print('---------------------------')
    print('MAKING RING PLOTS:\n')
    ########################
    ### the ring systems
    ########################

    for n, plot in enumerate(ax_rings_list):
        print("Making ring system #", n + 1)
        print('')
        in_fits = os.path.join(data, 'model' + str(n) + '.fits')
        draw_rings.plotting_rings(plot, in_fits, days_side, -y_range, y_range, data)  # this fills in the ring system grids

    # label the tips and tilts
    ax_light_1.text(-(days_side - 2), 0.96, '$\phi=0^{\circ}$\n$i=60^{\circ}$', fontsize=20)
    ax_light_2.text(-(days_side - 2), 0.96, '$\phi=30^{\circ}$\n$i=60^{\circ}$', fontsize=20)
    ax_light_3.text(-(days_side - 2), 0.96, '$\phi=60^{\circ}$\n$i=60^{\circ}$', fontsize=20)
    ax_light_4.text(-(days_side - 2), 0.96, '$\phi=0^{\circ}$\n$i=45^{\circ}$', fontsize=20)
    ax_light_5.text(-(days_side - 2), 0.96, '$\phi=30^{\circ}$\n$i=45^{\circ}$', fontsize=20)
    ax_light_6.text(-(days_side - 2), 0.96, '$\phi=60^{\circ}$\n$i=45^{\circ}$', fontsize=20)
    ax_light_7.text(-(days_side - 2), 0.96, '$\phi=0^{\circ}$\n$i=30^{\circ}$', fontsize=20)
    ax_light_8.text(-(days_side - 2), 0.96, '$\phi=30^{\circ}$\n$i=30^{\circ}$', fontsize=20)
    ax_light_9.text(-(days_side - 2), 0.96, '$\phi=60^{\circ}$\n$i=30^{\circ}$', fontsize=20)

    print('-------------------------------')
    print('PLOT SAVING TO FOLDER: ', data)
    print('-------------------------------')

    plt.show()
    fig_grid.savefig(os.path.join(data, current_exp + ".pdf"), bbox_inches='tight')

    #### Ending program
    end_time = time.time()
    print('Runtime:', str((end_time - start_time) / 60) + " min")
