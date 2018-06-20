"""
This is the ring plotting program for bring_disk_sim.py
Changes from previous versions: #ivanew

--------------------------------------------------------
USED BY grid_lightcurve_plots.py
--------------------------------------------------------
"""

# some specific changes should be reversed at time again - #ChMade

import sys, getopt
import numpy as np
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
from astropy.table import Table, Column
from astropy.io import ascii

import exomoons_py.modules.exorings as exorings
import exomoons_py.modules.j1407 as j1407


def plotting_rings(rplot, in_fits, days_side, y_min_ring, y_max_ring, targetdir):
    # --------------------------------------------------------------------
    # read in eclipse time line from the file created by "bring_disk_sim.py"
    time_line = Table.read(targetdir + '/time_params_bring.dat', format='ascii')

    eclipse_start = time_line['start'][0]  # [5 minute intervals]
    eclipse_end = time_line['end'][0]  # [5 minute intervals from time array]
    eclipse_center = 0.  # time_line['center'][0]     # As I am not in Julian Days or anything, I will keep the time array centered on zero and let everything be described as to how far away it is from the time of central eclipse, which is 0.
    eclipse_duration = time_line['duration'][0]  # [days]

    ####################################################setval
    trange = 400.  # [days] twice of how much left and right from the center point I want to go in plot (but with the data, not with plot size)
    ####################################################
    # --------------------------------------------------------------------

    # print ''
    # print "Hill eclipse start in imported data:",eclipse_start
    # print "Hill eclipse end in imported data:",eclipse_end
    # print "eclipse center in imported data:",eclipse_center
    print("Hill eclipse duration:", eclipse_duration)
    print("The moons' orbits cover only the very inner parts of the Hill sphere, this part will be shown by the ring system.")
    print("")

    # no scientific notation for numbers on plots
    mpl.rc('axes.formatter', limits=(-7, 7))

    # range of days that ring statistics should be considered
    rings_tmin = eclipse_center - (trange / 2.)
    rings_tmax = eclipse_center + (trange / 2.)

    # read in photometry
    (time, flux, flux_err) = j1407.j1407_photom_binned(targetdir + '/lightcurve_bring.dat', rings_tmin,
                                                       rings_tmax)  # the last two values are the limits in time units as used in the data file between which I want the photometry to be displayed. They are the begining and end time of the eclipse plus a bit of contingency, which here is the part of the data for the rest of the Hill sphere

    print('number of photometric points: %d' % time.size)

    # get gradients
    (grad_time, grad_mag, grad_mag_norm) = j1407.j1407_gradients('../../input_data/j1407_gradients.txt')

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hr:o:s:", ["rfile=", "ofile=", "f="])
    except getopt.GetoptError:
        print('%s -r <inputfile> -s <velocity in metres per second> -o <outputfile>' % sys.argv[0])
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print(help)
            sys.exit()
        elif opt in ("-r", "--rfile"):
            fitsin = arg
        elif opt in ("-o", "--ofile"):
            plotout = arg
        elif opt in ("-s", "--vstar"):
            v = np.array(float(arg))
    print("")
    # print 'Reading in ring and disk parameters from %s' % fitsin

    fitsin = in_fits
    fitsout = "testring.pdf"

    (res, taun_rings, rad_rings, dstar) = exorings.read_ring_fits(fitsin)
    # taun_rings = np.resize([0.,100.],len(rad_rings)) # want to be the rings completely non-transmissive #ivanew
    ###################
    ##################

    # printing  the ring radii and their tau values
    exorings.print_ring_tau(rad_rings, exorings.y_to_tau(taun_rings))

    print('\n RADII IN ARRAY FORM')
    print(rad_rings)

    print('\n Gap sizes in days:')
    for i in range(0, len(rad_rings) - 3, 2):
        print('free gap:', rad_rings[i + 1] - rad_rings[i], ' days')

    print('\n Gap sizes in hours:')
    for i in range(0, len(rad_rings) - 3, 2):
        print('free gap:', (rad_rings[i + 1] - rad_rings[i]) * 24., ' hours')

    # set up stellar disk
    kern = exorings.make_star_limbd(21, 0.8)

    # make the radius and projected gradient for the measured gradient points
    (ring_disk_fit, grad_disk_fit) = exorings.ring_grad_line(grad_time, res[0], res[1], res[2], res[3])

    # produce fine grained gradient and ring values  
    samp_t = np.arange(-100, 100, 0.001) + eclipse_center
    (samp_r, samp_g) = exorings.ring_grad_line(samp_t, res[0], res[1], res[2], res[3])
    hjd_minr = 0.  # samp_t[np.argmin(samp_g)] #pyadjust - it works with 0, so I'll keep it like that for now

    # times when there is no photometry
    (rstart, rend) = exorings.ring_mask_no_photometry(kern, dstar, time, res[0], res[1], res[2], res[3])
    rmask = (rstart < 40)
    rstart = rstart[rmask]
    rend = rend[rmask]

    # print disk parameters
    exorings.print_disk_parameters(res, hjd_minr, samp_r)

    # start making the figure
    # fig_ringsv = plt.figure('ringsv', figsize=(11, 7))   #ivanew

    # window for rings
    # rplot = plt.axes([0.1, 0.3, 0.85, 0.65])   #ivanew
    # rplot.axis('scaled')                       #ivanew

    # window for photometry curve
    # p2v = plt.axes([0.1, 0.11, 0.85, 0.20], sharex=rplot)

    # window for residuals
    # p3v = plt.axes([0.1, 0.06, 0.85, 0.05], sharex=rplot)

    # draw the rings
    exorings.draw_rings_vector(rad_rings, exorings.y_to_tau(taun_rings), res[1], res[2], res[3], rplot)

    # draw the no photometry rings
    exorings.draw_badrings(rstart, rend, res[1], res[2], res[3], rplot)

    # draw the path of the star behind the rings
    star_line = patches.Rectangle((hjd_minr - (trange / 2.), res[0] - dstar / 2.), trange, dstar, color='g',
                                      zorder=-15)  # width of the path should match size of the plot (found further below)
    rplot.add_patch(star_line)

    strip, convo, g = exorings.ellipse_strip(rad_rings, exorings.y_to_tau(taun_rings),
                                             res[0], res[1], res[2], res[3], kern, dstar)
    # strip and dummy never used again in this program - just ignore them

    # error bars on the photometry
    eb = dict(fmt='.', color='white', ecolor='red', capsize=0.0,
              marker='o', mfc='red', mec='red', ms=3, mew=0.001,
              elinewidth=0.5)
    # p2v.errorbar(time, flux, flux_err, zorder=10, **eb)

    fit_time = g[0]
    fit_flux = g[1]

    # p2v.plot(fit_time, fit_flux, linewidth=1, color='green')

    fluxcurve = interp1d(fit_time, fit_flux, kind='linear')
    flux_fit = fluxcurve(time)

    # p3v.errorbar(time, flux-flux_fit, flux_err, zorder=10, **eb)

    # adjust the ticks on the photometry plot
    '''   #ivanew
    for ax in fig_ringsv.axes: # go over all the subplots in the figure fig
        for i in ax.spines.itervalues(): # ... and go over all the axes too...
            i.set_linewidth(2)
        ax.minorticks_on() # switch on the minor ticks
        # set the tick lengths and tick widths
        ax.tick_params('both', length=5, width=2, which='major')
        ax.tick_params('both', length=3, width=1, which='minor')
    '''

    # p3v.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))

    # switch off plot box for rings
    rplot.set_axis_off()  # ChMade

    # one unit on the x-axis is of same length as one unit on the y-axis with the following:
    rplot.set_aspect('equal')

    (x1, x2, y1, y2) = rplot.axis()

    # new_y_rings:
    rplot.axis((-days_side, days_side, y_min_ring, y_max_ring))  # size of the ring plot gets set here   #ivanew

    # vertically zoom the lower plot centred on the splined fit
    # (x1, x2, y1, y2) = p2v.axis()
    # p2v.axis((x1, x2, 0., 1.1))

    # (x1, x2, y1, y2) = p3v.axis()
    # p3v.axis((x1, x2, -0.19, 0.19))

    majorLocator = MultipleLocator(0.1)
    # p3v.yaxis.set_major_locator(majorLocator)
    # p3v.axhline(y=0., color='k', ls='dashed')

    # p3v.set_xlabel('[days]')
    # p3v.set_ylabel('Error')
    # p2v.set_ylabel('Normalized Intensity')

    # fig_ringsv.savefig(plotout)   #ivanew

    print('\nDrawing figure done')
    print('---------------------\n')

