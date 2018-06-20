''' bring-ex - a set of routines used throughout my work on my
minor thesis at Leiden University, including everything I need
to do the BRING simulations
'''

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.io import ascii
import matplotlib as mpl
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from scipy.ndimage import convolve

import sys, getopt, os

sys.path.append('/Users/kenworthy/Dropbox/python_workbooks/lib')

# import pyfits
import exorings
import j1407

from scipy.optimize import fmin
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from matplotlib.patches import PathPatch
from scipy import stats

from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.gridspec import GridSpec

# constants taken from disk_sim.py
G = 6.673840e-11  # [m^3 kg^-1 s^-2] -> gravitational constant
yr = 365.242189669 * 86400  # [sec] -> year
msol = 1.98855e30  # [kg] -> solar mass
rsol = 6.5500e8  # [m] -> solar mass
mjup = 1.8986e27  # [kg] -> Jupiter mass
rjup = 6.9911e7  # [m] -> Jupiter radius
au = 1.49597870700e11  # [m] -> astronomical unit
pc = 3.0856e16  # [m] -> parsec

# beta Pic b parameters
a_Picb = 9.2 * au  # [m]
r_Picb = 1.65 * rjup  # [m]
m_Picb = 11 * mjup  # [kg]
M_beta_Pic = 1.61 * msol  # [kg] - this is the mass of the host star beta Pic, which is very, very uncertain
d_beta_Pic = 19.3 * pc  # distance to beta Pic [m]; from SIMBAD
R_beta_Pic = 4.77e8  # radius of beta Pic [m]; from SIMBAD


################
# Calculations #
################

def hill(a, m, M):
    """
    Calculates the Hill radius of a body
    From: https://en.wikipedia.org/wiki/Hill_sphere
    :param a: semi-major axis
    :param m: planet mass
    :param M: star mass
    :return: e
    """
    r = a * np.power(m / (3 * M), (float(1) / 3))
    return r


def vcirc(m1, m2, a):  # a = semi-major axis, m1,m2 = object masses
    """
    Calculates circular orbital velocity of m2 about m1 at distance a.
    From: http://en.wikipedia.org/wiki/Circular_orbit
    :param m1: object mass 1 [kg]
    :param m2: object mass 2 [kg]
    :param a: semi-major axis
    :return: vcirc [m/s]
    """
    mu = G * (m1 + m2)
    vcirc = np.power((mu / a), 0.5)
    return vcirc


def transit_time(vel, size_r):  # vel = velocity, size_r = transiting structure radius
    """
    Calculates the time it takes for a structire to transit the host star beta Pic.
    :param vel: velocity
    :param size_r: transiting structure radius
    :return: t [s]
    """
    # t_front_edge = 2*(R_beta_Pic/vel)
    # t_the object passing its own width = 2*(size_r/vel)
    # t = t_front_edge + t_object_width
    t = (2 / vel) * (R_beta_Pic + size_r)
    return t


def Ptoa(P, m1, m2):  # input not in SI
    """
    Calculates orbital radius from period
    :param P: period [yr]
    :param m1: mass of primary [Msol]
    :param m2: mass of secondary [Mjup]
    :return: semi-major axis [AU]
    """
    # a^3/P^2 = (G / 4 pi pi) (m1 + m2)
    c = G / (4. * np.pi * np.pi)
    mu = (m1 * msol) + (m2 * mjup)
    a3 = np.power(P * yr, 2.) * (c * mu)

    return np.power(a3, 1. / 3.) / au


def cross_section_to_angle(dist, r):
    """ calculates the projected cross-section as function of angle [degrees] from star, see Minor thesis fig. 9
	dist is distance between star and planet, theta is angle between planetary orbit and line of sight
    r is in meters
    """
    theta = np.arcsin(r / dist) * 180. / np.pi  # [degrees]
    return theta


def angle_to_cross_section(dist, theta):
    """ inverts cross_section_to_angle() """
    r = dist * np.sin(np.pi * theta / 180.)  # [m]
    return r


##############
# Data input #
##############

def moon_input(txt_file,
               table):  # not really using this anywhere - also because I changed the data form several times along the way
    """ Import supermoon data from .txt file generated elsewhere.
    Importing ID, mass, semimajor-axis and Hill radius."""
    super_table = Table.read(txt_file, format='ascii')
    super_names = np.array(["Sup-Io", "Sup-Europe", "Sup-Ganymede", "Sup-Callisto"])

    # defining arrays out of the supermoon table
    indices = np.array(super_table['ID'])
    masses = np.array(super_table['mass'])  # [kg]
    axes = np.array(super_table['axis'])  # [m]
    hill_radii = np.array(super_table['R_Hill'])  # [m]

    return super_names, indices, masses, axes, hill_radii


###########################################
# Import and interpolation of light curve #
###########################################

### Importing light curve data
def lightcurve_import(conv_curve, a, b, c):
    tcenter = 0  # centerpoint
    trange = 1500.  # width of plot window - just make it a lot bigger than the actual time range in data, so that it captures everything

    # range of days that ring statistics should be considered
    rings_tmin = tcenter - (trange / 2.)
    rings_tmax = tcenter + (trange / 2.)

    # print ''
    # import light curve data from .dat created by bring_disk_sim_v2.py
    light_table = Table.read(conv_curve, format='ascii')

    # defining arrays out of the light curve table
    time = np.array(light_table[a])
    flux = np.array(light_table[b])
    flux_errxx = np.array(light_table[c])
    # print ''
    # print 'time.size: ',time.size
    # print 'time:',time

    return time, flux, flux_errxx


### interpolation for 5min-sampling
def interpol_5min(time, flux):
    # this makes the function fluxcurve() that you can use to do interpolation with
    fluxcurve = interp1d(time, flux, kind='linear')

    # new_t is the new, higher time sampling you want
    # print ''
    new_time = np.arange(
        500. * 288.) / 288.  # np.arange(#of days * #of intervals per day)/#of intervals per day; #of days has to be small enough to be in the range of the data. If there is a wish to have a larger range, then the input data file has to be adjusted by bring_disk_sim.py
    new_time = new_time - (np.max(new_time) / 2.)  # centering
    # print 'new_time.size:',new_time.size
    # print 'new_time:', new_time

    new_flux = fluxcurve(new_time)  # the interpolation
    # print 'new_flux:',new_flux.size

    return new_time, new_flux


### introduce gaussian scatter
def gaussian_scatter(flux):
    mu, sigma = 0.0, 1.0
    gauss_noise = np.random.normal(mu, sigma, size=flux.shape)
    flux_with_noise = flux + 0.005 * gauss_noise

    return flux_with_noise


#########################################
# Forward scattering and phase function #
#########################################

def deg_to_days(dist, theta, v):
    length_in_days = angle_to_cross_section(dist, theta) / (24. * 60. * 60. * v)  # [days]

    return length_in_days


def days_to_deg(dist, day_input, v):
    invert = day_input * (24. * 60. * 60. * v)
    length_in_deg = cross_section_to_angle(dist, invert)

    return length_in_deg


def phasefunc(g, theta):
    """
    Calculate Henyey-Greenstein function.
    Source: https://www.astro.umd.edu/~jph/HG_note.pdf
    :param g: g parameter
    :param theta: angle
    :return:
    """
    phfct = (1. / 4. * np.pi) * (1. - np.sqrt(g)) / (1 + np.sqrt(g) - 2. * np.sqrt(g) * np.cos(theta)) ** (3. / 2.)

    return phfct


def make_phsfct_kernel(size_px, dpx, g_fac):
    """
    Make a kernel for phase function convolution
    :param size_px:
    :param dpx: [deg/px]
    :param g_fac:
    :return: ph_ker [deg]
    """
    ke = np.mgrid[:size_px, :size_px]
    half = (size_px - 1) / 2
    ke[0] -= half
    ke[1] -= half
    dist = np.sqrt(ke[0] * ke[0] + ke[1] * ke[1])
    dist_deg = dist * dpx
    ph_ker = phasefunc(g_fac, dist_deg)  # Fill radially with phase function
    # ph_ker = ph_ker/np.sum(ph_ker)
    ph_ker = ph_ker / (2. * np.pi)

    return ph_ker


def forward_scattering(i, phi, r, tau, y, theta, g_fac, v, star, width):  ### BROKEN ONE! ###
    """ includes forward scattering in a light curve
    theta in [degrees]
    v in [m/s]"""
    xrang = 800.  # [days] the length of the scattering-strip in days, make sure it's the same like for the convolution of the light curve in exorings.ellipse_strip_3()
    inter = angle_to_cross_section(a_Picb, theta)  # [m] the extent to which the phase function should be used
    yrang = np.max(inter / (
            24. * 60. * 60. * v))  # [days] the extent to which the phase function should be used                                             # (hours*minutes*seconds*velocity)

    star_x, star_y = np.shape(star)  # [px] size of star

    # calculate the effective pixel width and height for mgrid in days
    dt = width / (star_y - 1)  # [days/px] 1px entspricht dt Tagen
    nx = xrang / dt  # [px] length of ring strip
    ny = (
                 yrang / dt) + 1  # [px] width of ring strip (I added one because I want to make it an uneven number.) I do not know if this will hold up in terms of a working code later on. It will if neither the size of star in days or the star kernel size in pixels are changed. Or the orbital velocity.)

    yc = (ny - 1.) / 2.  # middle point of height of ring-strip
    xc = (nx - 1.) / 2.  # middle point of length of ring-strip

    agr = np.mgrid[:ny, :nx]  # THE RING-STRIP!!

    agr[0] -= yc  # centering it
    agr[1] -= xc  # centering it
    agr = agr * dt  # converting ring-strip to units of days

    # move up to the strip by the imported value y (= impact parameter b)
    agr[0] += y

    tau_disk, grad = exorings.ellipse_nest(agr, i, phi, r,
                                           tau)  # creates the values to fill the rings-strip (tau_disk) and the gradient of the rings, meaning the angles at which the star crosses the ring broders

    # test ring-strip
    # print "tau disk shape:",tau_disk.shape
    # print tau_disk[200,2000]
    # plt.imshow(tau_disk)
    # plt.show()

    # tau_disk has an outside region with transmission 0 which is there only because I defined the strip is bigger than the physical size of the disk. Here, I first set that right and set it to transmission 1.
    array_cp = np.asarray(tau_disk)
    zeros = (array_cp == 0)
    array_cp[zeros] = tau[1]
    ground_strip = array_cp

    # test adjusted ring-strip
    # print "value in center region:",ground_strip[500,5000]
    # print "value in outside region and gaps:",ground_strip[200,2000]
    # print "new strip shape:",ground_strip.shape
    # plt.imshow(ground_strip)
    # plt.show()

    # now I set all gaps to value 0 and all dust to value 1 to be able to overlay the phase function
    array_cp2 = np.asarray(ground_strip)
    be_one = array_cp2 < 0.95  # I chose 0.95 as limit because ~0.94 is where the dust absorbs
    be_zero = array_cp2 > 0.95  # and ~0.99 is where the gaps are, with full transmittance
    array_cp2[be_one] = 1  # 1 where the dust is, 0 where we have empty space (thus no ligth
    array_cp2[be_zero] = 0  # will be scattered from those parts)
    ring_strip = array_cp2

    # test final ring-strip
    # print "value in center region:",ring_strip[500,5000]
    # print "value in outside region and gaps:",ring_strip[200,2000]
    # print "ring strip shape:",ring_strip.shape
    # plt.imshow(ring_strip)
    # plt.show()

    # --- 2D solution start --- #

    # Generate phase function kernel
    ph_kernel = make_phsfct_kernel(ny, days_to_deg(a_Picb, dt, v), g_fac)

    # Test kernel
    # plt.imshow(ph_kernel)
    # plt.colorbar()
    # plt.show()
    kx, ky = np.shape(ph_kernel)
    # plt.plot(ph_kernel[(kx/2),:])
    # plt.show()

    # Convolve phasefunction kernel with ring strip
    ring_convolved = convolve(ring_strip, ph_kernel, mode='constant', cval=0.0)
    print('Shape of phase function convolved light curve is:'), ring_convolved.shape
    # plt.imshow(ring_convolved)
    # plt.show()

    # Sum the light up that gets into line of sight
    scatter_ph = np.sum(ring_convolved, axis=0)
    # plt.plot(scatter_ph)
    # plt.show()

    # --- 2D solution end --- #

    # --- 1D solution start --- #
    """
    # create array to hold the phase function
    single = ring_strip[:,:1]        # one pixel wide strip
    double = ring_strip[:,:2]        # two pixel wide strip - this is the one I use
    #single[-1,-1] = 2               # for testing
    #double[-1,-1] = 2               # for testing

    #print "double shape is:",double.shape

    # fill first of the pixel rows with equally spaced angle values
    npx_angle = 2*np.max(theta)/(double[:,0].shape[0])                      # degrees per pixel - used nowhere
    deg_in_pix = np.linspace(0.,2*np.max(theta),(double[:,0].shape[0]))     # !!! array that is as wide as ring strip but in units of degrees
    deg_in_pix = deg_in_pix - (np.max(deg_in_pix)/2.)                       # centering the angle array

    # generate phase function array
    phase_single = phasefunc(g_fac, deg_in_pix)
    double[:,0] = deg_in_pix                                 # one dimension holds degrees
    double[:,1] = phase_single                               # other dimension holds phase function values
    phase = double                                           # just renaming
    
    # test phase function array
    #print "phase shape is:",phase.shape
    #plt.imshow(phase, interpolation='none')
    #plt.show()

    #plt.plot(phase[:,0],phase[:,1])
    #plt.show()
    
    # testing the phase function array
    #print np.min(phase[:,0])
    #print np.max(phase[:,0])
    #print np.where(phase[:,1] == np.max(phase[:,1]))

    # renaming stuff makes it easier to use
    #################################
    # phase[:,0] - degrees
    # phase[:,1] - scattered light
    degrees = phase[:,0]
    scatter = phase[:,1]
    #################################
    
    #print phase[:,0].shape

    #testcol = ring_strip[:,4550:4551]         # picking a random pixel column in which structure is distinguishable
    #plt.imshow(testcol, interpolation='none') # not using this anywhere
    #plt.show()

    grow = np.resize(scatter,(ring_strip.shape[1],ring_strip.shape[0])) # I copy the single row array many times so that it fills an array of the same size like the ring strip
    phase_mask = np.transpose(grow)                                     # For some reasone the result is flipped, so I flip it back
    #plt.imshow(phase_mask, interpolation='none')
    #plt.show()

    #print phase_mask.shape
    #print ring_strip.shape

    overlay = phase_mask*ring_strip    # Using the 0 and 1 ring strip as a mask to determine which parts of phase function to take into account.
    #plt.imshow(overlay, interpolation='none')
    #plt.colorbar()
    #plt.show()

    light = np.sum(overlay, axis=0)     # Summing over y axis and collapsing light onto line of sight
    #print light.shape                   # Check that I summed over the correct axis
    #plt.plot(light)
    #plt.title("Scattered light")
    #plt.show()
    
    # --- 1D solution end --- #
    """
    return scatter_ph


def forward_scattering2(i, phi, r, tau, y, theta, g_fac, v, star, width):  ### working version ###
    """ includes forward scattering in a light curve
    theta in [degrees]
    v in [m/s]"""
    xrang = 800.  # [days] the length of the scattering-strip in days, make sure it's the same like for the convolution of the light curve in exorings.ellipse_strip_3()
    inter = angle_to_cross_section(a_Picb, theta)  # [m] the extent to which the phase function should be used
    yrang = np.max(inter / (
            24. * 60. * 60. * v))  # [days] the extent to which the phase function should be used                                             # (hours*minutes*seconds*velocity)

    star_x, star_y = np.shape(star)  # [px] size of star

    # calculate the effective pixel width and height for mgrid in days
    dt = width / (star_y - 1)  # [days/px] 1px entspricht dt Tagen
    nx = xrang / dt  # [px] length of ring strip
    ny = (
                 yrang / dt) + 1  # [px] width of ring strip (I added one because I want to make it an uneven number.) I do not know if this will hold up in terms of a working code later on. It will if neither the size of star in days or the star kernel size in pixels are changed. Or the orbital velocity.)

    yc = (ny - 1.) / 2.  # middle point of height of ring-strip
    xc = (nx - 1.) / 2.  # middle point of length of ring-strip

    agr = np.mgrid[:ny, :nx]  # THE RING-STRIP!!

    agr[0] -= yc  # centering it
    agr[1] -= xc  # centering it
    agr = agr * dt  # converting ring-strip to units of days

    # move up to the strip by the imported value y (= impact parameter b)
    agr[0] += y

    tau_disk, grad = exorings.ellipse_nest(agr, i, phi, r,
                                           tau)  # creates the values to fill the rings-strip (tau_disk) and the gradient of the rings, meaning the angles at which the star crosses the ring broders

    # test ring-strip
    # print "tau disk shape:",tau_disk.shape
    # print tau_disk[200,2000]
    # plt.imshow(tau_disk)
    # plt.show()

    # tau_disk has an outside region with transmission 0 which is there only because I defined the strip is bigger than the physical size of the disk. Here, I first set that right and set it to transmission 1.
    array_cp = np.asarray(tau_disk)
    zeros = (array_cp == 0)
    array_cp[zeros] = tau[1]
    ground_strip = array_cp

    # test adjusted ring-strip
    # print "value in center region:",ground_strip[500,5000]
    # print "value in outside region and gaps:",ground_strip[200,2000]
    # print "new strip shape:",ground_strip.shape
    # plt.imshow(ground_strip)
    # plt.show()

    # now I set all gaps to value 0 and all dust to value 1 to be able to overlay the phase function
    array_cp2 = np.asarray(ground_strip)
    be_one = array_cp2 < 0.95  # I chose 0.95 as limit because ~0.94 is where the dust absorbs
    be_zero = array_cp2 > 0.95  # and ~0.99 is where the gaps are, with full transmittance
    array_cp2[be_one] = 1  # 1 where the dust is, 0 where we have empty space (thus no ligth
    array_cp2[be_zero] = 0  # will be scattered from those parts)
    ring_strip = array_cp2

    # test final ring-strip
    # print "value in center region:",ring_strip[500,5000]
    # print "value in outside region and gaps:",ring_strip[200,2000]
    # print "ring strip shape:",ring_strip.shape
    # plt.imshow(ring_strip)
    # plt.show()

    # --- 2D solution start --- #

    # Generate phase function kernel
    ph_kernel = make_phsfct_kernel(ny, days_to_deg(a_Picb, dt, v), g_fac)
    kx, ky = np.shape(ph_kernel)

    # Test kernel
    # plt.imshow(ph_kernel)
    # plt.colorbar()
    # plt.show()
    # plt.plot(ph_kernel[(kx/2),:])
    # plt.show()

    # Convolve phasefunction kernel with ring strip
    ring_convolved = convolve(exorings.y_to_tau(ring_strip), ph_kernel, mode='constant', cval=0.0)
    print('Shape of phase function convolved light curve is:'), ring_convolved.shape
    # plt.imshow(ring_convolved)
    # plt.show()

    # Sum the light up that gets into line of sight
    scatter_ph = np.sum(ring_convolved, axis=0)
    # plt.plot(scatter_ph)
    # plt.show()

    # --- 2D solution end --- #

    x_tdc = agr[1, int(yc), :]  # time array
    return x_tdc, scatter_ph


#########################################
# Binning and plotting the light curves #
#########################################

### bin to sampling of 1 hour
def binned_1h(new_time, flux_noise):
    run = 1  # naming the plots
    fit = True  # do you want the fit in the plots?
    days_side = 200  # days to either side of zero point to be shown
    total = new_time.size

    bin_means_1h, bin_edges_1h, binnumber_1h = stats.binned_statistic(
        new_time[int((total / 2.) - (days_side * 288.)):int((total / 2.) + (days_side * 288.))],
        flux_noise[int((total / 2.) - (days_side * 288.)):int((total / 2.) + (days_side * 288.))], statistic='median',
        bins=days_side * 2 * 24)  # do the binning

    # print ''
    # print 'bin_means 1h sampling:',bin_means_1h
    # print 'bin means size 1h sampling:', bin_means_1h.size
    # print 'bin_edges 1h sampling:',bin_edges_1h
    # print 'bin edges size 1h sampling:', bin_edges_1h.size
    # print 'binnumber 1h sampling:',binnumber_1h
    # print 'binnumber size 1h sampling:', binnumber_1h.size

    # finding bin midpoints for 1h sampling
    binsize_1h = bin_edges_1h[1] - bin_edges_1h[0]  # figuring out the size of the bins (all are equal)
    bin_edges_1h_miss = np.delete(bin_edges_1h, bin_edges_1h[-1])  # deleting the last entry in bin_edges
    bin_time_1h = np.zeros_like(bin_means_1h)  # create array to hold the time values for bin_means_1h
    bin_time_1h = bin_edges_1h_miss + (binsize_1h / 2.)  # center the values of bin_time in the middle of the bins

    # print ''
    # print 'bin time size 1h sampling:',bin_time_1h.size
    # print 'bin means size 1h sampling:', bin_means_1h.size

    # get the standard deviation of the bins in 1hr sampling
    bin_err_1h = np.zeros_like(bin_means_1h)  # create array for errors, same length as bin mean values
    for mid in range(bin_means_1h.size):
        left = bin_edges_1h[mid]  # left bin edge
        right = bin_edges_1h[mid + 1]  # right bin edge
        indices = np.where((new_time > left) & (new_time < right))  # find the time indices of the current bin interval
        interval = flux_noise[indices]  # get the flux values of the current bin interval
        std = np.std(interval)  # std from data of the current bin interval
        bin_err_1h[mid] = std
        # print 'flux:',interval
        # print 'std:',std
    # print 'bin_err_1h.size:',bin_err_1h.size

    return bin_time_1h, bin_means_1h


def plot_1h_sampling(ax, bin_time_1h, bin_means_1h, new_time, new_flux, days_side):
    # print ''
    # print 'PLOTTING 1h SAMPLED DATA'
    # print ''

    # plotting 1 hour sampling
    ax.scatter(bin_time_1h, bin_means_1h, s=1)
    # plt.errorbar(bin_time_1h, bin_means_1h, yerr=bin_err_1h, fmt='k.')   # errorbars still very, very dense
    ax.plot(new_time, new_flux, c='r', linewidth=1.0)  # overplot the ideal curve
    # ax.axis((-days_side, days_side, 0.92, 1.02))              # original limits
    ax.axis((-days_side, days_side, 1.6, 1.8))
    # ax.axis((-days_side, days_side, 0.5, 10.))
    # ax.axis((-days_side, days_side, 60., 100.))
    # ax.set_xlabel('Time [Days]')
    # ax.set_ylabel('Transmission')
    # ax.set_title('set a title')    # describes the physical parameters used in bring_disk_sim.py to create the rings system which light curve we are presenting here
    # ax.tick_params(width=3, length=15)
    # ax.yticks(fontsize=40)
    # ax.xticks(fontsize=40)
    # plt.legend()


#########################################
# Stellar pulsations #
#########################################

### include stellar pulsations of all amplitudes and frequencies in light curve flux
def stellar_puls(time, flux, amplitude, frequency):
    """
    Generates a sine signal to overlay with the light curve.
    All of the input parameters are arays, even if A and f can be filled with a single entry.
    :param time:
    :param flux:
    :param amplitude:
    :param frequency:
    :return:
    """

    t = np.linspace(np.min(time), np.max(time), flux.size)  # Create array equivalent in size and min and max
    # on x-axis like the input flux.
    phi = np.pi / 2.  # Define a phase (can be anything, it just shifts
    # the sine wave on the x-axis).
    sine = np.zeros((t.size, amplitude.size))  # Create an empty array to hold all waves from all frequencies.


    for i in range(amplitude.size):  # loop over all frequencies and amplitudes
        sine[:, i] = amplitude[i] * np.sin(2 * np.pi * frequency[i] * t + phi)  # create each individual sine wave

    sine_tot = np.sum(sine, axis=1)  # sum up all the different sine waves
    flux_with_pulse = flux + sine_tot  # add the total sine wave to the input flux

    return flux_with_pulse, sine_tot


def cd_to_uHz(freq_cd):
    """
    Transforms the input frequency from cycles per day to microHerz.
    :param freq_cd:
    :return:
    """
    freq_uHz = freq_cd * 11.574074  # [microHz]

    return freq_uHz


def cd_to_periodDay(freq_cd):
    """
    Transforms the input frequency from cycles per day to a period in units of days.
    :param freq_cd:
    :return:
    """
    period_day = 1 / freq_cd  # [days]

    return period_day


###################################################
### Grids and axes for the light curves and rings #
###################################################

def name_axes(fig):
    """
    Puts axis names into the different parts of the grids. For testing and checking.
    :param fig:
    :return: none
    """
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i + 1), va="center", ha="center")


def make_x_ticklabels_invisible(fig, x_list):
    """
    Deletes the x-axis labels of the axes in the list "no_x_list".
    :param fig:
    :param x_list:
    :return: none
    """
    for i, ax in enumerate((x_list)):
        for tl in ax.get_xticklabels():
            tl.set_visible(False)


def make_y_ticklabels_invisible(fig, y_list):
    """
    Deletes the y-axis labels of the axes in the list "no_y_list".
    :param fig:
    :param y_list:
    :return: none
    """
    for i, ax in enumerate((y_list)):
        for tl in ax.get_yticklabels():
            tl.set_visible(False)


def make_x_labels(fig, axes):
    for i, ax in enumerate(axes):
        ax.set_xlabel('Days', fontsize=20)


def make_y_labels(fig, axes):
    for i, ax in enumerate(axes):
        ax.set_ylabel('Transmission', fontsize=20)


def x_tick_dist(fig, axes, interval):
    for i, ax in enumerate(axes):
        loc = mpl.ticker.MultipleLocator(base=interval)
        ax.xaxis.set_major_locator(loc)


def set_axes_limits(fig, xmin, xmax, ymin, ymax, ax_list):
    """
    Sets the limits on x- and y-axis for all the axes in ax_list.
    :param fig:
    :param xmin:
    :param xmax:
    :param ymin:
    :param ymax:
    :param ax_list:
    :return: none
    """
    for i, ax in enumerate((ax_list)):  # light curve limits
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
