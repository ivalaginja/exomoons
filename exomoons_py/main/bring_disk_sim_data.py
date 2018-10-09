###########################
#    Iva Laginja, 2016    #
###########################

# updated value for semi-major axis of beta Pic b!
"""
This program should create a ring system with tau values and subsequently create a light curve for different disk
parameters (tilt etc.) that can be put in by hand. In this program you sett the physical parameters of the system: moon
masses and distances, planetary and host star parameters and so on (full list further below). This then gets used by
'data_creation.py' to produce ring systems with different tips and tilts and then 'grid_lightcurve_plots.py' creates a
grid with ring systems and light curves.

Hashtags:
#setval - here values can be chaged maually
#pyadjust - here I adjusted the original code by commenting out big parts or somethig the like
#testing - testing parts
-Chekpoint - output for checkpoints

TO CHANGE NUMBER OF MOONS or MOON PROPERTIES or RING PROPERTIES: #moon-no
(also check these when something is not like you want it to be, chances are that you messed up somewhere in these)

- adjust main parameters in config_local.ini
- adjust ring radius values at #setvalradii
- adjust tau values


----------------------------------------------------------------
This program is not run on its own, it is only
USED by grid_lightcurve_plots.py (need to import it there)
----------------------------------------------------------------

Parts:

#-------------------------------------------------------------------------
# PART 1: Generating the supermoons and calculating planetary parameters |
#-------------------------------------------------------------------------
# PART 2: Generating ring parameters: radii and tau values               |
#-------------------------------------------------------------------------
# PART 3: Generating the light curve of the tilted system                |
#-------------------------------------------------------------------------

"""

import sys, getopt, os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.io import ascii
# from scipy.ndimage import convolve
from scipy.interpolate import interp1d

import exomoons_py.modules.exorings as exorings
import exomoons_py.modules.j1407 as j1407
import exomoons_py.modules.bring as bring
from exomoons_py.config import CONFIG_INI


# mpl.interactive(True)
# set sensible imshow defaults
mpl.rc('image', interpolation='nearest', origin='lower', cmap='gray')
mpl.rc('axes.formatter', limits=(-7, 7))

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


# switch - implemented from http://code.activestate.com/recipes/410692/

class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args:  # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False


def ringfunc(taun, *args):
    """Cost function for ring fit to photometric data."""
    (t, f, f_err, rad_ri, re, k, dst) = args
    # convolve and make smoothed ring photometry
    strip, dummy, g = exorings.ellipse_strip3(rad_ri, \
                                         exorings.y_to_tau(taun), re[0], re[1], re[2], re[3], k, dst)

    # interpolate the smoothed curve....
    ring_model = interp1d(g[0], g[1], kind='linear')

    # ... to the times of the photometry
    ring_model_phot = ring_model(t)

    # calculate the residuals and the chi squared
    diff = f - ring_model_phot
    chisq = np.sum(np.power(diff / f_err, 2))
    red_chisq = chisq / diff.size

    return red_chisq


def calc_ring_stats(taun, t, f, f_err, rad_ri, re, k, dst, tmin, tmax):
    """Full statistics function for ring fit to photometric data."""
    # convolve and make smoothed ring photometry
    strip, dummy, g = exorings.ellipse_strip(rad_ri, \
                                        exorings.y_to_tau(taun), re[0], re[1], re[2], re[3], k, dst)

    # select points within the rings_tmin/tmax range
    mask = (t > tmin) * (t < tmax)
    t_sel = t[mask]
    f_sel = f[mask]
    f_err_sel = f_err[mask]

    print('%d points in time range %.2f to %.2f' % (t_sel.size, tmin, tmax))

    # interpolate the smoothed curve....
    ring_model = interp1d(g[0], g[1], kind='linear')

    # ... to the times of the photometry
    ring_model_phot = ring_model(t_sel)

    # calculate the residuals and the chi squared
    diff = f_sel - ring_model_phot

    chisq = np.sum(np.power(diff / f_err_sel, 2))
    # degrees of freedom = number of photometry points - number of ring edges - 1
    dof = diff.size - taun.size - 1
    red_chisq = chisq / dof
    print('number of photometric = %d ' % diff.size)
    print('number of ring edges  = %d ' % taun.size)
    print('number of DOF         = %d ' % dof)
    print('chi squared           = %.2f' % chisq)
    print(' reduced chi squared  = %.2f' % red_chisq)

    # http://en.wikipedia.org/wiki/Bayesian_information_criterion
    # n - number of points in data
    # k - number of free parameters
    # BIC = chisquared + k . ln(n) + C
    # C is a constant which does not change between candidate models but is
    # dependent on the data points

    BIC = chisq + (taun.size) * np.log(diff.size)
    print(' BIC                  = %.2f' % BIC)
    return red_chisq


nn = 1


def costfunc(x, *args):
    (y, dt, i_deg, phi_deg) = x
    (grad_t, grad_mag_n, t0) = args
    # grad_time
    global nn

    (tmp, grad_disk_fit) = exorings.ring_grad_line(grad_t, y, dt, i_deg, phi_deg)

    # lazy way of calculating the time midpoint of the light curve
    rmintime = np.arange(np.min(grad_t), np.max(grad_t), 0.01)
    (tmp, rminline) = exorings.ring_grad_line(rmintime, y, dt, i_deg, phi_deg)
    rmint = rmintime[np.argmin(rminline)]

    # make a cost function
    delta = grad_disk_fit - grad_mag_n

    # if delta is positive, keep it
    # if delta is negative, make it positive and multiply by 50
    delta[np.where(delta < 0)] = -delta[np.where(delta < 0)] * 50.

    # dean is penalty to clamp rmint
    dean = np.abs(rmint - t0)

    cost = np.sum(delta) + (dean * 20)
    nn += 1
    return cost


def ind_ring(ring, r):
    """Returns index for closest ring r."""
    rdiff = ring - r
    # find index of smallest positive rdiff
    return np.argmin(np.abs(rdiff))


def ind_ring_big(ring, r):
    """Returns index for closest bigger ring r."""
    rdiff = ring - r
    rdiff[(rdiff < 0)] += 99999.
    # find index of smallest positive rdiff
    return np.argmin(rdiff)


print('')
print('-Checkpoint: BEGIN main program.')
print('')


########################################################################################
# BEGIN main program (I made it a single function, because it's used in another program)
########################################################################################
def bring_disk_data(impact, i_in, phi_in, output, targetdir, modelnumber, paramfile=False):

    # stellar and planetary parameters of beta Pic b
    a_Picb = CONFIG_INI.getfloat('beta_Pic', 'planet_axis') * au  # [m]  semi-major axis of planetary orbit around host star
    r_Picb = CONFIG_INI.getfloat('beta_Pic', 'planet_radius') * rjup  # [m]  planet radius
    m_Picb = CONFIG_INI.getfloat('beta_Pic', 'planet_mass') * mjup  # [kg] planet mass
    M_beta_Pic = CONFIG_INI.getfloat('beta_Pic', 'star_mass') * msol  # [kg] mass of the host star beta Pic, which is very, very uncertain
    d_beta_Pic = CONFIG_INI.getfloat('beta_Pic', 'star_distance') * pc  # [m]  distance to beta Pic; from SIMBAD
    R_beta_Pic = CONFIG_INI.getfloat('beta_Pic', 'star_radius') * rsol  # [m]  radius of beta Pic; from paper [http://iopscience.iop.org/article/10.1086/509912/pdf]

    # -----------------------------------------
    # PART 1: Generating the supermoons and calculating planetary parameters

    # Read in from config file which moon system you want to use.
    moon_system = CONFIG_INI.get('moon_parameters', 'system')
    moon_number = CONFIG_INI.getint(moon_system, 'number_of_moons')

    # Create empty arrays of size according to number of moons to hold masses and semi-major axes.
    moon_names = np.chararray(moon_number)
    index_array = np.zeros(moon_number)
    masses_array = np.zeros(moon_number)
    axes_array = np.zeros(moon_number)

    # Fill moon parameter arrays.
    for i in range(moon_number):
        moon_names[i] = CONFIG_INI.get(moon_system, 'name' + str(i+1))
        index_array[i] = i
        masses_array[i] = CONFIG_INI.getfloat(moon_system, 'mass' + str(i+1))
        axes_array[i] = CONFIG_INI.get(moon_system, 'axis' + str(i+1))

    # Blowing the system up by either scaling it up to beta Pic's measures or random scaling factors or both.
    scale_fac_a = CONFIG_INI.getfloat('moon_parameters', 'scale_fac_axes')
    scale_fac_m = CONFIG_INI.getfloat('moon_parameters', 'scale_fac_masses')
    scale_m_random = CONFIG_INI.get('moon_parameters', 'scale_masses_random')
    scale_a_random = CONFIG_INI.get('moon_parameters', 'scale_axes_random')
    scale_m_betaPic = CONFIG_INI.get('moon_parameters', 'scale_masses_betaPic')
    scale_a_betaPic = CONFIG_INI.get('moon_parameters', 'scale_axes_betaPic')

    # Convert masses to Jupiter masses and axes to Jupiter radii.
    m_sat_mjup = masses_array / mjup  # [mjup]!!!    # converting to Jupiter masses
    a_sat_rjup = axes_array / rjup  # [rjup]!!!      # converting to Jupiter radii [rjup]

    if scale_m_betaPic:
        m_sat_mjup = m_sat_mjup * m_Picb           # [kg] array filled with masses of the supermoons
    if scale_a_betaPic:
        a_sat_rjup = a_sat_rjup * r_Picb               # [m]  array filled with the sem.-maj. axes of the supermoons
    if scale_m_random:
        m_sat_mjup = m_sat_mjup * scale_fac_m      # [kg]
    if scale_a_random:
        a_sat_rjup = a_sat_rjup * scale_fac_a          # [m]

    moon_masses = m_sat_mjup
    moon_axes = a_sat_rjup

    # Calculate Hill radii of b Pic b and its supermoons.
    hill_betaPicb = bring.hill(a_Picb, m_Picb, M_beta_Pic)  # [m] Hill radius of beta Pic b (single number)
    hill_sat_Picb = bring.hill(moon_axes, moon_masses, m_Picb)  # [m] Hill radii of the supermoons (array)

    print('')
    print('------------------------------------------------')
    print('Mass and distance parameters of the moons')
    print('beta Pic b Hill sphere [m]:', hill_betaPicb)
    print('Moon masses in terms of m_Picb: ', moon_masses / m_Picb)
    print('Moon masses in terms of mearth:', moon_masses / mearth)
    print('Moon distances in terms of r_Picb:', moon_axes / r_Picb)
    print('Moon distances in terms of planet Hill sphere:', moon_axes / hill_betaPicb)
    print('')
    print('mearth/m_Picb:', mearth / m_Picb)
    print('mearth/mearth:', mearth / mearth)
    print('------------------------------------------------')

    # creating a .txt file of these parameters, for further use
    out0 = moon_names  # names of supermoons
    out1 = index_array  # ID's of supermoons
    out2 = moon_masses  # masses of supermoons [kg]
    out2c = moon_masses / mjup  # masses of supermoons [mjup]
    out2a = moon_masses / m_Picb  # masses of supermoons [m_Picb]
    out2b = moon_masses / mearth  # masses of supermoons [mearth]
    out3 = moon_axes  # semi-major axes of supermoons [m]
    out3a = moon_axes / hill_betaPicb  # semi-major axes of supermoons [R_H planet]
    out3b = moon_axes / r_Picb  # semi-major axes of supermoons [r_Picb]
    out4 = hill_sat_Picb  # Hill radii of supermoons [m]
    out4a = hill_sat_Picb / r_Picb  # Hill radii of supermoons [r_Picb]

    super_table = Table([out0, out1, out2, out2c, out2a, out2b, out3, out3a, out3b, out4, out4a],
                        names=['name', 'ID', 'mass[kg]', 'mass[mjup]', 'mass[m_Picb]', 'mass[mearth]', 'axis[m]',
                               'axis[R_H planet]', 'axis[r_Picb]', 'R_Hill[m]', 'R_Hill[r_Picb]'],
                        meta={'name': 'supermoons'})
    super_table.write(os.path.join(targetdir, output + '.txt'), format='ascii')

    print('')
    print("Writing supermoons' IDs, masses, sem.-maj. axes and Hill radii into 'model.txt'.")
    print('')

    # TEST print stuff   #testing
    print('masses of Galilean moons:', masses_array)
    print('axes of Galilean moons:', axes_array)
    print('')
    print('masses of Galilean moons in terms of mjup:', m_sat_mjup)
    print('axes of Galilean moons in terms of rjup:', a_sat_rjup)
    print('')
    print('masses of supermoons:', moon_masses)
    print('axes of supermoons:', moon_axes)
    print('')
    print('Hill radius of beta Pic b in SI unit meter:', hill_betaPicb)
    print('Hill radius of beta Pic b in terms of its radius:', hill_betaPicb / r_Picb)
    print('Hill radii of supermoons in terms of Pic bs radius:', hill_sat_Picb / r_Picb)
    print('')
    print('-Checkpoint: Creating moons and exporting the data worked.')
    print('')

    # calculate the orbital velocity of beta Pic b
    v_orb_calc = bring.vcirc(M_beta_Pic, m_Picb, a_Picb)  # [m/s]
    ###################################################
    # Decide which orbital velocity to use:            #setval
    # v_orb = v_orb_calc      # use calculated velocity
    v_orb = CONFIG_INI.getfloat('beta_Pic', 'planet_vcirc')  # use user value
    ###################################################

    # calculate transit times of beta Pic b's and the supermoons' Hill spheres
    t_Hill_Picb = bring.transit_time(v_orb, hill_betaPicb)  # [s]
    t_Hill_moons = bring.transit_time(v_orb, hill_sat_Picb)  # [s]

    # calculate star size in days
    d_star = (2 * R_beta_Pic / v_orb) / (60 * 60 * 24)  # diameter of star [days]

    # calculating total duration of planet's hill eclipse in days
    transit_days = (np.around(t_Hill_Picb) / 86400. + 1.)

    # -----------------------------------------
    # PART 2: Generating ring parameters: radii and tau values

    # creating a time array in days, as long as the Hill sphere transit lasts
    day_no = transit_days + 300.  # how many days the time arry should hold. here: extended to be longer than the entire Hill sphere transit
    days_array_start = np.linspace(0., day_no,
                                   day_no * 288.)  # [days], in 5min intervals (there are 288 5min intervals in a day)
    # empty_days = np.linspace(np.max(days_array_start),np.max(days_array_start)+50.,50*288) # extending the time array beyond the Hill sphere
    # days_array_start2 = np.append(days_array_start,empty_days)
    days_array = days_array_start - (day_no / 2.)  # centering it on 0

    # define center point of array, which will be where the array entry is 0
    ecl_cen = (days_array.size - 1.) / 2.

    # defining intervals of full rings and gaps over inner and outer edges of empty rings cleared out by the supermoons
    inner_edges = np.zeros_like(index_array)  # inner edges of gaps
    outer_edges = np.zeros_like(index_array)  # outer edges of gaps

    # creating correct values for inner and outer edges
    for i in range(len(index_array)):
        inner_edges[i] = (moon_axes[i] - hill_sat_Picb[i]) / v_orb  # [s]
        outer_edges[i] = (moon_axes[i] + hill_sat_Picb[i]) / v_orb  # [s]

    # transforming radii to units of days
    inner_edges = inner_edges / 86400.  # [days]
    outer_edges = outer_edges / 86400.  # [days]

    # create a list to hold both inner and outer radii
    radii = []
    for i in range(len(index_array)):
        radii.append(inner_edges[i])
        radii.append(outer_edges[i])

    Hill_param = "Half filled"
    radii.append(
        transit_days / 4.)  # adding an outermost radius in order to cut off the dust filled part far outside the moon gaps. This determines to what extent the Hill sphere is filled. If the argument is (transit_days/2.), then the full Hill sphere is filled. If it is (transit_days/4.), then only half of the Hill sphere (until 1/2 R_H) is filled.

    ######################################################
    ###### manual adjustment of ring radii #setvalradii

    # radii = [4., 5.4, 6.9, 7.4, 8., 8.7, 12.1, 13., 13.9, 14.8, 63.]   # Steven dyn sim 1 values # The last value is technically infinity, because that is outside of the (filled) Hill sphere.

    # radii = [23.7331829, 24.06185223, 24.39052157, 24.7191906, 37.3975677, 38.15855437, 38.91954103, 39.6805305, 59.8393059, 60.93222923, 62.02515257, 63.1180801, 78.0403366]   # values for manually blown up b=17% system with 3 supergalilean moons, with manually 1/3-gap dust lanes in the cleared gaps

    # radii = [24.0045147, 24.15229606666667, 24.300077433333332, 24.4478588, 38.0257874, 38.3679632, 38.710139, 39.0523148, 60.7415504, 61.23297753333333, 61.724404666666665, 62.2158318, 78.0403366]   # values for manually blown up b=17% system with 3 galilean moons, with manually 1/3-gap dust lanes in the cleared gaps

    ########## ORIGINAL ##############
    radii = [23.84712982, 24.09983444, 24.35253906, 24.60524368, 37.6613884, 38.246496836666665, 38.83160527333334, 39.41671371, 60.21819305, 61.058521266666666, 61.898849483333336, 62.7391777, 78.04033661] # values for manually blown up b=17% system with 3 galilean*5 moons, with manually 1/3-gap dust lanes in the cleared gaps
    ########## ORIGINAL END ##########

    # radii for solution (i) to fit Lecavelier light curve - straight thirds
    #radii = [60.21819305, 61.058521266666666, 61.898849483333336, 62.7391777]

    # radii for adapting the sizes to trying to fit curve by adjusting scattering rings strip sizes
    # outer strips are narrorwer by half and middle strip bigger by third (I think)
    # radii = [60.6383571583, 61.058521266666666, 62.3190135917, 62.7391777]

    # another try
    #radii = [61.06, 61.058521266666666, 61.898849483333336, 61.9]

    # radii = [23.66184044, 24.038071316666667, 24.41430219333333, 24.79053307, 37.23237228, 38.103491463333334, 38.97461064666667, 39.84572983, 59.60205841, 60.853146869999996, 62.10423533, 63.35532379, 78.04033661]   # values for manually blown up b=17% system with 3 supergalilean*1.5 moons, with manually 1/3-gap dust lanes in the cleared gaps

    # radii = [23.80214119, 24.084838233333333, 24.36753527666667, 24.65023232, 37.55720901, 38.21176528666667, 38.86632156333333, 39.52087784, 60.06858826, 61.00865428, 61.9487203, 62.88878632, 78.04033661]   # values for manually blown up b=17% system with 3 galilean*7 moons, with manually 1/3-gap dust lanes in the cleared gaps

    radii.append(300.)  # This value is technically infinity, because that is outside of the (filled) Hill sphere.

    ###### manual adjustment end
    ######################################################

    # create initial array for tau values, where gaps totally free and rings totally opaque
    taun_rings = np.resize([5.2, 100.],
                           len(radii))  # needs to have as many entries as radii list has, np.resize does that

    print('')
    print('- Checkpoint: radii and taun_rings same length/size?')
    print('radii length (list):', len(radii))
    print('taun_rings shape (array):', taun_rings.shape)
    print('tau_rings values:', taun_rings)
    print('')

    # taus for solution (i) to fit Lecavelier light curve
    #taun_rings = [100., 5.2, 3.0, 5.2, 100.]

    # ---
    # creating light curve plot for face-on situation
    t_gap = 1
    t_ring = 0
    I_0 = 1

    I_array = np.zeros_like(days_array)

    for n in range(len(index_array)):
        gap = np.where((np.abs(days_array) >= (inner_edges[n])) & (np.abs(days_array) <= (outer_edges[n])))
        print("Gap at moon #", n, "goes")
        print("from", inner_edges[n], "to", outer_edges[n])
        print(".")
        I_array[gap] = I_0 * t_gap

    # finding the area of interest, where the moons are in the light curve and zooming in there
    # and then plotting the light curve
    '''#COMMENT
    interest = np.where(I_array != 0)
    low_lim = np.min(interest)
    up_lim = np.max(interest)
    
    # extending the plotted part a bit
    plot_from = low_lim - 500.
    plot_to = up_lim + 500.
    
    plt.plot(days_array[plot_from:plot_to],I_array[plot_from:plot_to])
    
    plt.ylim(ymax = 1.1, ymin = 0.)
    plt.xlabel("Transit time in relation to the transit center (=0) [days]")
    plt.ylabel("Normalized intensity")
    plt.axvline(0,0,1.1,color='r',dashes=[1,1])   # marking the center of the transit
    #plt.show()
    '''
    # ---

    ### ---
    # data output for the other program that makes the drawings

    # create tables to export the photometry data
    flux_rms = np.ones_like(I_array)  # I need a placeholder for this column, so I set them all to one value
    flux_rms = flux_rms * 0.005

    light = Table([days_array, I_array, flux_rms], names=('time', 'flux', 'flux_rms'))

    eclipse_start = days_array[0]
    eclipse_end = days_array[-1]
    eclipse_center = days_array[
        int(ecl_cen)]  # it should be 0 here, but it's not exactly zero. Not sure if the result is close enough though.
    eclipse_duration = transit_days
    col_1 = Column([eclipse_start])
    col_2 = Column([eclipse_end])
    col_3 = Column([eclipse_center])
    col_4 = Column([eclipse_duration])

    time_points = Table([col_1, col_2, col_3, col_4], names=('start', 'end', 'center', 'duration'))

    # rad_radii = Table([radii_ar], names=('rad_radii'))  # not working at the moment

    # save the two tables to files
    ascii.write(light, os.path.join(targetdir, 'lightcurve_bring.dat'))
    ascii.write(time_points, os.path.join(targetdir, 'time_params_bring.dat'))

    ### ---

    # define disk parameters for initial, edge-on, ring system
    # DON'T CHANGE THESE!!! EVER!!! They're the starting point for the whole simulation.
    res = [0, 0, 0, 0]
    res[0] = 0  # impact parameter b [days] - setting it to 0 makes the star go through center of planet
    res[1] = ecl_cen  # date of central eclipse t_b = dt
    res[2] = 0  # inclination angle i_deg [degrees]: 0 means face-on, as seen from earth
    res[
        3] = 90  # phi_deg [degrees]: angle between the normal of the secondary companion's orbit and the normal of the ring plane

    # write radii, taus, disk parameters and star size into .fits file, to be used further
    exorings.write_ring_fits(os.path.join(targetdir, "face-on_disk.fits"), res, taun_rings, radii, d_star)
    print('')
    print("Writing ring radii, ring taus, disk parameters and star size into 'face-on_disk.fits'.")
    print('')

    # TEST print stuff again   #testing
    print('orbital velocity (calculated) is [m/s]:', v_orb_calc)
    print('orbita velocity (used):', v_orb)
    print("time for planet's Hill sphere to transit, in days:", transit_days)
    print("supermoons' Hill spheres transit in days:", t_Hill_moons / 86400.)
    print('star size in days:', d_star)
    print('radii in days:', radii)
    print('')
    print('-Checkpoint: Second round of input and prints works.')
    print('')

    # -----------------------------------------
    # PART 3: Generating the light curve of the tilted system

    # Here is where I copied my latest version of the disk_sim.py code and then adjusted it.
    # I'll now just comment out everything that has anything to do with importing a light curve - because I actually want to generate it, not to fit it.
    # COMMENT
    # pyadjust
    # following few lines can be ignored; things getting imported simply to make the program work
    tcenter = 0  # time of wiggle in curve - centerpoint

    ####################################setval
    trange = 400.  # width of photometry (and "useless" plotting window appearing)
    ####################################

    ring_offset = 0.  # days offset for the disk center

    # range of days that ring statistics should be considered
    rings_tmin = tcenter - (trange / 2.)
    rings_tmax = tcenter + (trange / 2.)

    (time, flux, flux_err) = j1407.j1407_photom_binned(os.path.join(targetdir, 'lightcurve_bring.dat'), rings_tmin, rings_tmax)

    # plot the data
    # plt.scatter(time, flux)

    print('restricting statistics of KIC to HJD range %.1f to %.1f' % (rings_tmin, rings_tmax))
    goodp_rings = (time > rings_tmin) * (time < rings_tmax)
    good_rings_npoints = goodp_rings.size
    print('number of points for statistics of KIC is %d' % (good_rings_npoints))
    # COMMENT
    # get gradients
    (grad_time, grad_mag, grad_mag_norm) = j1407.j1407_gradients('../../input_data/j1407_gradients.txt')

    # ignoring until here
    # ---

    # parse input options

    fitsin = 'ring002.fits'
    fitsout = 'ring003.fits'

    read_in_ring_parameters = False
    read_in_disk_parameters = False

    vstar = -1

    def print_help():
        print('bring_disk_sim.py -o <outputfile>')

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:r:o:t:s:",
                                   ["dfile=", "rfile=", "ofile=", "tx=", "vstar="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in ("-d", "--dfile"):
            fitsin_disk = arg
            read_in_disk_parameters = True
        elif opt in ("-r", "--rfile"):
            fitsin_ring = arg
            read_in_ring_parameters = True
        elif opt in ("-o", "--ofile"):
            fitsout = arg
        elif opt in ("-t", "--tx"):
            tx = np.array(float(arg))
            print('tx = time of central eclipse forced to be = ', tx)
        elif opt in ("-s", "--vstar"):
            vstar = np.array(float(arg))
    print('Output file is ', fitsout)

    # I won't be reading in a .fits file, but take the data from first part of program   #pyadjust
    fitsin_disk = os.path.join(targetdir, "face-on_disk.fits")  # this is for the array res, which we are changing further below anyway
    fitsin_ring = os.path.join(targetdir, "face-on_disk.fits")
    # fitsout = tilted_disk.fits        # name of the new .fits file with parameters for tilted ring system
    vstar = v_orb

    ###################################################
    # read in or create the ring system tau and radii
    # setval
    phi_deg = phi_in  # tilt disk in degrees
    i_deg = i_in  # inclination in degrees (0 = face-on, 90 = edge-on; use 1 insted of 0 for fancy math reasons)   # i_deg = 63.3 deg is Saturn's tilt, because in these programs obliquity = 90-i for some reason
    dt = 0.  # date of central eclipse (dr/dt)=0
    y = impact * (transit_days / 2.)  # impact parameter b (units of days)
    ###################################################

    # putting the disk parameters into res, which is used later on, because they're not read into program
    res = [y, dt, i_deg, phi_deg]

    # stellar and planetary parameters: (left in here to not have to change variable names)
    ##
    Rstar = R_beta_Pic / rsol  # [rsol] star radius - the division throuh rsol is here, because R_beta_Pic is given in [m]
    Mstar = M_beta_Pic / msol  # [msol] star mass - same as above
    Mb = m_Picb / mjup  # [mjup] planet mass - same as above
    Pb = 5.5  # [yr] planet period (not needed anywhere)

    # length of eclipse in days
    # t_ecl = trange      #pyadjust
    t_ecl = transit_days

    # convert to an orbital velocity
    a = bring.Ptoa(Pb, Mstar, Mb)  # [AU]; input not in SI is ok
    print('Primary mass     = %5.2f  msol' % Mstar)
    print('Primary radius   = %5.2f  rsol' % Rstar)
    print('Secondary mass   = %5.2f  mjup' % Mb)
    print('Orbital radius   = %5.2f  AU' % a)
    v = bring.vcirc(Mstar * msol, Mb * mjup, a * au)  # [m/s]; input needs to be in SI

    if vstar > 0:
        v = vstar
        print('manual velocity of star is %.1f km.s-1' % v)

    print('Orbital velocity = %5.2f  km/s (use option -s to set new velocity) - option not available in this version' % (
            v / 1.e3))

    # dstar = (Rstar * rsol * 2 / v) / 86400.   #pyadjust
    dstar = d_star  # from above

    print('Primary diameter = %5.2f  days' % dstar)

    # pyadjust
    ## no fitting to be done, I want to keep my above defined radii and tau values - force it to read in the .fits file, even when it didn't get called in the terminal
    # reading in radii and tau's from the .fits file input when program got called in the terminal
    # if read_in_ring_parameters:
    print('')
    print('Reading in rings from %s' % fitsin_ring)
    (resxx, taun_rings, rad_rings, xxxdstar) = exorings.read_ring_fits(fitsin_ring)

    # else:
    #    print "Starting with new rings...."
    #    rad_rings = np.array([59.0])
    #    taun_rings = np.array([0.0])
    #    rad_rings = np.append(rad_rings, (100.))
    #    taun_rings = np.append(taun_rings, (1000.))

    ## now it always reads in the .fits file generated in PART 2

    # printing list of radii and taus
    exorings.print_ring_tau(rad_rings, exorings.y_to_tau(taun_rings))

    print('\n CURRENTLY MAKING MODEL', modelnumber)

    # set up stellar disk
    kern = exorings.make_star_limbd(21, 0.8)
    # plt.imshow(kern)
    # plt.show()

    # ---
    ## ignoring from here... #pyadjust
    tcenter = 0  # pyadjust

    # make the radius and projected gradient for the measured gradient points
    (ring_disk_fit, grad_disk_fit) = exorings.ring_grad_line(grad_time, res[0], res[1], res[2], res[3])

    # produce fine grained gradient and ring values
    samp_t = np.arange(-50, 50, 0.001) + tcenter
    (samp_r, samp_g) = exorings.ring_grad_line(samp_t, res[0], res[1], res[2], res[3])
    hjd_minr = samp_t[np.argmin(samp_g)]

    hjd_to_ring = interp1d(samp_t, samp_r, kind='linear')

    (rstart, rend) = exorings.ring_mask_no_photometry(kern, dstar, time, res[0], res[1], res[2], res[3])

    ## ... to here
    # --- but I need hjd_minr, so I won't make a comment out of it

    ### RESULTS of fitting routine
    # keeping it in for display purposes

    print('')
    print('Disk parameters fitting to gradients')
    print('------------------------------------')
    print('')
    print(' impact parameter b   = %8.2f days' % res[0])
    print(' HJD min approach t_b = %8.2f days' % res[1])
    print(' disk inclination i   = %7.1f  deg' % res[2])
    print('        disk tilt phi = %7.1f  deg' % res[3])
    print(' HJD min gradient     = %8.2f days' % hjd_minr)
    print('             rmin     = %8.2f days' % np.min(samp_r))

    # testing
    print('')
    print('-Checkpoint: Ignoring parts of the program works.')
    print('')

    # http://en.wikipedia.org/wiki/Bayesian_information_criterion
    # n - number of points in data
    # k - number of free parameters
    # BIC = chisquared + k . ln(n) + C
    # C is a constant which does not change between candidate models but is
    # dependent on the data points

    ###########################################
    # pyadjust
    # Don't need all of these plots, as I don't intend to plot and/or fit light curve, but make one instead.

    # pyadjust
    # but I need some of the variables, so I'll leave it in the program now

    # plot folded light curve
    time0 = np.abs(time - hjd_minr)
    time0_grad = np.abs(grad_time - hjd_minr)

    # flux_color and flux_col
    # hold the color of the points for ingress and egress
    flux_color = np.chararray(time.shape, unicode=True)
    flux_color[:] = 'b'
    flux_color[(time > hjd_minr)] = 'r'

    # probably a better pythonic way to do this, but this works.
    flux_col = ''
    for b in flux_color.tolist():
        flux_col = str.join('', (flux_col, b))

    def plot_folded_phot(f):
        """Plot folded J1407 light curve."""

        # j1407 photometry
        h1.scatter(time0, flux, c=flux_col, s=20, edgecolors='none', zorder=-20)
        h1.errorbar(time0, flux, flux_err, zorder=-30, ls='none')

        # gradient measurements
        # h1.scatter(time0_grad,np.ones_like(time0_grad)*0.8)

    fig_fold = plt.figure(figsize=(16, 6))

    h1 = fig_fold.add_subplot(111)
    # plot_folded_phot(fig_fold)

    ###########################################

    # creating the light curve
    strip, dummy, g = exorings.ellipse_strip_3(rad_rings, exorings.y_to_tau(taun_rings), \
                                          res[0], res[1], res[2], res[3], kern, dstar)
    # strip and dummy never used again in this program - just ignore them
    print('')
    print('####################################################')
    print(days_array.size)

    # normalization by the part where there is nothing in front of the star, because transmission is normamlized to 1 there
    norm_place = np.where((g[0] < -50.))
    norm_set = np.max(norm_place)
    norm_val = np.max(g[1][norm_set])
    g[1] = g[1] / norm_val

    # the values I need are stored in the g array, which are:
    # g[0] = time (x-axis)
    # g[1] = stellar convolved tau - cenetral pixel line after convolution -> light curve
    # g[2] = stellar tau without convolution
    # g[3] = gradients

    ### save the fitted/convolved light curve to .dat file
    # creating fake errors
    g_err = np.ones_like(g[1])
    g_err = g_err * 0.005

    # export
    conv_light = Table([g[0], g[1], g_err], names=(
    'time', 'flux', 'flux_rms'))  # flux_rms still there cause I'll need it at the readout in bring_lightcurves.py
    ascii.write(conv_light, os.path.join(targetdir, output + '.dat'))

    ###########################################

    # creating the scattered light
    g_fac = 0.99
    theta_lim = 0.2  # limiting angle to either side of phase function; if you need a wider ring strip, for testing purposes for example, make this value bigger, e.g. to 3.
    step_size = 0.01  # [degrees]
    theta = np.linspace(0., 2 * theta_lim,
                        (2 * theta_lim / step_size) + 1)  # +1 because I need an uneven number of entries
    theta = theta - (np.max(theta) / 2.)  # centering the theta values on zero

    sc_time, sc_light = bring.forward_scattering2(res[2], res[3], rad_rings, exorings.y_to_tau(taun_rings), res[0], theta,
                                                  g_fac, v_orb, kern, dstar)

    # normalizing the scattered light is tricky. here's a trial
    scat_norm = 1.  # peak value of P(theta) = P(0) I want to normalize to
    sc_light = (sc_light / np.max(sc_light)) * scat_norm

    # plt.plot(sc_time,sc_light)
    # plt.title('in bring_disk_sim_data.py')
    # plt.show()

    # save scattered light to file, separate form convolved light curve
    scat_light = Table([sc_time, sc_light, g_err], names=('sc_deg', 'sc_light', 'pseudo_error'))
    ascii.write(scat_light, os.path.join(targetdir, output + '_scat.dat'))

    # save summed up light of both scattered and occulted light curve
    tot_light = Table([sc_time, g[1] + sc_light, g_err], names=('tot_deg', 'tot_light', 'pseudo_error'))
    ascii.write(tot_light, os.path.join(targetdir, output + '_tot.dat'))

    print('')
    print('-Checkpoint: Exporting the fitted/convolved curve works.')
    print('')
    print("Moon's Hill spheres in days:", outer_edges - inner_edges)
    print("Moon's Hill spheres in hours:", (outer_edges - inner_edges) * 24.)

    # rearanging the g parameters for the light curve fitting and plotting (in interactive mode)
    gt_abs = np.abs(g[0] - hjd_minr)
    g1 = g[1]
    gt_ingr = gt_abs[(g[0] <= hjd_minr)]
    gt_egr = gt_abs[(g[0] > hjd_minr)]
    g1_ingr = g1[(g[0] <= hjd_minr)]
    g1_egr = g1[(g[0] > hjd_minr)]
    '''#COMMENT
    # plotting the light curve
    plt.plot(np.abs(g[0]-hjd_minr), g[1])                  # total curve
    plt.plot(gt_ingr, g1_ingr, color='blue')               # ingress
    plt.plot(gt_egr, g1_egr, color='red')                  # egress
    plt.plot(np.abs(g[0]-hjd_minr), g[2], color='orange')  # unconvolved tau - would be like perfect fac-on
    plt.xlabel('Time from eclipse midpoint [days]')
    plt.ylabel('Transmission')
    
    plt.show() # - shows both the face-on light curve as well as the convolved plotting window
    '''  # COMMENT

    exorings.write_ring_fits(os.path.join(targetdir, output + '.fits'), res, taun_rings, rad_rings, dstar)

    # Creating a parameter text file in targetdir to know whats happening
    if paramfile == True:
        params = open(os.path.join(targetdir, "system_parameters.txt"), "w")
        params.write("Amount of moons: %d moons \n" % int(index_array.shape[0]))
        params.write("Moon masses in kg: %s times factor %s \n" % (str(masses_array), scale_fac_m))
        params.write("Moon distances in m: %s times factor %s \n" % (str(axes_array), scale_fac_a))
        params.write("Moon Hill sphere transit times in sec: %s \n" % (str(t_Hill_moons)))
        params.write("Hill sphere is: %s \n" % Hill_param)
        params.write("Hill sphere transit of beta Pic b in days: %s \n" % transit_days)
        params.write("\n")
        params.write("Radii of lines of different transmission:\n")
        params.write("%s \n" % str(radii))
        params.write("tau values:")
        params.write("%s \n" % str(taun_rings))
        params.write("\n")
        params.write("------------------------------------------------")
        params.write("\n")
        params.close()
