''' exorings - a set of routines to calculate ring geometries
 ring_to_sky() and sky_to_ring() are the fundamental transformations
 from the sky to the ring plane

 coordinate system is right-handed
 with rings originally in the X-Z plane centered at the origin

 Earth observers are along +ve Z-axis looking back at the origin

 rings are rotated i degrees about X axis
 then rotated phi degrees about z axis so that
 +ve phi is from +ve x axis towards +ve Y axis (i.e. CCW)

       Y
       ^
       |
       |
       +-----> X
      /
     /
    Z

 (Xr,Yr) are coordinates in the ring (X,Z) plane
 (Xs,Ys) are coordinates in the sky  (X,Y) plane
'''

import numpy as np
from astropy.io import fits
import matplotlib as mpl
from matplotlib import patches as pth
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from scipy.ndimage import convolve


def ring_to_sky(Xr, Yr, i_deg, phi_deg):
    """
    Go from ring coordinate frame to the sky coordinate frame.
    :param Xr:
    :param Yr:
    :param i_deg:
    :param phi_deg:
    :return:
    """

    i = np.pi * i_deg / 180.
    phi = np.pi * phi_deg / 180.

    Xs = (np.cos(phi) * Xr) + (np.sin(phi) * np.sin(i) * Yr)
    Ys = (np.sin(phi) * Xr) - (np.cos(phi) * np.sin(i) * Yr)
    Rs = np.sqrt((Xs * Xs) + (Ys * Ys))

    return Xs, Ys, Rs


def sky_to_ring(Xs, Ys, i_deg, phi_deg):
    """
    Take tilted ellipse geometry and get Xr,Yr,Rr in ring coordinates.
    :param Xs:
    :param Ys:
    :param i_deg:
    :param phi_deg:
    :return:
    """

    i = np.pi * i_deg / 180.
    phi = np.pi * phi_deg / 180.

    Xr = ((np.cos(phi) * Xs) + (np.sin(phi) * Ys))
    Yr = ((np.sin(phi) * Xs) - (np.cos(phi) * Ys)) / np.sin(i)
    Rr = np.sqrt((Xr * Xr) + (Yr * Yr))

    return Xr, Yr, Rr


def ring_tilt(i_deg, phi_deg):
    """
    Given the phi and i of the ring plane, return the total tilt.
    :param i_deg:
    :param phi_deg:
    :return: tilt
    """

    i_rad = np.pi * i_deg / 180.
    phi_rad = np.pi * phi_deg / 180.

    # convert the two angles to a total tilt
    ct = np.cos(i_rad) * np.cos(phi_rad)
    tilt = np.arccos(ct) * 180. / np.pi

    return tilt


def ellipse(im, i_deg, phi_deg):
    """
    Make an ellipse in an input array.
    :param im: output from mgrid with [0] containing y values [1] x values
    :param i_deg: inclination of rings in degrees (0 = face on)
    :param phi_deg: rotation in CCW direction from x-axis of rings in degrees
    :return: image of r and angle of tangent to local ellipse
    """

    Ys = im[0]
    Xs = im[1]
    # x,y are points in the sky plane that we transform to the ring plane
    # ellipse projection will occur automatically
    (x11, y11, ellipser) = sky_to_ring(Xs, Ys, 90. - i_deg, phi_deg)

    return ellipser


def ellipse_nest(im, i_deg, phi_deg, r=None, tau=None):
    """
    Make an ellipse in an input array.

    ASSUMES that r is ordered in increasing r!!!
    Source: http://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
    :param im: output from mgrid with [0] containing y values [1] x values
    :param i_deg: inclination of rings in degrees (0 = face on)
    :param phi_deg: rotation in CCW direction from x-axis of rings in degrees
    :param r: optional - vector with the radii of major axis of ellipses
    :param tau: optional - value to be put in ring r
    :return: radius r from the ring plane in the sky plane coordinates; tangent to the ellipse in radians
    """

    Ys = im[0]
    Xs = im[1]

    # im[0],im[1] are points in the sky plane that we transform to the ring plane
    # ellipse projection will occur automatically
    (Xr, Yr, ellipse) = sky_to_ring(Xs, Ys, 90. - i_deg, phi_deg)

    # in ring space, the tangent of a circle on point (X,Y) is (Y,-X)
    (Xs_tang, Ys_tang, ellipse2) = ring_to_sky(Yr, -Xr, 90. - i_deg, phi_deg)

    tan_out = np.arctan2(Ys_tang, Xs_tang)

    im_out = ellipse

    if r is not None:
        im_out = np.zeros_like(Xs) + 0.0
        r_inner = 0.0
        # loop from smallest r to biggest r
        for r_outer, curr_tau in zip(r, tau):
            sel = (ellipse >= r_inner) * (ellipse < r_outer)
            im_out[np.where(sel)] = curr_tau

            # must be last line in r_outer loop
            r_inner = r_outer

    return im_out, tan_out


def ellipse_para(a, b, phi, t):
    """
    Create a parametric ellipse.

    Source: http://en.wikipedia.org/wiki/Ellipse
    :param a: major axis
    :param b: minor axis
    :param phi: angle between the x and semimajor axis of the ellipse in a +ve phi = CCW about x axis
    :param t: parametric parameter in radians (0 to 2pi)
    :return: (X, Y, tang) in cartesian coordinates and tang is tangent to the ellipse
    """

    ct = np.cos(t)
    st = np.sin(t)
    cp = np.cos(phi)
    sp = np.sin(phi)

    # X and Y positions
    X = (a * ct * cp) - (b * st * sp)
    Y = (a * ct * sp) + (b * st * cp)

    # gradient at X,Y
    dY = -a * st * sp + b * ct * cp
    dX = -a * st * cp - b * ct * sp

    # grad = np.arctan2(dX,dY)
    tang = np.arctan2(dY, dX)

    return X, Y, tang


def ellipse_tangents(a, b, phi):
    """
    Create tangents to the ellipse.
    :param a: semimajor axis
    :param b: semiminor axis
    :param phi: +ve x axis CW to semimajor axis
    :return: dX0 - parametric parameter at which dX is zero, dy0 - parametric parameter at which dY is zero,
             txaxis - tangent angle along y=0 line for nonzero x
    """

    # cos i = b/a
    print('phi: ', phi)
    cosi = b / a
    iang = np.arccos(cosi)
    dX0 = np.arctan(-(cosi) * np.tan(phi))
    dY0 = np.arctan((cosi) / np.tan(phi))

    # what is dy/dx along y=0 for nonzero values of x?
    denom = np.sin(2 * phi) * np.sin(iang) * np.sin(iang)
    numer = 2 * (np.sin(phi) * np.sin(phi) + cosi * cosi * np.cos(phi) * np.cos(phi))

    tang = np.arctan(numer / denom)

    print('phi is %f rad' % phi)
    print('tang is %f rad' % tang)

    return dX0, dY0, tang


def ellipse_strip(r, tau, y, hjd_central, i, phi, star, width):
    """
    Calculate a strip of values suitable for star convolution.
    :param r:
    :param tau:
    :param y: y coordinate in days offset for star strip
    :param hjd_central:
    :param i:
    :param phi:
    :param star: 2D array with the image of the star
    :param width: the full width of the strip in days - the star's size in units of days
    :return:
    """

    xrang = 800.  # [days] the length of the ring-strip in days
    yrang = width  # [days] the height of the ring-strip

    star_x, star_y = np.shape(star)  # size of star

    ny = star_y
    # calculate the effective pixel width for mgrid in days
    dt = yrang / (ny - 1)  # [days/px] 1px entrpricht dt Tagen

    nx = np.floor(xrang / dt)

    yc = (ny - 1) / 2.  # middle point of height of ring-strip
    xc = (nx - 1) / 2.  # middle point of length of ring-strip

    agr = np.mgrid[:ny, :nx]  # THE RING-STRIP!!

    agr[0] -= yc  # centering it
    agr[1] -= xc  # centering it
    agr = agr * dt  # converting ring-strip to units of days

    # move up to the strip by the imported value y (impact parameter b)
    agr[0] += y

    tau_disk, grad = ellipse_nest(agr, i, phi, r,
                                  tau)  # creates the values to fill the rings-strip (tau_disk) and the gradient of the rings, meaning the angles at which the star crosses the ring broders

    tau_disk_convolved = convolve(tau_disk, star, mode='constant',
                                  cval=0.0)  # convolving the star with the ring-strip, general big convolution yielding a bigger strip

    # pick out central rows from convolution and the actual xvals - these form the "g" array used in main program
    x_tdc = agr[1, int(yc), :] + hjd_central  # time intervals; for x-axis
    tdc = tau_disk_convolved[int(yc), :]  # central pixel line - the convolved light curve (taus!!)
    td = tau_disk[int(yc), :]  # cetral pixel line of unconvolved ring strip
    g = grad[int(yc), :]  # gradients

    return tau_disk, tau_disk_convolved, np.vstack((x_tdc, tdc, td, g))


def ellipse_strip_3(r, tau, y, hjd_central, i, phi, star, width):   # edited duplicate
    """
    Calculate a strip of values suitable for star convolution. (Edited duplicate.)
    :param r:
    :param tau:
    :param y: y coordinate in days offset for star strip
    :param hjd_central:
    :param i:
    :param phi:
    :param star: 2D array with the image of the star
    :param width: the full width of the strip in days - the star's size in units of days
    :return:
    """
    xrang = 800.  # [days] the length of the ring-strip in days
    yrang = width  # [days] the height of the ring-strip

    star_x, star_y = np.shape(star)  # size of star

    ny = star_y
    # calculate the effective pixel width for mgrid in days
    dt = yrang / (ny - 1)  # [days/px] 1px entspricht dt Tagen

    print('')
    print('###########################################################')
    print('dt:', dt)
    nx = xrang / dt

    yc = (ny - 1) / 2.  # middle point of height of ring-strip
    xc = (nx - 1) / 2.  # middle point of length of ring-strip

    agr = np.mgrid[:ny, :nx]  # THE RING-STRIP!!

    agr[0] -= yc  # centering it
    agr[1] -= xc  # centering it
    agr = agr * dt  # converting ring-strip to units of days

    # move up to the strip by the imported value y (= impact parameter b)
    agr[0] += y

    tau_disk, grad = ellipse_nest(agr, i, phi, r,
                                  tau)  # creates the values to fill the rings-strip (tau_disk) and the gradient of the rings, meaning the angles at which the star crosses the ring broders

    tau_disk_convolved = convolve(tau_disk, star, mode='constant',
                                  cval=0.0)  # convolving the star with the ring-strip, general big convolution yielding a bigger strip

    # pick out central rows from convolution and the actual xvals - these form the "g" array used in main program
    x_tdc = agr[1, int(yc), :] + hjd_central  # time intervals; for x-axis
    # tdc = tau_disk_convolved[yc, :]       # central pixel line - the convolved light curve (taus!!)
    tdc = np.sum(tau_disk_convolved, axis=0)  # trying with sum because of scattered light trial
    td = tau_disk[int(yc), :]  # cetral pixel line of unconvolved ring strip
    g = grad[int(yc), :]  # gradients

    return (tau_disk, tau_disk_convolved, np.vstack((x_tdc, tdc, td, g)))


def ellipse_strip_2(r, tau, y, hjd_central, i, phi, star, width, transit_days):
    """ calculate a strip of values suitable for star convolution
        y - y coordinate in days offset for star strip
        star - 2D array with the image of the star
        width - the full width of the strip in days - the star's size in units of days
    """
    # xrang = time_size*2   # [days] the length of the ring-strip in pixels
    yrang = width  # [days] the height of the ring-strip

    star_x, star_y = np.shape(star)  # size of star

    ny = star_y
    # calculate the effective pixel width for mgrid in days
    dt = yrang / (ny - 1)  # 1px entrpricht dt Tagen
    min5 = 1. / 288.  # 5 min in Einheit von Tagen (1 Tag hat 288 5-min Intervalle)

    # nx = xrang/dt

    days_no = transit_days + 100.
    days_strip = np.linspace(0., days_no, days_no * 288.)
    # days_array = days_strip - (days_no/2.)
    nx = days_strip.size

    yc = (ny - 1) / 2.  # middle point of height of ring-strip
    xc = (nx - 1) / 2.  # middle point of length of ring-strip

    agr = np.mgrid[:ny, :nx]  # THE RING-STRIP!!

    agr[0] -= yc  # centering it
    agr[1] -= xc  # centering it
    agr[0] = agr[0] * dt  # converting ring-strip to units of days
    agr[1] = agr[1] * min5  # converting ring-strip to units of days

    # move up to the strip by the imported value y (impact parameter b)
    agr[0] += y

    tau_disk, grad = ellipse_nest(agr, i, phi, r,
                                  tau)  # creates the values to fill the rings-strip (tau_disk) and the gradient of the rings, meaning the angles at which the star crosses the ring broders

    tau_disk_convolved = convolve(tau_disk, star, mode='constant',
                                  cval=0.0)  # convolving the star with the ring-strip, general big convolution yielding a bigger strip

    # pick out central rows from convolution and the actual xvals - these form the "g" array used in main program
    x_tdc = agr[1, yc, :] + hjd_central  # time intervals; for x-axis
    tdc = tau_disk_convolved[yc, :]  # central pixel line - the convolved light curve (taus!!)
    td = tau_disk[yc, :]  # cetral pixel line of unconvolved ring strip
    g = grad[yc, :]  # gradients

    print('')
    print('######################################################################')
    print('nx: ', nx)
    print(agr[1])
    print(agr.shape)
    print('')

    return tau_disk, tau_disk_convolved, np.vstack((x_tdc, tdc, td, g))


def ring_grad_line(xt, yt, dt, i_deg, phi_deg):
    """
    Make a line of gradients using xt vector.
    :param xt:
    :param yt: vertical offset
    :param dt: xt offset for zero gradient
    :param i_deg:
    :param phi_deg:
    :return: (ellipse, grad_disk)
    """

    yy = np.ones_like(xt) * yt
    xx = xt - dt

    aa = np.vstack((yy, xx))

    (ellipse, gradient) = ellipse_nest(aa, i_deg, phi_deg)

    # convert gradient to abs projected along x-axis
    grad_disk = np.abs(np.sin(gradient))

    return ellipse, grad_disk


def ring_mask_no_photometry(star, width, xt, yt, dt, i_deg, phi_deg):
    """Find valid radii."""

    ring_radii = np.arange(0, 100, 0.001)
    ring_valid = np.zeros_like(ring_radii)

    star_x, star_y = np.shape(star)

    xc = (star_x - 1) / 2.
    yc = (star_y - 1) / 2.

    sgr = np.mgrid[:star_y, :star_x]

    tsize = width / (star_x - 1)

    # center the array in both x and y
    np.subtract(sgr[0], yc, out=sgr[0], casting="unsafe")  # This used to be "sgr[0] -= yc", but changed for new numpy
    np.subtract(sgr[1], xc, out=sgr[1], casting="unsafe")  # This used to be "sgr[1] -= xc", but changed for new numpy
    sgr = sgr * tsize  # transform to time

    masked = sgr[:, (star > 0)]
    # masked now is a flat array containing 2N elements
    # first N has the y positions
    # second N has the x positions

    for x in xt:
        # for each photometric time point in x, calculate the
        # max r and min r subtended by the star

        # pass those points to the routine to return r for all those points
        agr2 = np.copy(masked)
        agr2[0] += yt
        agr2[1] += (x - dt)

        rell = ellipse(agr2, i_deg, phi_deg)

        # get min r and max r in the area covered by the star at that
        # time
        minr = np.min(rell)
        maxr = np.max(rell)

        # mark that range of radii as valid
        ring_valid[(ring_radii > minr) * (ring_radii < maxr)] = 1.0

    # work up from zero r, find the boundaries for each bad masked ring
    drv = ring_valid[1:] - ring_valid[:-1]
    drv = np.append(drv,
                    [0.0])  # wild guess at how to fix this because of error with new numpy when this line is left out

    ringstart = ring_radii[(drv < -0.5)]
    ringends = ring_radii[(drv > 0.5)]

    # these two conditions make sure that the bad ring boundaries are
    # correctly offset
    # if the first ring_valid is zero. then the
    # pixels are bad from r=0 up to the first ring
    if ring_valid[0] == 0:
        ringstart = np.insert(ringstart, 0, 0.0)
    if ring_valid[-1] == 0:
        ringends = np.append(ringends, ring_radii[-1])

    return (ringstart, ringends)


def y_to_tau(y):
    return (np.arctan(y) / np.pi) + 0.5


def tau_to_y(tau):
    return np.tan(np.pi * (tau + 0.5))


########################
## DRAWING PRIMITIVES
########################

def ring_patch(r1, r2, i_deg, phi_deg, dr=([0, 0])):
    """
    Mmake a Patch in the shape of a tilted annulus.
    :param r1: inner radius
    :param r2: outer radius
    :param i_deg: inclination (degrees)
    :param phi_deg: tilt in CCW direction from x-axis (degrees)
    :param dr: (x,y) centre of the annulus
    :return:
    """

    i = np.cos(i_deg * np.pi / 180.)

    # get an Ellipse patch that has an ellipse
    # defined with eight CURVE4 Bezier curves
    # actual parameters are irrelevant - get_path()
    # returns only a normalised Bezier curve ellipse
    # which we then subsequently transform
    e1 = pth.Ellipse((1, 1), 1, 1, 0)

    # get the Path points for the ellipse (8 Bezier curves with 3
    # additional control points)
    e1p = e1.get_path()

    c = np.cos(phi_deg * np.pi / 180.)
    s = np.sin(phi_deg * np.pi / 180.)

    rotm = np.array([[c, s], [s, -c]])

    a1 = e1p.vertices * ([1., i])
    a2 = e1p.vertices * ([-1., i])

    e1r = np.dot(a1 * r2, rotm) + dr
    e2r = np.dot(a2 * r1, rotm) + dr

    new_verts = np.vstack((e1r, e2r))
    new_cmds = np.hstack((e1p.codes, e1p.codes))
    newp = Path(new_verts, new_cmds)

    return newp


##############################
## STELLAR FUNCTIONS
##############################

def make_star_limbd(dstar_pix, u=0.0):
    """Make a star disk with limb darkening."""
    ke = np.mgrid[:dstar_pix, :dstar_pix]
    dp2 = (dstar_pix - 1) / 2
    ke[0] -= dp2
    ke[1] -= dp2
    re = np.sqrt(ke[0] * ke[0] + ke[1] * ke[1])
    ren = re / dp2
    mask = np.zeros_like(ren)
    mask[(ren > 1.0000001)] = 1.
    ren[(ren > 1.0000001)] = 1.
    I = 1. - u * (1 - np.sqrt(1 - ren * ren))
    I[(mask > 0)] = 0.
    I /= np.sum(I)

    return I


###############################
## FILE INPUT AND OUTPUT
###############################


def write_ring_fits(fitsname, res, taun_rings, radii, dstar):
    """ Make Column objects with our output data."""
    col1 = fits.Column(name='taun', format='E', array=taun_rings)
    col2 = fits.Column(name='radius', format='E', array=radii)

    # create a ColDefs object for all the columns
    cols = fits.ColDefs([col1, col2])

    # create the binary table HDU object - a BinTableHDU
    tbhdu = fits.TableHDU.from_columns(cols)  # used to be fits.new_table(cols), but astropy changed for that

    prihdr = fits.Header()
    prihdr['TIMPACT'] = (res[0], 'Impact parameter (days)')
    prihdr['TMINR'] = (res[1], 'Time of minimum disk radius (days)')
    prihdr['DINCL'] = (res[2], 'Disk inclination (degrees)')
    prihdr['DTILT'] = (res[3], 'Disk tilt to orbital motion (degrees)')
    prihdr['DSTAR'] = (dstar, 'Diameter of star (days)')
    prihdr['HN'] = (0.907, 'Henweigh parameter')

    # open a PrimaryHDU object with no data (since you can't have TableHDU
    # in a PrimaryHDU) and append the TableHDU with the header
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(fitsname, clobber=True)
    print
    'write_ring_fits: wrote FITS file to %s' % fitsname


def read_ring_fits(fitsname):
    hdulist = fits.open(fitsname)
    prihdr = hdulist[0].header  # the primary HDU header

    # read in header keywords specifying the disk
    re = np.array((prihdr['TIMPACT'], \
                   prihdr['TMINR'], \
                   prihdr['DINCL'], \
                   prihdr['DTILT']))

    tbdata = hdulist[1].data
    return (re, tbdata['taun'], tbdata['radius'], prihdr['DSTAR'])


def print_ring_tau(rad, tau):
    """Pretty printing of ring radii and their tau values."""
    n = 0
    for (r, t) in zip(rad, tau):
        print
        'Ring %3d: tau = %5.3f out to radius %7.7f days' % (n, t, r)
        n += 1


def print_ring_tau_latex(rad, tau):
    """Pretty printing of ring radii and their tau values to a latex table."""
    n = 0
    from astropy.io import ascii
    for (r, t) in zip(rad, tau):
        print
        'Ring %3d: tau = %5.3f out to radius %7.3f days' % (n, t, r)
        n += 1
    from astropy.table import Table
    exptau = -np.log(tau)
    t = Table([rad, exptau], names=['Radius', 'Tau'])
    t['Radius'].format = '%.1f'
    t['Tau'].format = '%4.2f'
    ascii.write(t, output='ring_table1.tex', Writer=ascii.latex.AASTex, \
                col_align='ll', latexdict={'caption': \
                                               r'Table of ring parameters \label{tab:ring}', \
                                           'preamble': r'\tablewidth{0pt} \tabletypesize{\scriptsize}'})


def print_disk_parameters(res, minr_t, samp_r):
    print('')
    print('Disk parameters fitting to gradients')
    print('------------------------------------')
    print('')
    print(' impact parameter b   = %8.2f days' % res[0])
    print(' HJD min approach t_b = %8.2f days' % res[1])
    print(' disk inclination i   = %7.1f  deg' % res[2])
    print('        disk tilt phi = %7.1f  deg' % res[3])
    print(' HJD min gradient     = %8.2f days' % minr_t)
    print('             rmin     = %8.2f days' % np.min(samp_r))

    # make latex table with values
    ss = r'\begin{eqnarray*}b =& %8.2f \rm{d} \\ t_b =& %8.2f \rm{d} \\ i_{disk} =& %5.1f^o \\ \phi =& %5.1f^o \\ t_\parallel =& %8.2f \rm{d}\end{eqnarray*}' % (
    res[0], res[1], res[2], res[3], minr_t)

    return ss


def make_ring_grad_line(xt, yt, dt, i_deg, phi_deg):
    """
    Make a line of gradients using xt vector.
    :param xt:
    :param yt: vertical offset
    :param dt: xt offset for zero gradient
    :param i_deg:
    :param phi_deg:
    :return:
    """
    yy = np.ones_like(xt) * yt
    xx = xt - dt

    aa = np.vstack((yy, xx))

    (ellipse, gradient) = ellipse_nest(aa, i_deg, phi_deg)

    # convert gradient to abs projected along x-axis
    grad_disk = np.abs(np.sin(gradient))

    return ellipse, grad_disk


def draw_badrings(rstart, rend, xcen, incl, phi, p):
    for (bs, be) in zip(rstart, rend):
        path = ring_patch(bs, be, incl, phi, ([xcen, 0]))
        pnew = PathPatch(path, facecolor='#DDDDDD', edgecolor='none', zorder=-9)
        p.add_patch(pnew)


def draw_rings_vector(r, tau, xcen, incl, phi, p, ringcol='red'):
    ycen = 0.0
    xrang = 50.
    yrang = 50.

    p.set_xlim(xcen - xrang, xcen + xrang)
    p.set_ylim(ycen - yrang, ycen + yrang)

    for i in np.arange(0., r.size):
        if i == 0.:
            rin = 0.
            rout = r[0]
        else:
            rin = r[int(i - 1.)]
            rout = r[int(i)]

        ttau = 1. - tau[int(i)]

        path = ring_patch(rin, rout, incl, phi, ([xcen, 0]))
        pnew = PathPatch(path, facecolor=ringcol, ec='none', alpha=ttau, zorder=-10)
        p.add_patch(pnew)

    p.set_xlabel("Time [days]")
    p.ticklabel_format(axis='x', scilimits=(-2, 9))
    p.set_ylabel("Time [days]")


def draw_rings(r, tau, hjd_central, i, phi, p):
    """Draw pixel based ring image. Can be slow for large images."""
    xcen = hjd_central
    ycen = 0.0
    xrang = np.array(50.)
    yrang = np.array(50.)
    resol = 0.05

    nx = xrang / resol
    ny = yrang / resol

    xc = (nx - 1) / 2.
    yc = (ny - 1) / 2.

    agr = np.mgrid[:ny, :nx]

    agr[0] -= yc
    agr[1] -= xc
    agr = agr * resol
    # a is now in units of DAYS with the origin centered on J1407b

    tau_disk, grad = ellipse_nest(agr, i, phi, r, tau)

    ext = [np.min(agr[1]) + xcen, np.max(agr[1]) + xcen, np.min(agr[0]) + ycen, np.max(agr[0]) + ycen]
    p.imshow(tau_disk, extent=ext)
    p.scatter(hjd_central, 0)
    p.set_xlabel("Time [days]")
    p.ticklabel_format(axis='x', scilimits=(-2, 9))
    p.set_ylabel("Time [days]")


def draw_rings_lines(hjdm, r, p, dstar):
    """Plot line of stellar track across disk and point of rmin."""
    solar = mpl.patches.Rectangle((hjdm - 50, r[0] - dstar / 2.), 100, dstar, color='b', zorder=-15)
    p.add_patch(solar)
    p.scatter(hjdm, r[0], color='g', s=40, zorder=20)


def plot_tilt_faceon(radius, taun, r, tmin, dstar):
    """Draw tilted ring system and face on rings."""
    fig_rings = plt.figure('rings', figsize=(14, 7))
    p1 = fig_rings.add_subplot(121)
    draw_rings(radius, y_to_tau(taun), r[1], r[2], r[3], p1)
    p2 = fig_rings.add_subplot(122)
    draw_rings(radius, y_to_tau(taun), r[1], 0.0, 0.0, p2)
    draw_rings_lines(tmin, r, p1, dstar)


def plot_tilt_faceon_vect(radius, taun, r, tmin, rstart, rend, dstar):
    """Draw tilted ring system and face on rings."""
    fig_ringsv = plt.figure('ringsv', figsize=(14, 7))
    p1v = fig_ringsv.add_subplot(121, aspect='equal')
    draw_rings_vector(radius, y_to_tau(taun), r[1], r[2], r[3], p1v)
    p2v = fig_ringsv.add_subplot(122, aspect='equal')
    draw_rings_vector(radius, y_to_tau(taun), r[1], 0.0, 0.0, p2v)
    #    if (badrings):
    draw_badrings(rstart, rend, r[1], r[2], r[3], p1v)
    draw_badrings(rstart, rend, r[1], 0.0, 0.0, p2v)
    draw_rings_lines(tmin, r, p1v, dstar)


def plot_gradient_fit(t, f, fn, xt, yt, p):
    """
    Plot gradient fit.
    :param t:
    :param f: gradient of fit at points of measurement
    :param fn:
    :param xt:
    :param yt:
    :param p:
    :return: none
    """

    p.plot(xt, yt, lw=3.0, color='black', zorder=1)
    p.scatter(t, f, facecolor='1.0', s=60, color='black', zorder=2, lw=1)
    p.scatter(t, f, facecolor='None', s=60, color='black', zorder=3, lw=1)
    p.scatter(t, fn, facecolor='0.0', s=60, zorder=4, lw=1)
    p.vlines(t, f, fn, zorder=1, lw=2, color='0.5', linestyles='dotted')
    p.set_xlabel('HJD - 2450000 [Days]')
    p.set_ylabel('Normalized gradient (1/days)')
    p.set_title('Measured gradients and fit from J1407')


############################
## TESTING
############################

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib import interactive
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    interactive(True)
    import matplotlib as mpl

    mpl.rc('image', interpolation='nearest', origin='lower', cmap='gray')

    a = np.mgrid[:201, :401]
    yc = 101
    xc = 201
    a[0] -= yc
    a[1] -= xc

    i_test = 70.
    phi_test = 10.

    # testing ring_to_sky and the inverse

    fig1 = plt.figure("ring_to_sky()", figsize=(6, 6))
    pr = fig1.add_subplot(111)

    p = np.linspace(0, 2 * np.pi, 20)

    rad = 2.
    X = rad * np.cos(p)
    Y = rad * np.sin(p)

    pr.set_xlim(-3, 3)
    pr.set_ylim(-3, 3)

    plt.scatter(X, Y)
    (Xp, Yp, Rp) = ring_to_sky(X, Y, 90. - i_test, phi_test)
    (X2, Y2, R) = sky_to_ring(Xp, Yp, 90. - i_test, phi_test)

    plt.scatter(Xp, Yp, color='red')
    plt.plot(X2, Y2, color='blue')

    # testing ellipse()
    fig2 = plt.figure("ellipse()")
    axre = fig2.add_subplot(211)
    axre2 = fig2.add_subplot(212)
    (im) = ellipse(a, i_test, phi_test)
    axre.imshow(im, origin='lower', aspect='equal', interpolation='nearest')
    axre2.imshow((im < 50), origin='lower', aspect='equal', interpolation='nearest')

    # testing ellipse_nest() simply
    fig = plt.figure("ellipse_nest() 3 args")
    axre = fig.add_subplot(211)
    axre2 = fig.add_subplot(212)
    (im, tan) = ellipse_nest(a, i_test, phi_test)
    axre.imshow(im, origin='lower', aspect='equal', interpolation='nearest')
    axre2.imshow((im < 50), origin='lower', aspect='equal', interpolation='nearest')

    # testing ellipse_nest()

    fig = plt.figure("ellipse_nest() 5 args")
    axre = fig.add_subplot(211)
    axre2 = fig.add_subplot(212)
    (im, tangent) = ellipse_nest(a, i_test, phi_test, (50, 100, 150), (0.1, 0.2, 0.3))
    axre.imshow(im, origin='lower', aspect='equal', interpolation='nearest')
    axre2.imshow(tangent, origin='lower', aspect='equal', interpolation='nearest')

    length = 20

    # pick out points from the tangent grid and draw vector arrows
    yt, xt = np.mgrid[:10, :20] * 20

    eliang = tangent[yt, xt]

    xarr = length * np.cos(eliang)
    yarr = length * np.sin(eliang)

    axre.scatter(xt, yt, c='green', s=75)
    axre.quiver(xt, yt, xarr, yarr, facecolor='white', edgecolor='black', linewidth=0.5)

    # testing the parametric version

    t = np.linspace(0, 2 * np.pi, num=100)

    phit = phi_test * np.pi / 180.
    a = 100.
    b = a * np.cos(i_test * np.pi / 180.)
    (X, Y, grad) = ellipse_para(a, b, phit, t)

    xarrp = length * np.cos(grad)
    yarrp = length * np.sin(grad)

    axre.plot(X + xc, Y + yc, c='white', linewidth=3)

    (dX0, dY0, tangang) = ellipse_tangents(a, b, phit)

    (x1, y1, grad1) = ellipse_para(a, b, phit, dX0)
    (x2, y2, grad2) = ellipse_para(a, b, phit, dY0)

    print("grad check: %f" % (grad1 * 180. / np.pi))
    print("grad check: %f" % (grad2 * 180. / np.pi))

    # now move up an Y=const. line to the dY0 position
    a_test = np.mgrid[:1, :401]
    xx = np.arange(401)
    a_test[0] += int(y2)
    a_test[1] -= xc

    # test the angle along y=0
    tlen = 20
    tx = 30
    ty = 0
    tdx = tlen * np.cos(tangang)
    tdy = tlen * np.sin(tangang)

    # we should plot xc and see test_tan
    (test_im, test_tan) = ellipse_nest(a_test, i_test, phi_test)
    fig5 = plt.figure("ellipse_nest()")
    ax5 = fig5.add_subplot(211)
    ax6 = fig5.add_subplot(212)
    ax5.scatter(xx, test_im)
    ax6.scatter(xx, test_tan * 180. / np.pi)
    ax6.scatter(xx, np.sin(test_tan) * 180)

    ttx = np.array((tx, tx + tdx))
    tty = np.array((ty, ty + tdy))
    axre.plot(ttx + xc, tty + yc, c='blue')
    axre.scatter(tx + xc, ty + yc, c='blue')

    # vertical tangent
    axre.scatter(x1 + xc, y1 + yc, c='red', s=65)

    # horizontal tangent
    axre.scatter(x2 + xc, y2 + yc, c='orange', s=65)

    axre.quiver(X + xc, Y + yc, xarrp, yarrp, facecolor='blue', edgecolor='black', linewidth=0.5)

    # test drawing translucent rings using Patches
    fig2 = plt.figure("ring_patch()")
    axre = fig2.add_subplot(111)

    axre.set_xlim(-1, 1)
    axre.set_ylim(-1, 1)

    # plot random grid of background points to highlight alpha
    xx = np.random.random(100)
    yy = np.random.random(100)

    plt.scatter(xx, yy)

    for r in np.arange(0, 1, 0.1):
        path = ring_patch(r, r + 0.1, i_test, phi_test)
        newp = PathPatch(path, color='red', alpha=r, ec='none')
        axre.add_patch(newp)

    # test converting ring angles to ring tilt
    print(ring_tilt(i_test, phi_test))

    # test
    fig = plt.figure("make_star_limbd()")
    import matplotlib.gridspec as gridspec

    gs = gridspec.GridSpec(2, 2)
    (ax) = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[0, 1])
    ax4 = plt.subplot(gs[1, 1])
    fig.add_subplot(ax)
    fig.add_subplot(ax2)
    fig.add_subplot(ax3)
    fig.add_subplot(ax4)
    ax.imshow(make_star_limbd(15))
    ax2.imshow(make_star_limbd(25))
    ax3.imshow(make_star_limbd(15, 0.6))
    ax4.imshow(make_star_limbd(25, 0.6))
    print(make_star_limbd(25, 0.6).shape)

    plt.show()
    input('press return to continue')
