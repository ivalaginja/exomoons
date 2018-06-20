"""
This is an example program for how to build a grid with GridSpec. First step in creating grid_lightcurve_plots.py, but not used anymore.
"""

# from http://matplotlib.org/users/gridspec.html

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def name_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")

def make_x_ticklabels_invisible(fig):
    for i, ax in enumerate((ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax15,ax17)):
        for tl in ax.get_xticklabels():
            tl.set_visible(False)

def make_y_ticklabels_invisible(fig):
    for i, ax in enumerate((ax1,ax3,ax4,ax5,ax6,ax7,ax9,ax10,ax11,ax12,ax13,ax15,ax16,ax17,ax18)):
        for tl in ax.get_yticklabels():
            tl.set_visible(False)


f = plt.figure()

plt.suptitle("Building a grid")

gs1 = GridSpec(2, 3)
gs1.update(top=0.98, bottom=0.68, hspace=0.05)
ax1 = plt.subplot(gs1[0,0])
ax2 = plt.subplot(gs1[1,0])
ax3 = plt.subplot(gs1[0,1])
ax4 = plt.subplot(gs1[1,1])
ax5 = plt.subplot(gs1[0,2])
ax6 = plt.subplot(gs1[1,2])

gs2 = GridSpec(2, 3)
gs2.update(top=0.64, bottom=0.36, hspace=0.05)
ax7 = plt.subplot(gs2[0,0])
ax8 = plt.subplot(gs2[1,0])
ax9 = plt.subplot(gs2[0,1])
ax10 = plt.subplot(gs2[1,1])
ax11 = plt.subplot(gs2[0,2])
ax12 = plt.subplot(gs2[1,2])

gs3 = GridSpec(2, 3)
gs3.update(top=0.32, bottom=0.03, hspace=0.05)
ax13 = plt.subplot(gs3[0,0])
ax14 = plt.subplot(gs3[1,0])
ax15 = plt.subplot(gs3[0,1])
ax16 = plt.subplot(gs3[1,1])
ax17 = plt.subplot(gs3[0,2])
ax18 = plt.subplot(gs3[1,2])

name_axes(plt.gcf())
make_x_ticklabels_invisible(plt.gcf())
make_y_ticklabels_invisible(plt.gcf())


plt.show()

