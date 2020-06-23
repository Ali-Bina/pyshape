#! /home/ajan/anaconda/bin/python
from pylab import *
from mayavi.mlab import *


#pulse_area = loadtxt('script_results_axes_pulse_area.txt', skiprows=1)
#pulse_area = transpose(pulse_area)
##chirp = loadtxt('script_results_axes_chirp.txt', skiprows=1)
#chirp = transpose(chirp)
#fidelity = loadtxt('script_results.txt', skiprows=1)
#print pulse_area
#X, Y = meshgrid(chirp/1000000.0, pulse_area)
#contour3d(fidelity, contours=4, transparent=True)
#show()

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = Axes3D(fig)
pulse_area = loadtxt('script_results_axes_pulse_area.txt', skiprows=1)
pulse_area = transpose(pulse_area)
chirp = loadtxt('script_results_axes_chirp.txt', skiprows=1)
chirp = transpose(chirp)
X, Y = meshgrid(chirp/1000000.0, pulse_area)

fidelity = loadtxt('script_results.txt', skiprows=1)
Z = fidelity

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
ax.set_zlim(0.00, 1.00)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

ax.view_init(40,-150)
plt.show()

