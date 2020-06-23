from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import os
from itertools import product, combinations
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)
        
def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)



folder = './'

os.chdir(folder)

S=np.loadtxt("S.txt")
rho_1 = np.loadtxt("rho_1.txt", complex)
rho_2 = np.loadtxt("rho_2.txt", complex)

#x1 = 2.0*rho_1[:, 6]
#y1 = -2.0*rho_1[:, 7]
#z1 = rho_1[:, 2] - rho_1[:, 0]

x1 = S[0,0]
y1 = S[0,1]
z1 = S[0,2]

x2 = S[-1,0]
y2 = S[-1,1]
z2 = S[-1,2]


dot1_desired = np.loadtxt("rhod_1.txt", float)
dot2_desired = np.loadtxt("rhod_2.txt", float)

s_x1 = 2.0*dot1_desired[6]
s_y1 = -2.0*dot1_desired[7]
s_z1 = dot1_desired[2] - dot1_desired[0]

s_x2 = 2.0*dot2_desired[6]
s_y2 = -2.0*dot2_desired[7]
s_z2 = dot2_desired[2] - dot2_desired[0]

#---------------------------------------------------------------------
matplotlib.rcParams.update({'font.size': 17})

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
(xs,ys,zs) = drawSphere(0,0,0,1)
ax.plot_wireframe(xs,ys,zs, color = "0.8")
a = Arrow3D([0, x1], [0, y1], [0,z1], mutation_scale=20,lw=1, arrowstyle="-|>", color="k")
ax.add_artist(a)
b= Arrow3D([0, x2], [0,y2], [0,z2], mutation_scale=20,lw=1, arrowstyle="-|>", color="r")
ax.add_artist(b)

#ax.plot(x1, y1, z1, '-k')
#ax.plot(x2, y2, z2, '-r')

#ax.scatter(x1[-1], y1[-1], z1[-1], color = 'k')#Final point, simulation
#ax.scatter(x2[-1], y2[-1], z2[-1], color = 'r')#final point, simulation
#ax.scatter(s_x1, s_y1, s_z1, color = 'g')#final point, desired
#ax.scatter(s_x2, s_y2, s_z2, color = 'b')#final point, desired

ax.set_xlabel(r'$S_x$')
ax.set_ylabel(r'$S_y$')
ax.set_zlabel(r'$S_z$')

ax.set_xlim3d(-1.01, 1.01)
ax.set_ylim3d(-1.01, 1.01)
ax.set_zlim3d(-1.01, 1.01)
ax.view_init(elev=25, azim=45)
subplots_adjust(hspace=0.25, wspace=0.40, top=0.94, right=0.98, bottom=0.1, left=0.30)


plt.title(r"$ \Delta \phi = 1.0000 \pi$, $F = 0.249999$", fontsize = 18)


plt.show()

