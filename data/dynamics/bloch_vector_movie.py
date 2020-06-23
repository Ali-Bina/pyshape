#! /home/ajan/anaconda/bin/python


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import os
from itertools import product, combinations
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from subprocess import call


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
detuning=np.loadtxt("/home/ajan/Documents/pyshape_edited/data/pulse/detuning.txt",float)
rabi_freq=np.loadtxt("/home/ajan/Documents/pyshape_edited/data/pulse/rabi_freq.txt",float)


#x1 = 2.0*rho_1[:, 6]
#y1 = -2.0*rho_1[:, 7]
#z1 = (rho_1[:, 2] - rho_1[:, 0])

x1 = S[:,0]
y1 = S[:,1]
z1 = S[:,2]

#x2 = 2.0*rho_2[:, 6]
#y2 = -2.0*rho_2[:, 7]
#z2 = rho_2[:, 2] - rho_2[:, 0]

dot1_desired = np.loadtxt("rhod_1.txt", float)
#dot2_desired = np.loadtxt("rhod_2.txt", float)

Lambda_x=-1*rabi_freq[0::10]
Lambda_y=np.zeros_like(Lambda_x)
Lambda_z=detuning[0::10]
control_x=np.zeros_like(Lambda_x)
control_y=np.zeros_like(Lambda_x)
control_z=np.zeros_like(Lambda_x)

for i in range(len(Lambda_x)):
    control_x[i]=Lambda_x[i]/(Lambda_x[i]**2+Lambda_z[i]**2)**0.5
    control_z[i]=-1*Lambda_z[i]/(Lambda_x[i]**2+Lambda_z[i]**2)**0.5
    

s_x1 = 2.0*dot1_desired[6]
s_y1 = -2.0*dot1_desired[7]
s_z1 = dot1_desired[2] - dot1_desired[0]

#s_x2 = 2.0*dot2_desired[6]
#s_y2 = -2.0*dot2_desired[7]
#s_z2 = dot2_desired[2] - dot2_desired[0]

#---------------------------------------------------------------------
count=0

print "creating images...."
for i in range(1800,3000):
#    print i,count
    fig = plt.figure()
    matplotlib.rcParams.update({'font.size': 17})
    
    ax = fig.add_subplot(111, projection='3d')
    (xs,ys,zs) = drawSphere(0,0,0,1)
    ax.plot_wireframe(xs,ys,zs, color = "0.8")    
    ax.quiver(x1[i], y1[i], z1[i],0,0,0)
    ax.quiver(Lambda_x[i]/(Lambda_x**2+Lambda_z**2)**0.5,Lambda_y[i]/(Lambda_x**2+Lambda_z**2)**0.5, Lambda_z[i]/(Lambda_x**2+Lambda_z**2)**0.5,0,0,0)
    a = Arrow3D([0, x1[i]], [0, y1[i]], [0,z1[i]], mutation_scale=20,lw=1, arrowstyle="-|>", color="k")
    ax.add_artist(a)
    b= Arrow3D([0, control_x[i]], [0,control_y[i]], [0,control_z[i]], mutation_scale=20,lw=1, arrowstyle="-|>", color="r")
    ax.add_artist(b)
    #ax.plot(x2, y2, z2, '-r')
    #ax.plot(Lambda_x,Lambda_y,Lambda_z,'-r')
    
#    ax.scatter(x1[-1], y1[-1], z1[-1], color = 'b')#Final point, simulation
#    ax.scatter(x2[-1], y2[-1], z2[-1], color = 'r')#final point, simulation
#    ax.scatter(s_x1, s_y1, s_z1, color = 'g')#final point, desired
#    ax.scatter(s_x2, s_y2, s_z2, color = 'b')#final point, desired
    
    ax.set_xlabel(r'$S_x$')
    ax.set_ylabel(r'$S_y$')
    ax.set_zlabel(r'$S_z$')
    
    ax.set_xlim3d(-1.01, 1.01)
    ax.set_ylim3d(-1.01, 1.01)
    ax.set_zlim3d(-1.01, 1.01)
    ax.view_init(elev=25, azim=45)
    subplots_adjust(hspace=0.25, wspace=0.40, top=0.94, right=0.98, bottom=0.1, left=0.30)
    
    
    plt.title(str(z1[i]), fontsize = 18)
    plt.savefig("/home/ajan/Documents/pyshape_edited/data/dynamics/Untitled Folder/"+'{0:04d}'.format(count)+".png")
    count=count+1
    plt.close()

#call('echo $PWD',shell=True)
call('cd /home/ajan/Documents/pyshape_edited/data/dynamics/Untitled\ Folder/ && ffmpeg -y -start_number 0 -i %4d.png video.mp4',shell=True )
