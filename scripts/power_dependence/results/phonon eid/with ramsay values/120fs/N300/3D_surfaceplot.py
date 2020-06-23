import numpy as np
import glob
from mayavi import mlab
from mayavi.modules.axes import Axes

data_files=glob.glob("*.txt")
data_files.sort()


wl=np.zeros(len(data_files))
pulse_area=np.loadtxt(data_files[0],usecols=(0,))
occupation=np.zeros((len(pulse_area),len(data_files)))

for i in range(len(data_files)):
    wl[i]=float(data_files[i].rsplit('_')[2].rstrip("nm.txt"))
    occupation[:,i]=0.5*(np.loadtxt(data_files[i],usecols=(3,))+1)



mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

# Visualize the points
pts = mlab.points3d(wl, pulse_area, occupation, scale_mode='none', scale_factor=0.2)

# Create and visualize the mesh
mesh = mlab.pipeline.delaunay2d(pts)
surf = mlab.pipeline.surface(mesh)

#mlab.view(47, 57, 8.2, (0.1, 0.15, 0.14))
mlab.show()