#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:52:57 2016

@author: ajan
"""

#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 
import glob
import os

os.chdir("/home/ajan/Documents/pyshape_edited/scripts/power_dependence/results/20180601/more_resolved/N6140/")
data_files=glob.glob("*.txt")
data_files.sort()


wl=np.zeros(len(data_files))
pulse_area=np.loadtxt(data_files[0],usecols=(0,))
occupation_N=np.zeros((len(pulse_area),len(data_files)))

for i in range(len(data_files)):
    wl[i]=float(data_files[i].rsplit('_')[2].rstrip("nm.txt"))
    occupation_N[:,i]=0.5*(np.loadtxt(data_files[i],usecols=(3,))+1)



plt.figure()
A,B=np.meshgrid(wl,pulse_area)
plt.contourf(A, B,occupation_N,300)
plt.clim(0,1)
plt.colorbar(pad=0.05)
plt.xlabel("Wavelength(nm)")
plt.ylabel(r'$\sqrt{Power(\mu W)}$')
#plt.rc('xtick',labelsize=12)
#plt.rc('ytick',labelsize=12)
plt.ticklabel_format(useOffset=False)
plt.show()
plt.savefig("TL_N.png")  

