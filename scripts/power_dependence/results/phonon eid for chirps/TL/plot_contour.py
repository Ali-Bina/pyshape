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

data_files=glob.glob("*.txt")
data_files.sort()


wl=np.zeros(len(data_files))
pulse_area=np.loadtxt(data_files[0],usecols=(0,))
occupation=np.zeros((len(pulse_area),len(data_files)))

for i in range(len(data_files)):
    wl[i]=float(data_files[i].rsplit('_')[2].rstrip("nm.txt"))
    occupation[:,i]=np.loadtxt(data_files[i],usecols=(2,))



plt.figure()
A,B=np.meshgrid(wl,pulse_area)
plt.contourf(A, B,occupation,300)
#plt.clim(np.amin(BX+Y),np.amax(BX+Y))
plt.colorbar(pad=0.12)
plt.xlabel("Wavelength(nm)")
plt.ylabel(r'$\sqrt{Power(\mu W)}$')
plt.show()
plt.savefig("TL.png")  

