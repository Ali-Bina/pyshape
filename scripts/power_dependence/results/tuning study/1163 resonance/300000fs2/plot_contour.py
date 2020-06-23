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
    occupation[:,i]=(np.loadtxt(data_files[i],usecols=(3,)))



energy_ev=1239.84/wl
detuning = (energy_ev-1.06607)*1000

plt.style.use('physrev')
plt.figure()
A,B=np.meshgrid(detuning,pulse_area)
plt.contourf(A, B,occupation,300)
plt.clim(0,1)
cb=plt.colorbar(pad=0.05,ticks=np.arange(0,1,0.2))
cb.ax.tick_params(labelsize=fs)
plt.xlabel("Detuning (meV)",fontsize=fs)
plt.ylabel(u'\u221A'+'Power('+r'$\mu$'+'W)',fontsize=fs)
#plt.ylabel(r'$\mathrm{\sqrt{Power(\mu W)}}$')
plt.xticks(np.linspace(-15,15,7),fontsize=fs)
plt.yticks(np.arange(0,26,5),fontsize=fs)
plt.ylim(0,2.5)
plt.show()

