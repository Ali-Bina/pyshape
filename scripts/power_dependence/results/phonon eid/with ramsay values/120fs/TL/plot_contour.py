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
    occupation[:,i]=0.5*(np.loadtxt(data_files[i],usecols=(3,))+1)



#plt.figure()
#A,B=np.meshgrid(wl,pulse_area)
#plt.contourf(A, B,occupation,300)
#plt.clim(0,1)
#plt.colorbar(pad=0.05)
#plt.xlabel("Wavelength(nm)")
#plt.ylabel(r'$\sqrt{Power(\mu W)}$')
#plt.rc('xtick',labelsize=12)
#plt.rc('ytick',labelsize=12)
#plt.ticklabel_format(useOffset=False)
#plt.show()
##plt.savefig("TL.png")  
count=0
plt.figure(figsize=(3,25))
for j in range(0,30,3):
    count=count+1
    plt.subplot(10,1,count)
    plt.rc('xtick',labelsize=8)
    plt.rc('ytick',labelsize=8)
    plt.plot(wl,occupation[j,:])
    plt.ylim(0,1.1)
    plt.locator_params(nbins=4)
    
    plt.text(1170,0.9,"Pulse Area ="+str(round(pulse_area[j],2)),fontsize=8)
#    plt.legend(loc='best',frameon=False, prop={'size': 8})
    plt.tight_layout()
plt.show()
plt.savefig("pulse_area_slice.svg")