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


data_files=glob.glob("occupation*")
data_files.sort()
fs=18

filename= 'Efield_intensity.txt'
filename2='freq.txt'
var=np.genfromtxt(filename) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename2,skiprows = 0)
var2=1240.0/var2     

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
#plt.colorbar(pad=0.12)
#plt.xlabel("Wavelength(nm)")
#plt.ylabel(r'$\sqrt{Power(\mu W)}$')
#plt.rc('xtick',labelsize=12)
#plt.rc('ytick',labelsize=12)
#plt.ticklabel_format(useOffset=False)
#plt.show()
#plt.savefig("TL.png")  
#plt.style.use('physrev')
#count=0
#plt.figure(figsize=(3,25))
#for j in range(0,30,3):
#    count=count+1
#    plt.subplot(10,1,count)
#    plt.xticks([1130.0,1160.0,1190.0])
#    plt.rc('xtick',labelsize=12)
#    plt.rc('ytick',labelsize=12)
#    plt.plot(wl,occupation[j,:])
#    plt.ylim(0,1.1)
#    plt.locator_params(nbins=4)
#    
#    plt.text(1135,0.9,r"$\Theta $="+str(round(pulse_area[j],2)),fontsize=12)
##    plt.legend(loc='best',frameon=False, prop={'size': 8})
#    plt.tight_layout()
#plt.show()
#plt.savefig("pulse_area_slice.svg")

plt.plot(wl,occupation[27,:],label=r"$\phi ''<0$",linestyle='--')
#plt.plot(var2,var/max(var),label="spectrum")
plt.xlim(1130,1190)
plt.ylim(0,1.1)
plt.legend(loc='best',prop={'size': 20})
plt.xlabel("Wavelength(nm)",fontsize=20)
plt.ylabel(r"$\rho_{11}$",fontsize=20)
plt.show()