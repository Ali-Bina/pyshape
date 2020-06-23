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
#var2=1240.0/var2        

wl=np.zeros(len(data_files))
pulse_area=np.loadtxt(data_files[0],usecols=(0,))
occupation=np.zeros((len(pulse_area),len(data_files)))

for i in range(len(data_files)):
    wl[i]=float(data_files[i].rsplit('_')[2].rstrip("nm.txt"))
    occupation[:,i]=0.5*(np.loadtxt(data_files[i],usecols=(3,))+1)

#
#plt.style.use('physrev')
##plt.figure()
#fig, ax1 = plt.subplots()
#A,B=np.meshgrid(wl,pulse_area)
#plt.contourf(A, B,occupation,300)
##plt.clim(np.amin(BX+Y),np.amax(BX+Y))
##plt.clim(0,4250)
#cb=plt.colorbar(pad=0.05,ticks=np.arange(0,1.1,0.2))
##plt.clim(0,4250)
#cb.ax.tick_params(labelsize=fs)
#plt.xlabel("Wavelength(nm)",fontsize=fs)
#plt.ylabel('Pulse Area ($\pi $)',fontsize=fs)
##plt.ylabel(r'$\mathrm{\sqrt{Power(\mu W)}}$')
#plt.xlim(1148,1183)
#plt.xticks(np.arange(1150,1181,10),fontsize=fs)
#plt.yticks(fontsize=fs)
#plt.ylim(0.1,2.5)
#plt.show()




energy_ev=1239.84/wl
detuning = (energy_ev-1.06607)*1000
det=(var2-1.06607)*1000


plt.style.use('physrev')
fig, ax1 = plt.subplots()
A,B=np.meshgrid(detuning,pulse_area)
plt.contourf(A, B,occupation,300)
plt.clim(0,1)
cb=plt.colorbar(pad=0.05,ticks=np.arange(0,1.1,0.25))
cb.ax.tick_params(labelsize=fs)
plt.xlabel("Detuning (meV)",fontsize=fs)
plt.ylabel(r'Pulse Area ($\pi$)',fontsize=fs)
#plt.ylabel(r'$\mathrm{\sqrt{Power(\mu W)}}$')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    right=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True,
    )
plt.xticks(np.linspace(-15,15,7),fontsize=fs)
plt.yticks([0,0.5,1,1.5,2,2.5],fontsize=fs)
plt.ylim(0.1,2.5)
plt.xlim(-15,15)
plt.show()
#
##plt.figure()
##plt.plot(det,var)
ax2=ax1.twinx()                          
ax2.plot(det,var,'w',linewidth=2.5)
##ax2.set_xlim(min(Wavelength),max(Wavelength))
##ax2.set_ylabel('T_0 (mV)')
plt.yticks([])
ax2.set_xlim(-15,15)
#plt.savefig("Rabi_theory.png")  
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
#    
#plt.show()
#plt.savefig("pulse_area_slice.svg")
#plt.plot(wl,occupation[9,:],label=r"$\phi ''=0$ ",linestyle='-.')
##plt.plot(var2,var/max(var),label="spectrum")
#plt.xlim(1130,1190)
#plt.ylim(0,1.1)
#plt.legend(loc='best',prop={'size': 20})
#plt.xlabel("Wavelength(nm)",fontsize=20)
#plt.ylabel(r"$\rho_{11}$",fontsize=20)
#plt.show()

