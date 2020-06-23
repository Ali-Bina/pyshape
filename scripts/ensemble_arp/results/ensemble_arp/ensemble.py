#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:44:49 2020

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 
import glob
#from file_loader import *
import os


dipole_moments=np.linspace(16,34,18,dtype=int)
omega_0=np.around(np.arange(1.04107,1.09107,0.002),decimals=6)

omega_fwidth=30
dipole_hwidth=2.0
centre_dipole=25
centre_omega=1.06607

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

def gauss(t,b,c):                                                                #gaussian with no background
    return (np.exp(-1*(1/(2*(b**2/2.355**2)))*(t-c)**(2)))

def lorrentzian(x,b,c):  #b: half-width
    return (b/np.pi)*((x-c)**2+b**2)**(-1)/max((b/np.pi)*((x-c)**2+b**2)**(-1))

weight=np.zeros((len(dipole_moments),len(omega_0)))
ensemble=np.zeros((len(dipole_moments),len(omega_0)),dtype=object)
final_sum=np.zeros((len(dipole_moments),len(omega_0)))

def array_element(x,y):
    loc=os.path.abspath(os.getcwd())
    os.chdir(loc+'/'+x+"/")
    data_files=glob.glob("*"+y+".txt")
    data_files.sort()
    chirp=np.zeros(len(data_files))
    pulse_area=np.loadtxt(data_files[0],usecols=(0,),skiprows = 0,)
    occupation=np.zeros((len(pulse_area),len(data_files)))
    for i in range(len(data_files)):
#        print data_files[i]
        chirp[i]=float(data_files[i].rsplit('_')[1].rstrip("fs2"))
        occupation[:,i]=0.5*(np.loadtxt(data_files[i],usecols=(3,))+1)
    os.chdir(ROOT_DIR)
    return occupation


def get_values(x,y):
    loc=os.path.abspath(os.getcwd())
    x=str(x)+"D"
    y=str(y)+'0'
    os.chdir(loc+'/'+x+"/")
    data_files=glob.glob("*"+y+".txt")
    data_files.sort()
    chirp=np.zeros(len(data_files))
    pulse_area=np.loadtxt(data_files[0],usecols=(0,),skiprows = 0,)
#    occupation=np.zeros((len(pulse_area),len(data_files)))
    for i in range(len(data_files)):
#        print data_files[i]
        chirp[i]=float(data_files[i].rsplit('_')[1].rstrip("fs2"))
#        occupation[:,i]=0.5*(np.loadtxt(data_files[i],usecols=(3,))+1)
    os.chdir(ROOT_DIR)
    return pulse_area,chirp
    

pulse_area=np.linspace(0,4.0,40)#get_values(dipole_moments[5],omega_0[5])[0]
chirp=get_values(dipole_moments[9],omega_0[12])[1]

for k in range(len(dipole_moments)):
    for l in range(len(omega_0)):
        weight[k,l]=lorrentzian(dipole_moments,dipole_hwidth,centre_dipole)[k]*gauss(omega_0,omega_fwidth*1E-3,centre_omega)[l]

for i in range(len(dipole_moments)):
    for j in range(len(omega_0)):
        os.chdir(ROOT_DIR)
        a=str(dipole_moments[i])+"D"
        b=str(omega_0[j])+'0'
        ensemble[i,j]=array_element(a,b)
#        print x,y
weighted_ensemble=np.multiply(ensemble,weight)

final_result=np.sum(np.sum(weighted_ensemble,axis=0),axis=0)/np.amax(np.sum(np.sum(weighted_ensemble,axis=0),axis=0))
        

        

##plt.plot(omega_0,gauss(omega_0,2*1E-3,omega_0[5]))
#plt.plot(dipole_moments,lorrentzian(dipole_moments,dipole_moments[5]/3.0,dipole_moments[5]))
#
plt.figure()
A,B=np.meshgrid(pulse_area,chirp)
plt.contour(A, B,np.transpose(final_result),levels=[0.95],colors=('k',),linewidths=(2,),linestyles=('--',))
plt.contourf(A, B,np.transpose(final_result),100,cmap='plasma')
plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    right=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False,
    labelleft=False)
#plt.clim(np.amin(BX+Y),np.amax(BX+Y))
plt.xticks([0,1,2,3,4],fontsize=20)
plt.yticks([0,10000,20000,30000],fontsize=20)
cbar = plt.colorbar()
cbar.set_ticks([0,0.25,0.5,0.75,1])
cbar.ax.tick_params(labelsize=20)
#plt.xlabel("Pulse Area ($\pi$ Units)",fontsize=20)
#plt.ylabel('Chirp ($fs^2$)',fontsize=20)
plt.tight_layout(pad=0.1)
plt.show()
#plt.savefig("ensemble_30meV.png") 