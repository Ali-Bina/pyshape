#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:14:04 2019

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 
import os
import glob
from ThetaT_finder_n import chirp_data_n as chirp_n

#Tl_path="/home/ajan/Documents/pyshape_edited/scripts/power_dependence/results/tuning study/TL/"

#TL_files=glob.glob(Tl_path+"*.txt")
chirp_files=glob.glob("*.txt")

#TL_files.sort()
chirp_files.sort()


#Tl_data=np.zeros([30,len(TL_files)])
chirp_data_p=np.zeros([40,len(chirp_files)])
pulse_area=np.loadtxt(chirp_files[0],usecols=(0,),skiprows = 0,)
wl=[]
omega_c=[]
chirp=[]
theta_T=np.zeros([len(chirp_files)])

for i in range(len(chirp_files)):
#    Tl_data[:,i]=np.loadtxt(TL_files[i],usecols=(3,),skiprows = 0,)
    chirp_data_p[:,i]=0.5*(np.loadtxt(chirp_files[i],usecols=(3,),skiprows = 0,)+1)
    wl.append(str(chirp_files[i].rsplit('_')[2].rstrip("nm")))
    omega_c.append(float(chirp_files[i].rsplit('_')[3].rstrip(".txt")))
    chirp.append(str(chirp_files[i].rsplit('_')[1].rstrip("fs2")))

omega_c=np.array(omega_c)
difference=chirp_data_p-chirp_n

for j in range(len(chirp_files)):
    v=np.where(difference[:,j]>0.01)
    try:
        theta_T[j]=pulse_area[int(v[-1][-1]+1)]
    except IndexError:
        theta_T[j]=0.0
    
plt.plot(6.727/(1000*omega_c*2**0.5),theta_T,marker='o',label='Sphere_P')
plt.xlabel(r'$\hbar \omega_c$')
plt.ylabel(r'$\theta_T$')