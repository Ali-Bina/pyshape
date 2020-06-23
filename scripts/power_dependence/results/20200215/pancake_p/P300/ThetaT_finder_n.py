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

#Tl_path="/home/ajan/Documents/pyshape_edited/scripts/power_dependence/results/tuning study/TL/"

#TL_files=glob.glob(Tl_path+"*.txt")
path="../N300/"
chirp_files=glob.glob(path+"*.txt")

#TL_files.sort()
chirp_files.sort()


#Tl_data=np.zeros([30,len(TL_files)])
chirp_data_n=np.zeros([40,len(chirp_files)])
pulse_area=np.loadtxt(chirp_files[0],usecols=(0,),skiprows = 0,)
wl=[]
omega_c=[]
chirp=[]

for i in range(len(chirp_files)):
#    Tl_data[:,i]=np.loadtxt(TL_files[i],usecols=(3,),skiprows = 0,)
    chirp_data_n[:,i]=0.5*(np.loadtxt(chirp_files[i],usecols=(3,),skiprows = 0,)+1)
    wl.append(str(chirp_files[i].rsplit('_')[2].rstrip("nm")))
    omega_c.append(str(chirp_files[i].rsplit('_')[3].rstrip(".txt")))
    chirp.append(str(chirp_files[i].rsplit('_')[1].rstrip("fs2")))
    
    
#for j in range(20):
#    plt.subplot(5,4,j+1)
#    plt.plot(pulse_area,chirp_data_n[:,j],label='N300',linewidth=2,linestyle='--',color='r')
#    plt.ylim(0,1.1)
#    plt.locator_params(nbins=4)
##    plt.text(6.0,0.8,str(round(float(omega_c[j])*1000,5))+' meV',fontsize=8)
#    plt.legend(loc='lower right',frameon=False, prop={'size': 8})
#    plt.tight_layout()
#plt.show()