#! /home/ajan/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:05:16 2017

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 
import os
import glob

#Tl_path="/home/ajan/Documents/pyshape_edited/scripts/power_dependence/results/tuning study/TL/"

#TL_files=glob.glob(Tl_path+"*.txt")
chirp_files=glob.glob("*.txt")

#TL_files.sort()
chirp_files.sort()


#Tl_data=np.zeros([30,len(TL_files)])
chirp_data=np.zeros([30,len(chirp_files)])
pulse_area=np.loadtxt(chirp_files[0],usecols=(0,),skiprows = 0,)
wl=[]
omega_c=[]
chirp=[]

for i in range(len(chirp_files)):
#    Tl_data[:,i]=np.loadtxt(TL_files[i],usecols=(3,),skiprows = 0,)
    chirp_data[:,i]=0.5*(np.loadtxt(chirp_files[i],usecols=(3,),skiprows = 0,)+1)
    wl.append(str(chirp_files[i].rsplit('_')[2].rstrip("nm")))
    omega_c.append(str(chirp_files[i].rsplit('_')[3].rstrip(".txt")))
    chirp.append(str(chirp_files[i].rsplit('_')[1].rstrip("fs2")))

plt.figure(figsize=(10,20))
for j in range(7):
    plt.subplot(3,3,j+1)
    plt.plot(pulse_area,chirp_data[:,j],label='P300',linewidth=2,color='b')
    plt.ylim(0,1.1)
    plt.locator_params(nbins=4)
#    plt.text(2.0,0.2,str(omega_c[j]*1000)+" meV",fontsize=8)
    plt.legend(loc='best',frameon=False, prop={'size': 8})
    plt.tight_layout()
plt.show()
#plt.savefig("chirp_figure1.svg")