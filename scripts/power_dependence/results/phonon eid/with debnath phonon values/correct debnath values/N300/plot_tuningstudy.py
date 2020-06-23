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

for i in range(len(chirp_files)):
#    Tl_data[:,i]=np.loadtxt(TL_files[i],usecols=(3,),skiprows = 0,)
    chirp_data[:,i]=0.5*(np.loadtxt(chirp_files[i],usecols=(3,),skiprows = 0,)+1)
    wl.append(str(chirp_files[i].rsplit('_')[2].rstrip("nm.txt")))

plt.figure(figsize=(3,25))
for j in range(10):
    plt.subplot(10,1,j+1)
    plt.plot(pulse_area,chirp_data[:,j])
    plt.ylim(0,1.1)
    plt.locator_params(nbins=4)
    plt.text(3.2,0.5,wl[j],fontsize=15)
    plt.tight_layout()
plt.show()
plt.savefig("chirp_figure.svg")