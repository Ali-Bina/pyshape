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
#chirp_files.sort()


#Tl_data=np.zeros([30,len(TL_files)])
chirp_data_x=np.zeros([30,len(chirp_files)])
chirp_data_y=np.zeros([30,len(chirp_files)])
pulse_area=np.loadtxt(chirp_files[0],usecols=(0,),skiprows = 0,)
wl=[]
chirp=[]
prefix=chirp_files[0].rsplit('_')[0]+"_"
suffix="fs2_1170.004nm.txt"

for i in range(len(chirp_files)):
#    Tl_data[:,i]=np.loadtxt(TL_files[i],usecols=(3,),skiprows = 0,)
#    chirp_data_x[:,i]=np.loadtxt(chirp_files[i],usecols=(2,),skiprows = 0,)
#    chirp_data_y[:,i]=np.loadtxt(chirp_files[i],usecols=(3,),skiprows = 0,)
    wl.append(float(chirp_files[i].rsplit('_')[2].rstrip("nm.txt")))
    chirp.append(float(chirp_files[i].rsplit('_')[1].rstrip("fs2")))

chirp.sort()

for j in range(len(chirp_files)):
#    Tl_data[:,i]=np.loadtxt(TL_files[i],usecols=(3,),skiprows = 0,)
#    chirp_data_x[:,j]=np.loadtxt(prefix+str(chirp[j])+suffix,usecols=(2,),skiprows = 0,)
    chirp_data_y[:,j]=np.loadtxt(prefix+str(chirp[j])+suffix,usecols=(3,),skiprows = 0,)


plt.figure(figsize=(3,20))
for j in range(6):
    plt.subplot(6,1,j+1)
#    plt.plot(pulse_area,chirp_data_x[:,j],label="x")
    plt.plot(pulse_area,chirp_data_y[:,j],label="y")
    plt.ylim(0,1.1)
    plt.locator_params(nbins=4)
    plt.text(5,0.5,str(chirp[j])+"fs2")
    plt.tight_layout()
#plt.legend(bbox_to_anchor=(0.5, -0.1), loc=9, ncol=4)
plt.show()
plt.savefig("chirp_figure.svg")