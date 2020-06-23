#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.signal import chirp, find_peaks, peak_widths

Title="M25_S4 "
filename= 'Efield_intensity.txt'
filename2='freq.txt'

var1_name="name"
var2_name="name"
var3_name="name"

var=np.genfromtxt(filename) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename2,skiprows = 0) #2nd variable -yaxis1 check usecol,delimiter
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter
#Int=(np.absolute(var))**2

detuning = (var2-1.06607)*1000

var=var/max(var)
x=np.where(var==max(var))

peaks, _ = find_peaks(var)
results_half = peak_widths(var, x[0], rel_height=0.5)
"""
#Delay position to ps
speed_of_light = 299792458.0
zero_delay = -13800.0  #zero_delay
var[:,0] = -1E12*2*(var[:,0] - zero_delay)*1E-6/(speed_of_light)
"""
#plt.figure()
plt.plot(detuning,var/max(var),'r',linestyle='--',linewidth=3.0,label='laser') #add 3rd variable 
#plt.axis([1200, 1320, 0, 1.1])
#plt.xlabel(var1_name)
#plt.ylabel(var2_name)
plt.xlim(-25,25)
#plt.ylim(0,max(Int))
#plt.title(Title)
#plt.grid(True)
plt.legend()
plt.show()    
