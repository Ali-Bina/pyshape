#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title="M25_S4 "
filename= 'phase.txt'
filename2='freq.txt'
filename3 = 'Efield_intensity.txt'

var1_name="name"
var2_name="name"
var3_name="name"

var=np.genfromtxt(filename) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename2,skiprows = 0) #2nd variable -yaxis1 check usecol,delimiter
var3=np.genfromtxt(filename3)
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter
#Int=(np.absolute(var))**2

#var2=1240.0/var2
"""
#Delay position to ps
speed_of_light = 299792458.0
zero_delay = -13800.0  #zero_delay
var[:,0] = -1E12*2*(var[:,0] - zero_delay)*1E-6/(speed_of_light)
"""
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(var2,var/max(var),'r',linestyle='--') #add 3rd variable 
ax2.plot(var2,var3/max(var3)) 
#plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name)
ax1.set_ylabel(var2_name)

plt.xlim(0.95,1.05)
ax1.set_ylim(-5,5)
ax2.set_ylim(0,1)

#plt.ylim(0,max(Int))
plt.title(Title)
plt.grid(True)
plt.show()    
