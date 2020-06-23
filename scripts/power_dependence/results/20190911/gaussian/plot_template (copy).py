#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title="M25_S4 "
filename= 'occupation_300000.0fs2_1148.002nm.txt' 

var1_name="name"
var2_name="name"
var3_name="name"
fs=36
pulse_area=np.loadtxt(filename,usecols=(0,)) #1st variable -x axis; usecols=(0,1)
var2=0.5*(np.loadtxt(filename,usecols=(3,))+1) #2nd variable -yaxis1 check usecol,delimiter
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter

#x = np.arange(0, 10, 0.1)
#y1 = 0.05 * x**2
#y2 = -1 *y1
#
#fig, ax1 = plt.subplots()
#
#ax2 = ax1.twinx()
#ax1.plot(x, y1, 'g-')
#ax2.plot(x, y2, 'b-')
#
#ax1.set_xlabel('X data')
#ax1.set_ylabel('Y1 data', color='g')
#ax2.set_ylabel('Y2 data', color='b')

#plt.show()
#
plt.style.use('physrev')
plt.plot(pulse_area,var2,linewidth=3.0) #add 3rd variable 
plt.axis([0, 8, 0, 1.1])
plt.xlabel('Pulse Area ($\pi $ units)',fontsize=fs)
plt.xticks([2,4,6],fontsize=fs)
plt.yticks([0,0.5,1.0],fontsize=fs)
#plt.ylabel(var2_name)
#plt.xlim(0,1024)
#plt.title(Title)
#plt.grid(True)
plt.show()    
