#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title=" "
filename= 'Theta_t.txt' 

var1_name="d (nm)"
var2_name=r"$\Theta_T$"
var3_name="name"

#var=np.loadtxt(filename,skiprows = 0,) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename,usecols=(1,),skiprows = 0) #2nd variable -yaxis1 check usecol,delimiter
var3=np.loadtxt(filename,usecols=(2,),skiprows = 0) #3rd variable -yaxis2 check usecol,delimiter

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
plt.plot(var2,var3,label="Pancake_p",marker='D',linewidth=2.0,markeredgecolor='g',color='g',markeredgewidth=2,markersize=10,markerfacecolor='None') #add 3rd variable 
#plt.axis([1200, 1320, 0, 1.1])
plt.xticks([4,6,8,10,12])
plt.yticks([1.0,2.0,3.0,4.0])
#plt.xlabel(var1_name)
#plt.ylabel(var2_name)
plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    right=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False,
    labelleft=False)
plt.legend()
plt.ylim(0,4)
plt.title(Title)
#plt.grid(True)
plt.show()    
