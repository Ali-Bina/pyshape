#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title="M25_S4 "
filename= 'occupation_500000.0fs2_1162.861nm.txt' 
filename2="occupation_-500000.0fs2_1162.861nm.txt"

l1=filename.rsplit("_")[1]
l2=filename2.rsplit("_")[1]
var1_name="Pulse Area"
var2_name="Counts"
var3_name="name"

var=0.5*(np.loadtxt(filename,usecols=(3,),skiprows = 0,)+1) #1st variable -x axis; usecols=(0,1)
var2=0.5*(np.loadtxt(filename2,usecols=(3,),skiprows = 0,)+1) #2nd variable -yaxis1 check usecol,delimiter
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter
pulse_area=np.loadtxt(filename,usecols=(0,),skiprows = 0,)

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
plt.plot(pulse_area,var,label=l1)
plt.plot(pulse_area,var2,label=l2) 
#plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name)
plt.ylabel(var2_name)
plt.ylim(0,1.1)
plt.legend(loc=4)
#plt.title(Title)
#plt.grid(True)
plt.show()    
plt.savefig(l1+".svg")
