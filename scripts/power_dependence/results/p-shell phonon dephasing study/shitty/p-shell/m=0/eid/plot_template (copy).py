#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title=""
filename= '/home/ajan/Documents/pyshape_edited/data/eid/K_real.txt' 
filename2= '/home/ajan/Documents/pyshape_edited/data/eid/Omega_eid.txt' 

var1_name="Omega_eid"
var2_name="K_real"
var3_name="name"

var=np.loadtxt(filename,skiprows = 1,) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename2,skiprows = 1) #2nd variable -yaxis1 check usecol,delimiter
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

plt.plot(var2,var,'r') #add 3rd variable 
#plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name)
plt.ylabel(var2_name)
#plt.xlim(0,1024)
plt.title(Title)
plt.grid(True)
plt.show()    
