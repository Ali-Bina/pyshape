#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title="M25_S4 "
filename= 'occupation_300000.0fs2_1163.002nm_0.000720.txt' 
#location2='/home/ajan/Documents/pyshape_edited/scripts/power_dependence/results/20191110/N300/'
filename2='occupation_-300000.0fs2_1163.002nm_0.000720.txt' 


var1_name=r"Pulse Area ($\pi$ units) "
var2_name=r"$\rho_{11}$"
var3_name="name"
fs=20

var=np.loadtxt(filename,skiprows = 0,) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename2,skiprows = 0) #2nd variable -yaxis1 check usecol,delimiter
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter
var3=(0.5*(var[:,3]+1))-(0.5*(var2[:,3]+1))
#v=np.where(min(var-var2, key=lambda x:abs(x-0.001)))

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
#plt.figure()
plt.plot(var[:,0],0.5*(var[:,3]+1),label=r'$\phi^{\prime \prime} >0$',linewidth=2.5,linestyle='-',color='r') #add 3rd variable
plt.plot(var[:,0],0.5*(var2[:,3]+1),label=r'$\phi^{\prime \prime} >0$',linewidth=2.5,linestyle='--',color='b')
plt.legend(loc='best',fontsize=fs-2) 
#plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name,fontsize=fs)
plt.ylabel(var2_name,fontsize=fs)
plt.ylim(0,1.1)
plt.xticks(np.linspace(0,4,5))
#plt.title(Title)
plt.rc('xtick',labelsize=fs-2)
plt.rc('ytick',labelsize=fs-2)

#plt.grid(True)
plt.show()    
