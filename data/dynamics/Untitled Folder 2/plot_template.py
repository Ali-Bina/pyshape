#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title="M25_S4 "
filename= 'Efield_intensity.txt'
filename2='time.txt'
filename3="detuning.txt"

plt.style.use('physrev')

var1_name="time (fs)"
var2_name="Intensity"
var3_name="Detuning (meV)"

var=np.genfromtxt(filename) #1st variable -x axis; usecols=(0,1)
var2=np.loadtxt(filename2,skiprows = 0) #2nd variable -yaxis1 check usecol,delimiter
var3=np.loadtxt(filename3,skiprows = 0)
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter
#Int=(np.absolute(var))**2
#time= var2-25000
#var2=1240.0/var2
#y=[]
x=np.absolute(var-0.5)
id1=np.where(x==x[x>0].min())[0]
id2=np.where(x==x[x>x[x>0].min()].min())[0]
print var3[id1],var3[id2],np.absolute(var3[id1]-var3[id2])
#y.append(var3[id1]-var3[id2])
"""
#Delay position to ps
speed_of_light = 299792458.0
zero_delay = -13800.0  #zero_delay
var[:,0] = -1E12*2*(var[:,0] - zero_delay)*1E-6/(speed_of_light)
"""
fig, ax1 = plt.subplots()
plt.plot(var2,var,'g') #add 3rd variable 
#plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name)
plt.ylabel(var2_name)
#plt.xlim(-4000,4000)
#plt.ylim(0,max(Int))
#plt.title(Title)
#plt.grid(True)

hbar=1.0546E-34
e=1.6022E-19
chirp_s=10000.0
chirp_t=chirp_s/(chirp_s**2+120.0**4)
meV_conv=2*np.pi*1E15*hbar*1000/e
#
#ax2=ax1.twinx()                          
#ax2.plot(time,var3,'c')
##ax2.plot(time,var3,'b')
#ax2.set_ylim(-50,50)
##ax2.set_xlim(-4000,4000)
#ax2.set_ylabel(var3_name)
#plt.show()    
