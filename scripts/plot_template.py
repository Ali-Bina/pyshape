#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:52:57 2016

@author: ajan
"""

#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

import numpy as np
import matplotlib.pyplot as plt 

Title="M25_S4 "
filename= 'occupation.txt' 
#filename2= 'tarray.txt' 

var1_name="name"
var2_name="name"
var3_name="name"

var=np.loadtxt(filename,usecols=(0,3),skiprows = 0,) #1st variable -x axis; usecols=(0,1)
#var2=np.linspace(0,100,num=len(var[:,2]))
#var2=np.loadtxt(filename,usecols=(1,),skiprows = 0, unpack=True) #2nd variable -yaxis1 check usecol,delimiter
#var3=np.loadtxt(filename,delimiter=',',usecols=(2,),skiprows = 0, unpack=True) #3rd variable -yaxis2 check usecol,delimiter

"""
#Delay position to ps
speed_of_light = 299792458.0
zero_delay = -13800.0  #zero_delay
var[:,0] = -1E12*2*(var[:,0] - zero_delay)*1E-6/(speed_of_light)
"""
plt.figure()
plt.plot(var[:,0],var[:,1],'r') #add 3rd variable 
#plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name)
plt.ylabel(var2_name)
#plt.xlim(0,1024)
plt.title(Title)
plt.grid(True)
plt.show()    
