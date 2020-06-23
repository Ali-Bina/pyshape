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
import glob

data_files=glob.glob("*.txt")
data_files.sort()

Title=" $\omega_0$ =1174.393 nm, $\omega_L$ = 1174.393 nm  "
#filename= 'occupation_30000.0fs2_1.07514eV.txt' 
#filename2= 'tarray.txt' 

#chirp=float(data_files[0].rsplit('_')[1].rstrip("fs2"))
chirp=np.zeros(len(data_files))

var1_name="Pulse Area ($\pi$ Units)"
var2_name="Occupation"
var3_name="name"

plt.figure()
for i in range(len(data_files)):
    chirp[i]=float(data_files[i].rsplit('_')[1].rstrip("fs2"))

    var=np.loadtxt(data_files[i],usecols=(0,3),skiprows = 0,) #1st variable -x axis; usecols=(0,1)
    plt.plot(var[:,0],var[:,1],label=str(chirp[i])+"fs2") #add 3rd variable 
    #plt.axis([1200, 1320, 0, 1.1])
plt.xlabel(var1_name)
plt.ylabel(var2_name)
plt.legend(bbox_to_anchor=(0.5, -0.1), loc=9, ncol=3)
plt.ylim(0,1.1)
plt.title(Title)
plt.grid(True)
plt.show()  
plt.savefig("Dot3_arp.svg")  

#for i in range(len(data_files)):
#    plt.figure()
#    wl[i]=float(data_files[i].rsplit('_')[2].rstrip("nm.txt"))
#
#    var=np.loadtxt(data_files[i],usecols=(0,3),skiprows = 0,) #1st variable -x axis; usecols=(0,1)
#    plt.plot(var[:,0],var[:,1],label=str(wl[i])+"nm") #add 3rd variable 
#    #plt.axis([1200, 1320, 0, 1.1])
#    plt.xlabel(var1_name)
#    plt.ylabel(var2_name)
#    plt.legend(bbox_to_anchor=(0.5, -0.1), loc=9, ncol=4)
#    plt.ylim(0,1.1)
#    plt.title(Title)
#    plt.grid(True)
#    plt.savefig("Dot2_TL"+str(wl[i])+"nm.svg") 
#    plt.close()