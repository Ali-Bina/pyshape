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

#Title="$\omega_0$ =1180.802 nm, "
#filename= 'occupation_30000.0fs2_1.07514eV.txt' 
#filename2= 'tarray.txt' 

#chirp=float(data_files[0].rsplit('_')[1].rstrip("fs2"))
wl=np.zeros(len(data_files))
chirp=np.zeros(len(data_files))
pulse=[]

var1_name="Pulse Area ($\pi$ Units)"
var2_name="Occupation"
var3_name="name"

#plt.figure()
#for i in range(len(data_files)):
#    chirp[i]=float(data_files[i].rsplit('_')[1].rstrip("fs2"))
#    pulse.append(data_files[i].rsplit('_')[3].rstrip(".txt"))
#
#    var=np.loadtxt(data_files[i],usecols=(0,3),skiprows = 0,) #1st variable -x axis; usecols=(0,1)
#    plt.plot(var[:,0],var[:,1],label=pulse[i]) #add 3rd variable 
#    #plt.axis([1200, 1320, 0, 1.1])
#plt.xlabel(var1_name)
#plt.ylabel(var2_name)
#plt.legend(bbox_to_anchor=(0.5, -0.1), loc=9, ncol=4)
#plt.ylim(0,1.1)
##plt.title(Title)
#plt.grid(True)
#plt.show()  
#plt.savefig("Dot3_arp.svg")  

for i in range(len(data_files)):
    plt.figure()
    wl[i]=float(data_files[i].rsplit('_')[2].rstrip("nm.txt"))

    var=np.loadtxt(data_files[i],usecols=(0,2,3,4),skiprows = 0,) #1st variable -x axis; usecols=(0,1)
    plt.plot(var[:,0],var[:,1],label="XX") #add 3rd variable 
    plt.plot(var[:,0],var[:,2],label="YY")
    plt.plot(var[:,0],var[:,3],label="BX")
    plt.plot(var[:,0],var[:,1]+var[:,2]+var[:,3],label="sum")
    #plt.axis([1200, 1320, 0, 1.1])
    plt.xlabel(var1_name)
    plt.ylabel(var2_name)
    plt.legend(bbox_to_anchor=(0.5, -0.1), loc=9, ncol=4)
    plt.ylim(-0.1,1.1)
    plt.title(wl[i])
    plt.grid(True)
    plt.savefig("Dot2_"+str(wl[i])+"nm.svg") 
    plt.show()