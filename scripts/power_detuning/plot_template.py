#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 10:32:34 2014

@author: ajan
"""

from pylab import *


Title="ARP "
filename1= 'script_results_axes_detuning.txt' 
filename2= 'script_results_axes_pulse_area.txt' 
filename3= 'occupation.txt' 



pulse_area=loadtxt(filename1,skiprows = 1,) #1st variable -x axis; usecols=(0,1)
detun=loadtxt(filename2,skiprows = 1,)#var2=np.loadtxt(filename,usecols=(1,),skiprows = 0, unpack=True) #2nd variable -yaxis1 check usecol,delimiter
occupation=loadtxt(filename3,skiprows = 0,)
"""
#Delay position to ps
speed_of_light = 299792458.0
zero_delay = -13800.0  #zero_delay
var[:,0] = -1E12*2*(var[:,0] - zero_delay)*1E-6/(speed_of_light)
"""

X, Y = meshgrid(pulse_area,detun)
contourf(X, Y, occupation, np.linspace(0.0, 1.0, 80))
v = linspace(0.0, 1.0, 11, endpoint=True)
cbar = colorbar(ticks=v, orientation='horizontal')
cbar.ax.set_xlabel('Occupation')
rotation='horizontal'
xlabel(r'Detuning $\Delta$ (meV)')
ylabel(r'Pulse Area $\Theta$ ($\pi$ rad)')
show()

