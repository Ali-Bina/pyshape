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

omega_c= np.array([0.561,0.612,0.673,0.747,0.841,0.961,1.121,1.345,1.682])
omega_d= np.array([0.561,0.612,0.673,0.747,0.841,0.961,1.121,1.345,1.682])

Theta_T=np.loadtxt('omega_grid_1.txt')


plt.figure()
A,B=np.meshgrid(omega_c,omega_d)
plt.contourf(A, B,np.transpose(Theta_T),100)
#plt.clim(np.amin(BX+Y),np.amax(BX+Y))
plt.colorbar(pad=0.12)
plt.xlabel(r'$\omega_d$')
plt.ylabel(r'$\omega_L$')
plt.show()