#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 14:31:28 2016

@author: ajan
"""

import numpy as np

#from pulse import *
#import scipy.interpolate
#
#from pulse import *
#from convert import *
#from read_config import *
#from phasemasks import *
#from constants import *
import matplotlib.pyplot as plt

from time import sleep
import sys


speed_of_light = 299792458.0
time_range=np.linspace(-5e-11,5e-11,1000000) #fs
lamda=1160e-9                  #Centre Wavelength of the pulse
omega_0=2*np.pi*(speed_of_light/lamda) #Angular frequency of centre Wavelength of the pulse

pulse_area=2
pulse_width=100E-15

Efield_envelope=np.sign(np.sin(2*np.pi*(time_range)/pulse_width))
Efield = Efield_envelope*(np.exp(1j*omega_0*(time_range)))
plt.plot(time_range,Efield)