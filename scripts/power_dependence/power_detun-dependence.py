#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 17:07:52 2016

@author: ajan
"""

#! /home/ajan/anaconda/bin/python
# Filename: run_script.py


import subprocess
from configobj import ConfigObj
from validate import Validator
import cPickle
import glob				# used to list directory contents
import os				# used for os dependent functionality
import operator		# used for sorting lists
from pylab import *
import csv
import asciitable
import numpy



config_file = '../../config.ini'
config = ConfigObj(config_file, configspec='../../configspec.ini')
validator = Validator()
result = config.validate(validator)

mask_name = 'phasemask_poly'  #phasemask_poly for ARP  phasemask_none for TL
apply_chirp=0    #chirp to apply in fs2


center_wl=array([1153.192,1160.540,1163.4,1166.513,1170.981,1174.396,1179.310])
center_omega=1239.84193/center_wl

dot_reson=config['qdot']['omega_yo']
config['run']['phasemask'] = mask_name
config['masks'][mask_name]['phi_2']=apply_chirp
config.write()
if mask_name == 'phasemask_poly':
    chirp=str(config['masks'][mask_name]['phi_2'])
elif mask_name =='phasemask_none':
    chirp="TL"

for i in range(len(center_omega)):
    print "Run parameters: omega_o="+str(center_omega[i])+" mask="+mask_name+" chirp="+str(chirp)
    config['pulse']['omega_o']=round(center_omega[i],5)
    config['pulse']['detun']=round(abs(round(center_omega[i],5)-dot_reson),5)
    config.write()
    
    subprocess.check_call('python power_dep_one_dot.py', shell=True)