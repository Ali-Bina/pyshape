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
import gc


config_file = '../../config.ini'
config = ConfigObj(config_file, configspec='../../configspec.ini')
validator = Validator()
result = config.validate(validator)

mask_name = 'phasemask_poly'  #phasemask_poly for ARP  phasemask_none for TL
#apply_chirp=0    #chirp to apply in fs2


center_wl=array([1163.0])
#array([1148.0,1153.0,1156.0,1159.0,1163.0,1166.4,1170.0,1174.4,1178.0,1181.0])#array([1171.0])#linspace(1167,1175,num=20)#array([1148.0,1153.0,1156.0,1159.0,1163.0,1166.4,1170.0,1174.4,1178.0,1181.0])
center_omega=1239.84193/center_wl
chirp=array([300000.0,-300000.0])
omega_eid=array([0.000841,0.000961,0.001121,0.001345,0.001682])
#array([0.000561,0.000612,0.000673,0.000747,0.000841,0.000961,0.001121,0.001345,0.001682])
#array([0.00056,0.00071,0.0009,0.001,0.00115,0.0013,0.00145,0.0016])

dot_reson=config['qdot']['omega_yo']
#config['run']['phasemask'] = mask_name
#config['masks'][mask_name]['phi_2']=apply_chirp
config.write()
#if mask_name == 'phasemask_poly':
#    chirp=str(config['masks'][mask_name]['phi_2'])
#elif mask_name =='phasemask_none':
#    chirp="TL"

#for i in range(len(center_omega)):
#    print "Run parameters: omega_o="+str(center_omega[i])+" mask="+mask_name+" chirp="+str(chirp)
#    config['pulse']['omega_o']=round(center_omega[i],5)
#    config['pulse']['detun']=round(abs(round(center_omega[i],5)-dot_reson),5)
#    config.write()
#    
#    subprocess.check_call('python power_dep_one_dot.py', shell=True)

#
#for i in range(len(center_omega)):
#    config['pulse']['omega_o']=round(center_omega[i],5)
#    config['pulse']['detun']=round(abs(round(center_omega[i],5)-dot_reson),5)
#    config.write()
#    for j in range(len(chirp)):
#        config['masks'][mask_name]['phi_2']=chirp[j]
#        config.write()
#        print "Run parameters: omega_o="+str(center_omega[i])+" mask="+mask_name+" chirp="+str(chirp[j])
#
#
#
#        subprocess.check_call('python power_dep_one_dot.py', shell=True)
#        gc.collect()
for i in range(len(omega_eid)):
    config['laphononeid']['omega_c_eid']=omega_eid[i]
    config['pulse']['omega_o']=round(center_omega[0],5)
    config['pulse']['detun']=round(abs(round(center_omega[0],5)-dot_reson),5)
    config.write()
    for j in range(len(chirp)):
        config['masks'][mask_name]['phi_2']=chirp[j]
        config.write()
        print "Run parameters: omega_o="+str(center_omega)+" mask="+mask_name+" chirp="+str(chirp[j])+" omega_eid= "+str(omega_eid[i])
#
#
#
        subprocess.check_call('python power_dep_one_dot.py', shell=True)
        gc.collect()
