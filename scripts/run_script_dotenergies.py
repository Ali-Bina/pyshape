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


config_file = 'config.ini'
config = ConfigObj(config_file, configspec='configspec.ini')
validator = Validator()
result = config.validate(validator)

f = open("rawdata.txt","w")
f.close()

NITER = 4

qdot_3 = 1.07147
qdot_3_low = qdot_3 - 0.0025
qdot_3_high = qdot_3 + 0.0025

qdot_2 = 1.07459
qdot_2_low = qdot_2 - 0.0025
qdot_2_high = qdot_2 + 0.0025


config['run']['optimize'] = 'False'
for i in range(NITER):
	config['qdot']['omega_yo'] = qdot_3_low + i*(qdot_3_high - qdot_3_low)/NITER
	for j in range(NITER):
		print i, j
		config['qdot2']['omega_yo'] = qdot_2_low + j*(qdot_2_high - qdot_2_low)/NITER
		config.write()
		subprocess.check_call('python main.py', shell=True)


config['qdot']['omega_yo'] = qdot_3
config['qdot2']['omega_yo'] = qdot_2
config.write()

f = open("rawdata.txt","rb")
for line in f:
	print line

f.close()

a, f1, f2 = loadtxt("rawdata.txt", delimiter = " ", unpack=True)
print f1
fidelity = (1.0-f1)*(1.0-f2)
x = 1000.0*(linspace(qdot_3_low, qdot_3_high, NITER) - qdot_3)
y = 1000.0*(linspace(qdot_2_low, qdot_2_high, NITER) - qdot_2)
X, Y = meshgrid(y, x)

f = fidelity.reshape((NITER, NITER))

contourf(X, Y, f, linspace(0.0, 1.0, 80))
cbar = colorbar(orientation='horizontal')
cbar.ax.set_ylabel('Fidelity')
xlabel(r'Dot 2 Energy (meV)')
ylabel(r'Dot 3 Energy (meV)')
show()

