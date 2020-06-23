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

f = open("test.txt","w")
f.close()

NITER = 5tom

qdot_1_dipole = 38.006
qdot_1_low = qdot_1_dipole - 20
qdot_1_high = qdot_1_dipole + 20

qdot_2_dipole = 45.791
qdot_2_low = qdot_2_dipole - 20
qdot_2_high = qdot_2_dipole + 20

param = 'd_yo'
config['run']['optimize'] = 'False'
for i in range(NITER):
	config['qdot'][param] = qdot_1_low + i*(qdot_1_high - qdot_1_low)/NITER
	for j in range(NITER):
		print i, j
		config['qdot2'][param] = qdot_2_low + j*(qdot_2_high - qdot_2_low)/NITER
		config.write()
		#subprocess.check_call('python main.py >> errors.txt', shell=True)
		subprocess.check_call('python main.py', shell=True)



config['qdot']['d_yo'] = qdot_1_dipole
config['qdot2']['d_yo'] = qdot_2_dipole
config.write()

f = open("test.txt","rb")
for line in f:
	print line

f.close()

a, b, fidelity = loadtxt("test.txt", delimiter = " ", unpack=True)

x = linspace(qdot_1_low, qdot_1_high, NITER)
y = linspace(qdot_2_low, qdot_2_high, NITER)
X, Y = meshgrid(x, y)
print x
print y
f = fidelity.reshape((NITER, NITER))
print f
contourf(X, Y, f, linspace(0.0, 1.0, 80))
cbar = colorbar(orientation='horizontal')
cbar.ax.set_ylabel('Fidelity')
xlabel(r'Dot 1 Dipole Moment (Debye)')
ylabel(r'Dot 2 Dipole Moment (Debye)')
show()

