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

NITER = 100

cos_phase_low = 0.0
cos_phase_high = 0.5*pi
phase = linspace(cos_phase_low, cos_phase_high, NITER)


config['run']['optimize'] = 'False'

for i in range(NITER):
	print "THIS IS RUN ", i
	config['masks']['phasemask_cos']['phi_delta'] = cos_phase_low + i*(cos_phase_high - cos_phase_low)/NITER
	config.write()
	#subprocess.check_call('python main.py >> errors.txt', shell=True)
	subprocess.check_call('python main.py', shell=True)


f = open("rawdata.txt","rb")
for line in f:
	print line

f.close()

EO, fidelity1, fidelity2  = loadtxt("rawdata.txt", delimiter = " ", unpack=True)

plot(phase/(pi), fidelity1, 'r', label=r'$\rho_{11}$ dot 1')
plot(phase/(pi), 1-fidelity2, 'g', label=r'$\rho_{11}$ dot 2')
plot(phase/(pi), fidelity1*fidelity2, 'b--',  label='fidelity')

gca().legend()
xlabel(r'Pulse Phase (PI radians)')
ylabel(r'Fidelity')
grid(True)
show()

