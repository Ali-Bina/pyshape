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

pulse_width_low = 100
pulse_width_high = 200
tau = linspace(pulse_width_low, pulse_width_high, NITER)


config['run']['optimize'] = 'False'

for i in range(NITER):
	print "THIS IS RUN ", i
	config['pulse']['width'] = pulse_width_low + i*(pulse_width_high - pulse_width_low)/NITER
	config.write()
	#subprocess.check_call('python main.py >> errors.txt', shell=True)
	subprocess.check_call('python main.py', shell=True)


f = open("rawdata.txt","rb")
for line in f:
	print line

f.close()

EO, fidelity1, fidelity2  = loadtxt("rawdata.txt", delimiter = " ", unpack=True)

plot(tau, fidelity1, 'r', label=r'$\rho_{11}$ dot 1')
plot(tau, 1-fidelity2, 'g', label=r'$\rho_{11}$ dot 2')
plot(tau, fidelity1*fidelity2, 'b--',  label='fidelity')

gca().legend()
xlabel(r'Pulse Width (fs)')
ylabel(r'Fidelity')
grid(True)
show()

