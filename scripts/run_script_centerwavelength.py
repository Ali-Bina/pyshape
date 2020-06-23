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
from constants import *

config_file = 'config.ini'
config = ConfigObj(config_file, configspec='configspec.ini')
validator = Validator()
result = config.validate(validator)

f = open("rawdata.txt","w")
f.close()

NITER = 100

pulse_center_low = 1.0642
pulse_center_high = 1.0828
center = linspace(pulse_center_low, pulse_center_high, NITER)


config['run']['optimize'] = 'False'

for i in range(NITER):
	print "THIS IS RUN ", i
	config['pulse']['omega_o'] = pulse_center_low + i*(pulse_center_high - pulse_center_low)/NITER
	config.write()
	#subprocess.check_call('python main.py >> errors.txt', shell=True)
	subprocess.check_call('python main.py', shell=True)


f = open("rawdata.txt","rb")
for line in f:
	print line

f.close()

EO, fidelity1, fidelity2  = loadtxt("rawdata.txt", delimiter = " ", unpack=True)
center_lambda = 1E9*(2*pi*H_BAR*C)/(center*E_CHARGE)
plot(center_lambda, fidelity1, 'r', label=r'$\rho_{11}$ dot 1')
plot(center_lambda, 1-fidelity2, 'g', label=r'$\rho_{11}$ dot 2')
plot(center_lambda, fidelity1*fidelity2, 'b--',  label='fidelity')

gca().legend()
xlabel(r'Pulse Center (nm)')
ylabel(r'Fidelity')
grid(True)
show()

