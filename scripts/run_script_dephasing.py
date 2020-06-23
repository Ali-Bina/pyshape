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

NITER = 20

dephasing_low = 0.5
dephasing_high = 10.0
dephasing = linspace(dephasing_low, dephasing_high, NITER)


config['run']['optimize'] = 'False'

for i in range(NITER):
	print "THIS IS RUN ", i
	config['qdot']['T_yo'] = dephasing_low + i*(dephasing_high - dephasing_low)/NITER
	config.write()
	#subprocess.check_call('python main.py >> errors.txt', shell=True)
	subprocess.check_call('python main.py', shell=True)


f = open("rawdata.txt","rb")
for line in f:
	print line

f.close()

EO, fidelity1, fidelity2  = loadtxt("rawdata.txt", delimiter = " ", unpack=True)

plot(dephasing, fidelity1, 'r', label=r'$\rho_{11}$ dot 1')
plot(dephasing, 1-fidelity2, 'g', label=r'$\rho_{11}$ dot 2')
plot(dephasing, fidelity1*fidelity2, 'b--',  label='fidelity')

gca().legend()
xlabel(r'Pulse Phase (PI radians)')
ylabel(r'Fidelity')
grid(True)
show()

