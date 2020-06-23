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

folder = "../../"
os.chdir(folder)

config_file = 'config.ini'
config = ConfigObj(config_file, configspec='configspec.ini')
validator = Validator()
result = config.validate(validator)
print config['masks']['phasemask_cos']['pulse_area']


f = open("scripts/power_dependence/script_results.txt","w")
f.close()


# number of iterations
NITER = 4

# low and high values for parameter to be changed
low = 0.1
high = 3.5
delta = (high - low)/NITER

pulse_area = zeros(NITER)
fidelity1 = zeros(NITER)
fidelity2 = zeros(NITER)
overallFid = zeros(NITER)

config['run']['optimize'] = 'False'
config['run']['show_plot'] = 'False'
config['run']['gate'] = 'threedot'


for i in range(NITER):
	print "THIS IS RUN ", i
	config['masks']['phasemask_cos']['pulse_area'] = low + i*delta
	config.write()
	pulse_area[i] = low + i*delta
	subprocess.check_call('python main.py', shell=True)
	
	t = asciitable.read('data/results.txt', Reader=asciitable.CommentedHeader)
	# extract results from results.txt
	a = t['Parameter']
	fid1_index = [item for item in range(len(a)) if a[item] == 'Fidelity DOT 1']
	fidelity1[i] = t['Value'][fid1_index]
	fid2_index = [item for item in range(len(a)) if a[item] == 'Fidelity DOT 2']
	fidelity2[i] = t['Value'][fid2_index]
	fid3_index = [item for item in range(len(a)) if a[item] == 'Fidelity DOT 3']
	fidelity3[i] = t['Value'][fid3_index]
	overallfid_index = [item for item in range(len(a)) if a[item] == 'Overall Fidelity']
	overallFid[i] = t['Value'][overallfid_index]

config['run']['show_plot'] = 'True'

# save data to ifile
savetxt('scripts/power_dependence/script_results.txt', transpose((pulse_area, fidelity1, fidelity2, fidelity3, overallFid)), header='Pulse Area (PI rad), Fidelity DOT 1, Fidelity DOT 2, Fidelity DOT 3, Overall Fidelity')

# load data and plot
pulse_area, fidelity1, fidelity2, fidelity3, overall_fid = loadtxt("scripts/power_dependence/script_results.txt", delimiter = " ", unpack=True, skiprows=1)
plot(pulse_area, fidelity1, 'r', label=r'Fidelity dot 1')
plot(pulse_area, fidelity2, 'g', label=r'Fidelity dot 2')
plot(pulse_area, fidelity3, 'k', label=r'Fidelity dot 3')
plot(pulse_area, overall_fid, 'b--',  label='Overall fidelity')

gca().legend(loc=4)
xlabel(r'Pulse Area (PI radians)')
ylabel(r'Fidelity')
grid(True)
show()

