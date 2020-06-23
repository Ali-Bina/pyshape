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
import gc

folder = "../../"
os.chdir(folder)

config_file = 'config.ini'
config = ConfigObj(config_file, configspec='configspec.ini')
validator = Validator()
result = config.validate(validator)

f = open("scripts/power_dependence/script_results.txt","w")
f.close()


# number of iterations
NITER = 40

# low and high values for parameter to be changed
low = 0.0
high = 4.8
delta = (high - low)/(NITER-1)


pulse_area = zeros(NITER)
fidelity = zeros(NITER)
mask_name = 'phasemask_poly'
#config['qdot']['eid_c'] = 150.0
#config['laphononeid']['omega_c_eid'] = 0.010
#config['masks'][mask_name]['phi_2'] = 132934.0
#config['pulse']['width'] = 396.0
#config['run']['optimize'] = 'False'
#config['run']['show_plot'] = 'True'
#config['run']['gate'] = 'onedot'
#config['run']['mask'] = mask_name
if mask_name == 'phasemask_poly':
    chirp=config['masks'][mask_name]['phi_2']
elif mask_name =='phasemask_none':
    chirp="TL"

if config['run']['gate'] == 'twolevel':
    occupation=zeros((NITER,4))
else:
    occupation=zeros((NITER,17))


pulse_omega=config['pulse']['omega_o']
pulse_wl=str(round(1239.84193/pulse_omega,3))
omega_eid=config['laphononeid']['omega_c_eid']
dot_reson=config['qdot']['omega_yo']

for i in range(NITER):
	print "THIS IS RUN ", i ,"goes from", low,"to",high
	config['masks'][mask_name]['pulse_area'] = low + i*delta
	config.write()
	pulse_area[i] = low + i*delta
	print pulse_area[i]," wl="+str(pulse_wl)+"nm"
	 
	subprocess.check_call('python main.py', shell=True)
	
	t = asciitable.read('data/results.txt', Reader=asciitable.CommentedHeader)
	# extract results from results.txt
	a = t['Parameter']
	fidelity_index = [item for item in range(len(a)) if a[item] == 'Fidelity']
	fidelity[i] = t['Value'][fidelity_index]
	occupation[i,1:]=loadtxt('data/occupation.txt')
	occupation[i,0]=pulse_area[i]
#config['run']['show_plot'] = 'True'
config.write()
# save data to ifile
savetxt('scripts/power_dependence/script_results.txt', transpose((pulse_area, fidelity)), header='Pulse Area (PI rad), Fidelity')
if mask_name == 'phasemask_poly':
    savetxt('scripts/power_dependence/results/occupation_'+'%08.1f' %chirp+'fs2_'+pulse_wl+'nm_'+'%08f' %dot_reson+".txt", occupation)
else:
    savetxt('scripts/power_dependence/results/occupation_'+chirp+'fs2_'+pulse_wl+"nm.txt", occupation)
# load data and plot
pulse_area, fidelity = loadtxt("scripts/power_dependence/script_results.txt", delimiter = " ", unpack=True, skiprows=1)

#plot(pulse_area, fidelity, 'b-o',  label='fidelity')
#ylim(0.0, 1.0)
#xlabel(r'Pulse Area (PI radians)')
#ylabel(r'Fidelity')
#grid(True)
#show(block=False)
#close()
gc.collect()
