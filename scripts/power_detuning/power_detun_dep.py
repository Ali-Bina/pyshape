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

f = open("scripts/power_detuning/script_results.txt","w")
f.close()

# number of iterations
NITER_power = 11
# low and high values for parameter to be changed
low_power = 0.01
high_power = 5.0
pulse_area = linspace(low_power, high_power, num=NITER_power, endpoint=True)

NITER_detun =11
# low and high values for parameter to be changed
low_detun = -2.0
high_detun = 2.0
detun = linspace(low_detun, high_detun, num=NITER_detun, endpoint=True)
print detun
fidelity = zeros((NITER_power, NITER_detun))

occupation=zeros((NITER_power, NITER_detun))

config['run']['optimize'] = 'False'
config['run']['show_plot'] = 'False'
config['run']['gate'] = 'onedot'
config['run']['optimize'] = 'False'

omega_o = config['pulse']['omega_o']

for i in range(NITER_power):
	config['masks']['phasemask_poly']['pulse_area'] = pulse_area[i]
	for j in range(NITER_detun):
		print i, j
		config['pulse']['omega_o'] = omega_o + detun[j]/1000.0
		config.write()
		subprocess.check_call('python main.py', shell=True)

		t = asciitable.read('data/results.txt', Reader=asciitable.CommentedHeader)
		# extract results from results.txt
		a = t['Parameter']
		fidelity_index = [item for item in range(len(a)) if a[item] == 'Fidelity']
		fidelity[i, j] = t['Value'][fidelity_index]
		
		occupation[i,j]=loadtxt('data/occupation.txt',usecols=(2,))
		print  detun[j], pulse_area[i], fidelity[i,j]
		
config['run']['show_plot'] = 'True'
config['pulse']['omega_o'] = omega_o
config.write()

# save data to file
savetxt('scripts/power_detuning/script_results_axes_pulse_area.txt', transpose(pulse_area), header='Pulse Area (PI rad)')
savetxt('scripts/power_detuning/script_results_axes_detuning.txt', transpose(detun), header='Detuning (meV)')
savetxt('scripts/power_detuning/occupation.txt', occupation)
savetxt('scripts/power_detuning/script_results.txt', transpose(fidelity), header='Fidelity')

X, Y = meshgrid(detun, pulse_area)
contourf(X, Y, occupation, linspace(0.0, 1.0, 80))
v = linspace(0.0, 1.0, 11, endpoint=True)
cbar = colorbar(ticks=v, orientation='horizontal')
cbar.ax.set_xlabel('Occupation')
rotation='horizontal'
xlabel(r'Detuning $\Delta$ (meV)')
ylabel(r'Pulse Area $\Theta$ ($\pi$ rad)')
show()



