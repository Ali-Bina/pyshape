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

NITER_power = 50
# low and high values for parameter to be changed
low_power = 0.1
high_power = 5.0
pulse_area = linspace(low_power, high_power, num=NITER_power, endpoint=True)

NITER_chirp = 50
# low and high values for parameter to be changed
low_chirp = -10000000.0
high_chirp = 10000000.0
chirp = linspace(low_chirp, high_chirp, num=NITER_chirp, endpoint=True)

fidelity = zeros((NITER_power, NITER_chirp))


config['run']['optimize'] = 'False'
config['run']['show_plot'] = 'False'
config['run']['gate'] = 'twolevel'

config['run']['optimize'] = 'False'
for i in range(NITER_power):
	config['masks']['phasemask_poly']['pulse_area'] = pulse_area[i]
	for j in range(NITER_chirp):
		print i, j
		config['masks']['phasemask_poly']['phi_2'] = chirp[j]
		config.write()
		subprocess.check_call('python main.py', shell=True)

		t = asciitable.read('data/results.txt', Reader=asciitable.CommentedHeader)
		# extract results from results.txt
		a = t['Parameter']
		fidelity_index = [item for item in range(len(a)) if a[item] == 'Fidelity']
		fidelity[i, j] = t['Value'][fidelity_index]
		
config['run']['show_plot'] = 'True'
config.write()
# save data to file
savetxt('scripts/power_chirp/script_results_axes_pulse_area.txt', transpose(pulse_area), header='Pulse Area (PI rad)')
savetxt('scripts/power_chirp/script_results_axes_chirp.txt', transpose(chirp), header='Chirp (fs2)')

savetxt('scripts/power_chirp/script_results.txt', transpose(fidelity), header='Fidelity')

X, Y = meshgrid(chirp/1000000.0, pulse_area)
contourf(X, Y, fidelity, linspace(0.0, 1.0, 80))
v = linspace(0.0, 1.0, 11, endpoint=True)
cbar = colorbar(ticks=v, orientation='horizontal')
cbar.ax.set_xlabel('Fidelity')
xlabel(r'Chirp $\phi_2$ (ps$^2$)')
ylabel(r'Pulse Area $\Theta$ ($\pi$ rad)')
show()



