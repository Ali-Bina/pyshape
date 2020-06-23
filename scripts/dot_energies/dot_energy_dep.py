#!/usr/bin/python
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

def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)


folder = "../../"
os.chdir(folder)

config_file = 'config.ini'
config = ConfigObj(config_file, configspec='configspec.ini')
validator = Validator()
result = config.validate(validator)

f = open("scripts/power_dependence/script_results.txt","w")
f.close()
#1.06883
NITER_energy_dot_1 = 20
# low and high values for parameter to be changed (eV)
low_energy_dot_1 = 1.06883 - 0.005
high_energy_dot_1 = 1.06883 + 0.005
energy_dot_1 = linspace(low_energy_dot_1, high_energy_dot_1, num=NITER_energy_dot_1, endpoint=True)

NITER_energy_dot_2 = 20
# low and high values for parameter to be changed (eV)
low_energy_dot_2 = 1.06883 - 0.005
high_energy_dot_2 = 1.06883 + 0.005
energy_dot_2 = linspace(low_energy_dot_2, high_energy_dot_2, num=NITER_energy_dot_2, endpoint=True)


fidelity = zeros((NITER_energy_dot_1, NITER_energy_dot_2))


config['run']['optimize'] = 'True'
config['run']['show_plot'] = 'False'
config['run']['gate'] = 'twodot'

for i in range(NITER_energy_dot_1):
	config['qdot']['omega_yo'] = energy_dot_1[i]
	for j in range(NITER_energy_dot_2):
		print i, j
		config['qdot2']['omega_yo'] = energy_dot_2[j]
		config.write()
		subprocess.check_call('python main.py', shell=True)

		t = asciitable.read('data/results.txt', Reader=asciitable.CommentedHeader)
		# extract results from results.txt
		a = t['Parameter']
		fidelity_index = [item for item in range(len(a)) if a[item] == 'Overall Fidelity']
		fidelity[i, j] = t['Value'][fidelity_index]
		
config['run']['show_plot'] = 'True'
config.write()
# save data to file
savetxt('scripts/dot_energies/script_results_axes_energy_dot_1.txt', transpose(energy_dot_1), header='Dot 1 transition energy (eV)')
savetxt('scripts/dot_energies/script_results_axes_energy_dot_2.txt', transpose(energy_dot_2), header='Dot 2 transition energy (eV)')

savetxt('scripts/dot_energies/script_results.txt', transpose(fidelity), header='Fidelity')

X, Y = meshgrid(energy_dot_2, energy_dot_1)
contourf(X, Y, fidelity, linspace(0.0, 1.0, 80))
v = linspace(0.0, 1.01, 11, endpoint=True)
cbar = colorbar(ticks=v, orientation='horizontal')
cbar.ax.set_xlabel('Fidelity')
gca().xaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
gca().yaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
xlabel(r'$\hbar\omega_2$ (eV)')
ylabel(r'$\hbar\omega_1$ (eV)')
show()


