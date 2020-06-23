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

data_folder = 'scripts/two_dot_energy_dipole/'										#***********

# define parameters of dot1
energy_dot_1_init = 1.07
dipole_dot_1_init = 25.0
config['qdot']['omega_yo'] = energy_dot_1_init
config['qdot']['d_yo'] = dipole_dot_1_init

config['qdot']['initial_state_co_r'] = 1.0	# modulus r of state o
config['qdot']['initial_state_co_phi'] = 0.0	# phase phi of state o (as a fraction of PI)
config['qdot']['initial_state_cy_r'] = 1.0	# modulus r of state o
config['qdot']['initial_state_cy_phi'] = 0.0	# phase phi of state o (as a fraction of PI)

config['qdot']['desired_state_co_r'] = 1.0	# modulus r of state o
config['qdot']['desired_state_co_phi'] = 0.0	# phase phi of state o (as a fraction of PI)
config['qdot']['desired_state_cy_r'] = 0.0	# modulus r of state o
config['qdot']['desired_state_cy_phi'] = 0.0	# phase phi of state o (as a fraction of PI)

# define parameters of dot2
energy_dot_2_init = 1.071
dipole_dot_2_init = 27.0
config['qdot2']['omega_yo'] = energy_dot_2_init
config['qdot2']['d_yo'] = dipole_dot_2_init

config['qdot2']['initial_state_co_r'] = 1.0	# modulus r of state o
config['qdot2']['initial_state_co_phi'] = 0.0	# phase phi of state o (as a fraction of PI)
config['qdot2']['initial_state_cy_r'] = 1.0	# modulus r of state o
config['qdot2']['initial_state_cy_phi'] = 0.0	# phase phi of state o (as a fraction of PI)

config['qdot2']['desired_state_co_r'] = 1.0	# modulus r of state o
config['qdot2']['desired_state_co_phi'] = 0.0	# phase phi of state o (as a fraction of PI)
config['qdot2']['desired_state_cy_r'] = 0.0	# modulus r of state o
config['qdot2']['desired_state_cy_phi'] = 0.0	# phase phi of state o (as a fraction of PI)

# define pulse parameters
config['run']['optimize'] = 'True'
config['run']['show_plot'] = 'False'
config['run']['gate'] = 'twodot'

config.write()

# set up axes
var1_init = 0.0														#***********
var1_delta = 1.0																		#***********
NITER_var1 = 61		# make sure NITER is odd to get the value at 0		#***********
var2_init = 0.0														#***********
var2_delta = 1.0																		#***********
NITER_var2 = 61		# make sure NITER is odd to get the value at 0		#***********

# low and high values for parameter to be changed 
low_var1 = var1_init - var1_delta
high_var1 = var1_init + var1_delta
var1 = linspace(low_var1, high_var1, num=NITER_var1, endpoint=True)

# low and high values for parameter to be changed
low_var2 = var2_init - var2_delta
high_var2 = var2_init + var2_delta
var2 = linspace(low_var2, high_var2, num=NITER_var2, endpoint=True)

fidelity = zeros((NITER_var1, NITER_var2))

# save x and y axis values, and results to file
savetxt(data_folder + 'script_results_axes_w_diff.txt', transpose(var1), header='w2 - w1')		#***********
savetxt(data_folder + 'script_results_axes_phase_diff.txt', transpose(var2), header='Phase Difference (PI rad)')	#***********
savetxt(data_folder + 'script_results.txt', transpose(fidelity), header='Fidelity')


i_initial = 0							#***********
i_final = 2					#***********
j_initial = 0							#***********
j_final = NITER_var2					#***********
for i in range(i_initial, i_final):
	w1 = -var1[i]/2.0								#***********
	co_r_1 = sqrt((w1 + 1.0)/2.0)
	cy_r_1 = sqrt((1.0 - w1)/2.0)

	w2 = var1[i]/2.0
	co_r_2 = sqrt((w2 + 1.0)/2.0)
	cy_r_2 = sqrt((1.0 - w2)/2.0)

	config['qdot']['desired_state_co_r'] = co_r_1
	config['qdot']['desired_state_cy_r'] = cy_r_1
	config['qdot2']['desired_state_co_r'] = co_r_2
	config['qdot2']['desired_state_cy_r'] = cy_r_2

	for j in range(j_initial, j_final):
		print "i,j = {0}, {1}".format(i, j)
		config['qdot']['desired_state_co_phi'] = 0.0					#***********
		config['qdot']['desired_state_cy_phi'] = -var2[j]/2.0
		config['qdot2']['desired_state_co_phi'] = 0.0
		config['qdot2']['desired_state_cy_phi'] = var2[j]/2.0

		config.write()
		subprocess.check_call('python main.py', shell=True)

		t = asciitable.read('data/results.txt', Reader=asciitable.CommentedHeader)
		# extract results from results.txt
		a = t['Parameter']
		fidelity_index = [item for item in range(len(a)) if a[item] == 'Overall Fidelity']
		fidelity[i, j] = t['Value'][fidelity_index]
		# save data to file
		savetxt(data_folder + 'current_iteration.txt', array([i,j]), header='i, j')
		savetxt(data_folder + 'script_results.txt', transpose(fidelity), header='Fidelity')

config['run']['show_plot'] = 'True'
config.write()


X, Y = meshgrid(var1, var2)
contourf(X, Y, transpose(fidelity), linspace(0.0, 1.0, 80))
v = linspace(0.0, 1.0, 11, endpoint=True)
cbar = colorbar(ticks=v, orientation='horizontal')
cbar.ax.set_xlabel('Fidelity')
gca().xaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
gca().yaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
xlabel(r'w2-w1')		#***********
ylabel(r'phi2-phi1 (PI rad)')	#***********
show()


