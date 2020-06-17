import sys
import os
import sobol_lib
from scipy import optimize
import random

import laphononeid

from convert import *
from read_config import *
from phasemasks import *
from obj_twodot_two_level_dm import *
from obj_twolevel import *
from obj_threelevel import *
from obj_onedot import *
from obj_twodot import *
from obj_threedot import *
from obj_crot import *

def main(args):
    # read in config parameters and validate against spec file
    config_file = 'config.ini'
    configspec_file = 'configspec.ini'
    params = read_config(config_file, configspec_file)

    # generate the LA phonon eid kernel if required
    if(params['run_params'].run_phonon_eid):
        laphononeid.eidkernel(params)

    # choose objective function from objective function dictionaries
    obj_gates = {'twodottwoleveldm': obj_twodot_two_level_dm,  'twolevel': obj_twolevel, 'threelevel': obj_threelevel, 'onedot': obj_onedot, 'twodot': obj_twodot, 'threedot': obj_threedot, 'crot': obj_crot}
    obj_func = obj_gates[params['run_params'].gate]

    # choose amplitude and phase mask function from the mask function dictionary
    ampmasks = {'ampmask_chen': chenAmpMask, 'mask_slmcc': slmccAmpMask, 'ampmask_none': noAmpMask}
    phasemasks = {'phasemask_cos': cosinePhaseMask, 'phasemask_poly': polyPhaseMask, 'mask_slmcc': slmccPhaseMask, 'phasemask_none': noPhaseMask}
    key = params['run_params'].run_ampmask
    if key in ampmasks:
        ampmaskfunc = ampmasks[key]
    else:
        print "amp mask not valid"
    key = params['run_params'].run_phasemask
    if key in phasemasks:
        phasemaskfunc = phasemasks[key]
    else:
        print "phase mask not valid"

    # read in configuration file and convert parameters to atomic Rhydberg units
    x = read_params(config_file)
    x = convert_aru(x, config_file, TO_ARU)

    # run optimization if required and get time dependence
    run_optimize = params['run_params'].optimize
    if run_optimize==True:
        x = nloptimize(obj_func, ampmaskfunc, phasemaskfunc, params, config_file)
        fidelity = timedep(obj_func, ampmaskfunc, phasemaskfunc, params, config_file)
    else:
        fidelity = timedep(obj_func, ampmaskfunc, phasemaskfunc, params, config_file)

    # clean up directory
    os.system(" rm *.pyc")



def nloptimize(obj_func, ampmask, phasemask, params, config_file):
    ''' optimize obj_func gate using ampmask and phasemask '''

    print "optimizing {0} gate fidelity using {1} and {2}".format(params['run_params'].gate, params['run_params'].run_ampmask, params['run_params'].run_phasemask)

    NOPTS = params['mask_params'].NOPTS                            # number of parameters to optimize
    NITER = params['run_params'].NITER                            # number of initial vectors

    x_start = zeros(NOPTS)                                                  #  array to hold initial values for parameters
    x_bounds = params['mask_params'].x_bounds                               # upper and lower bounds for parameters

    # varaiables to store optimal fidelity and parameter values
    optimal_fidelity = 0.0
    optimal_x = zeros(NOPTS)
    cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - x_bounds[0][0]},
            {'type': 'ineq', 'fun': lambda x:  x_bounds[0][1] - x[0]},
            {'type': 'ineq', 'fun': lambda x:  x[1] - x_bounds[1][0]},
            {'type': 'ineq', 'fun': lambda x:  x_bounds[1][1] - x[1]},
            {'type': 'ineq', 'fun': lambda x:  x[2] - x_bounds[2][0]},
            {'type': 'ineq', 'fun': lambda x:  x_bounds[2][1] - x[2]},
            {'type': 'ineq', 'fun': lambda x:  x[3] - x_bounds[3][0]},
            {'type': 'ineq', 'fun': lambda x:  x_bounds[3][1] - x[3]})

    # optimize for NITER different initial values
    seed = params['run_params'].sobol_seed   # initialize sobol sequence seed
    run_optimize = True
    tol = 1.0e-3
    for i in range(NITER):
        # choose random initial values for parameters from allowed phase space
        [ ran_var, seed_out ] = sobol_lib.i4_sobol(NOPTS, seed)  # use sobol sequence to generate NOPT random variables
        seed = seed_out     # update seed for sobol sequence
        for j in range(NOPTS):
            #ran_var[j] = random.random()
            x_start[j] = (x_bounds[j][1] - x_bounds[j][0])*ran_var[j] + x_bounds[j][0]

        # tuple of arguments for the minimization routine - first set of arguments is for the obj_func
        args = ((params, ampmask, phasemask, run_optimize))

        # optimize for single instance of x_start, use options={'disp': True} to display more info
        output = optimize.minimize(obj_func, x_start, args=args, method='SLSQP', constraints=cons, tol=tol)
        x = output.x
        fidelity = 1 - output.fun

        if fidelity > optimal_fidelity:
            optimal_x = x
            optimal_fidelity = fidelity
        print "Iteration {0}, Current Result: {1}, Best Result: {2}".format(i, fidelity, optimal_fidelity)

    # convert optimal result to normal units and write to config file
    converted_optimal_x = convert_aru(optimal_x, config_file, FROM_ARU)
    # update config file with optimal parameters
    write_params(converted_optimal_x, config_file)

    # print results to screen
    print "Optimal results: "
    print "Fidelity =", optimal_fidelity
    for i in range(NOPTS):
        print "{0} = {1} {2}".format(params['mask_params'].param_list[i], optimal_x[i], params['mask_params'].x_units[i])
    print "\n"

def timedep(obj_func, ampmaskfunc, phasemaskfunc, params, config_file):
    ''' read in parameters from config file and carry out integration using small time steps '''

    print "time dependence for {0} gate using {1} and {2}".format(params['run_params'].gate, params['run_params'].run_ampmask, params['run_params'].run_phasemask)

    # read in parameters from file and convert to atomic Rhydberg units
    x = read_params(config_file)
    x = convert_aru(x, config_file, TO_ARU)

    run_optimize = False

    # run objective function for single instance and convert results back to normal units
    data = obj_func(x, params, ampmaskfunc, phasemaskfunc, run_optimize)
    x_converted = convert_aru(x, config_file, FROM_ARU)

    # print results to screen
    print "\nFidelity =", data["fidelity"]
    NOPTS = params['mask_params'].NOPTS
    for i in range(NOPTS):
        print "{0} = {1} {2}".format(params['mask_params'].param_list[i], x_converted[i], params['mask_params'].x_units[i])

    # save results to the data folder
    os.system("cp config.ini data/")
    savetxt("data/dynamics/tarray.txt", data["t"])
    if (params['run_params'].gate == 'onedot' or params['run_params'].gate == 'crot'):
        savetxt("data/dynamics/rho_1.txt", data["rho_1"])
        savetxt("data/dynamics/rho0_1.txt", data["rho0_1"])
        savetxt("data/dynamics/rhod_1.txt", data["rhod_1"])
    elif params['run_params'].gate == 'twodot':
        savetxt("data/dynamics/rho_1.txt", data["rho_1"])
        savetxt("data/dynamics/rho_2.txt", data["rho_2"])
        savetxt("data/dynamics/rho0_1.txt", data["rho0_1"])
        savetxt("data/dynamics/rhod_1.txt", data["rhod_1"])
        savetxt("data/dynamics/rho0_2.txt", data["rho0_2"])
        savetxt("data/dynamics/rhod_2.txt", data["rhod_2"])
    elif params['run_params'].gate == 'twodottwoleveldm':
        savetxt("data/dynamics/rho_1.txt", data["rho_1"])
        savetxt("data/dynamics/rho_2.txt", data["rho_2"])
        savetxt("data/dynamics/rho0_1.txt", data["rho0_1"])
        savetxt("data/dynamics/rhod_1.txt", data["rhod_1"])
        savetxt("data/dynamics/rho0_2.txt", data["rho0_2"])
        savetxt("data/dynamics/rhod_2.txt", data["rhod_2"])
    elif params['run_params'].gate == 'threedot':
        savetxt("data/dynamics/rho_1.txt", data["rho_1"])
        savetxt("data/dynamics/rho_2.txt", data["rho_2"])
        savetxt("data/dynamics/rho_3.txt", data["rho_3"])
        savetxt("data/dynamics/rho0_1.txt", data["rho0_1"])
        savetxt("data/dynamics/rhod_1.txt", data["rhod_1"])
        savetxt("data/dynamics/rho0_2.txt", data["rho0_2"])
        savetxt("data/dynamics/rhod_2.txt", data["rhod_2"])
        savetxt("data/dynamics/rho0_3.txt", data["rho0_3"])
        savetxt("data/dynamics/rhod_3.txt", data["rhod_3"])
    elif params['run_params'].gate == 'twolevel':
        savetxt("data/dynamics/S.txt", data["S"])

    with open("data/results.txt", "w") as resultsfile:
        resultsfile.write("#Parameter, Value, Units\n")
        for i in range(NOPTS):
            resultsfile.write("%s, %.12f, %s\n" % (params['mask_params'].param_list[i], x_converted[i], params['mask_params'].x_units[i]))
        if (params['run_params'].gate == 'onedot' or params['run_params'].gate == 'crot' or params['run_params'].gate == 'twolevel' or params['run_params'].gate == 'threelevel'):
            resultsfile.write("Fidelity, %.12f, arb. units\n" % (data["fidelity"]))
        elif params['run_params'].gate == 'twodot':
            resultsfile.write("Fidelity DOT 1, %.12f, arb. units\n" % (data["fidelity dot1"]))
            resultsfile.write("Fidelity DOT 2, %.12f, arb. units\n" % (data["fidelity dot2"]))
            resultsfile.write("Overall Fidelity, %.12f, arb. units\n" % (data["fidelity"]))
        elif params['run_params'].gate == 'twodottwoleveldm':
            resultsfile.write("Fidelity DOT 1, %.12f, arb. units\n" % (data["fidelity dot1"]))
            resultsfile.write("Fidelity DOT 2, %.12f, arb. units\n" % (data["fidelity dot2"]))
            resultsfile.write("Overall Fidelity, %.12f, arb. units\n" % (data["fidelity"]))
        elif params['run_params'].gate == 'threedot':
            resultsfile.write("Fidelity DOT 1, %.12f, arb. units\n" % (data["fidelity dot1"]))
            resultsfile.write("Fidelity DOT 2, %.12f, arb. units\n" % (data["fidelity dot2"]))
            resultsfile.write("Fidelity DOT 3, %.12f, arb. units\n" % (data["fidelity dot3"]))
            resultsfile.write("Overall Fidelity, %.12f, arb. units\n" % (data["fidelity"]))

    return data["fidelity"]


if __name__ == "__main__":
    sys.exit(main(sys.argv))



