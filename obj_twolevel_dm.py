#! /home/ajan/anaconda/bin/python
from __future__ import division
from pulse import *
from bloch_two_level_dm import *
from matmanip import *
from scipy.integrate import ode
#import math
import numpy
import matplotlib.pyplot as plt


def obj_twolevel_dm(x, params, ampmaskfunc, phasemaskfunc, optimize):

    # unpack x and assign to variables
    param_list = params['mask_params'].param_list
    NOPTS = params['mask_params'].NOPTS
    kwargs = {}
    
    # create keyword arguments
    for i in range(NOPTS):
        kwargs[param_list[i]] = x[i]

    # assign pulse area if it is being optimized
    if 'pulse_area' in kwargs:
        params['pulse_params'].pulse_EO = kwargs['pulse_area']

    pulse = Pulse(params)
    pulse.setShape()
    pulse.fft()
    delta=abs(pulse.omega_o-params['qdot_params'].qdot_omega_yo)

    # apply mask to pulse
    pulse.applyMask(ampmaskfunc, phasemaskfunc, params, **kwargs)
    pulse.ifft()
    pulse.rabifreq(params['qdot_params'].qdot_d_yo)
#    pulse.gen_rabifreq(params['qdot_params'].qdot_d_yo,delta)
    rabi=pulse.rabifreq(params['qdot_params'].qdot_d_yo)
    pulse.plotefield()

    pulse.writePulse()
    # define bloch equations with shaped electric field
    bloch = BlochEqns(params, pulse.efield, pulse.t, DOT1)

    # define integration parameters
    rtol = params['rkdp_params'].rkdp_rtol
    atol = params['rkdp_params'].rkdp_atol
    t_start = params['run_params'].run_t_start
    t_end = params['run_params'].run_t_end
    if optimize:
        nsteps = params['rkdp_params'].rkdp_nsteps_opt
    else:
        nsteps = params['rkdp_params'].rkdp_nsteps_tdep

    # define initial conditions
    initial_dmatrix = state_to_dmatrix(params['qdot_params'].qdot_initial_state)
    rho0_1 = dmatrix_to_array(initial_dmatrix)
    #print initial_dmatrix

    # define desired final state
    qdot1_desired_dmatrix = state_to_dmatrix(params['qdot_params'].qdot_desired_state)
    rhod_1 = dmatrix_to_array(qdot1_desired_dmatrix)

    # integrate bloch equations
    r = ode(bloch.operator)
    r.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)
    r.set_initial_value(rho0_1, t_start)

    dt = params['rkdp_params'].rkdp_tdep_step
    NSTEPS = int(numpy.ceil((t_end - t_start)/dt))+1
    dt = (t_end - t_start)/NSTEPS

    if optimize:
        rho_end = r.integrate(t_end)

    else:
        t = numpy.zeros([NSTEPS])
        rho = numpy.zeros([NSTEPS, 16])

        # integrate the Bloch equations
        i = 0
        while r.successful() and r.t < t_end:
	    
            r.integrate(r.t + dt)
            rho[i-1,:] = r.y
            t[i-1] = r.t
            i =  i + 1
        rho_end = rho[-1,:]

        # plot the results
        if params['run_params'].show_plot:
            plt.subplot(2,2,1)
            pulse.plotIntensity()
            #pulse.plotinstfreq()
            plt.subplot(2,2,2)
            plt.plot(convert(t, ARU_TO_FEMTO), rho[:,0], 'r')
            plt.plot(convert(t, ARU_TO_FEMTO), rho[:,1], 'b')
            plt.plot(convert(t, ARU_TO_FEMTO), rho[:,2], 'g')
            plt.plot(convert(t, ARU_TO_FEMTO), rho[:,3], 'k')
            plt.xlabel('time (fs)')
            plt.ylabel('occupation')
            plt.grid(True)
            plt.axis('tight')
            plt.subplot(2,2,3)
            #pulse.plotEfieldAmp()
            #pulse.plotrabifreq()
            #plt.xlim(1.0, 1.15)
            plt.subplot(2,2,4)
            pulse.plotEfieldPhase()
            #plt.xlim(1.0, 1.15)
            #plt.ylim(-10.0, 10.0)
            plt.show(block=False)

        # write pulse data to file
        pulse.intensityAC()
        pulse.writePulse()

    # calculate the fidelity
    final_dmatrix = array_to_dmatrix(rho_end)
    fidelity = numpy.trace(numpy.dot(final_dmatrix,qdot1_desired_dmatrix)).real

    if optimize:
        return 1.0 - fidelity
    else:
        data = {"fidelity": fidelity, "t": convert(t, ARU_TO_FEMTO), "rho_1": rho, "pulse_data": pulse.data(), "rho0_1": rho0_1, "rhod_1": rhod_1}
        return data













