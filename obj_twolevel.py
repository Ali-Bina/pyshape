#! /home/ajan/anaconda/bin/python
from pulse import *
from numpy import *
from bloch_two_level import *
from matmanip import *
from scipy.integrate import ode

def obj_twolevel(x, params, ampmaskfunc, phasemaskfunc, optimize):

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

    # create pulse
    pulse = Pulse(params)
    pulse.setShape()
    pulse.fft()

    # apply mask to pulse
    pulse.applyMask(ampmaskfunc, phasemaskfunc, params, **kwargs)
    pulse.ifft()
    pulse.rabifreq(params['qdot_params'].qdot_d_yo)
    #print "main",pulse.omega_o

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
    S_0 = state_to_blochvector(params['qdot_params'].qdot_initial_state)
    dmatrix_0 = 1/2.0 * np.array([[1-S_0[2], S_0[0]-1j*S_0[1]] ,[S_0[0] + 1j*S_0[1], 1+S_0[2]]], dtype=complex)

    print "Looky hererere", dmatrix_0

    print "The desired state is ", params['qdot_params'].qdot_desired_state
    desired_S = state_to_blochvector(params['qdot_params'].qdot_desired_state)
    desired_state = blochvector_to_state(desired_S)
    desired_dmatrix = state_to_dmatrix(desired_state)
    desired_dmatrix = 1/2.0 * np.array([[1-desired_S[2], desired_S[0]-1j*desired_S[1]] ,[desired_S[0] + 1j*desired_S[1], 1+desired_S[2]]], dtype=complex)


    # integrate bloch equations
    r = ode(bloch.operator)
    r.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)
    r.set_initial_value(S_0, t_start)

    # determine number of integration steps
    dt = params['rkdp_params'].rkdp_tdep_step
    NSTEPS = int(ceil((t_end - t_start)/dt)) + 1
    dt = (t_end - t_start)/NSTEPS
    

    # integrate the Bloch equations
    if optimize:
        S_end = r.integrate(t_end)
    else:
        t = zeros([NSTEPS])
        t = linspace(t_start, t_end, NSTEPS, endpoint=True)
        S = zeros([NSTEPS, 3])
        i = 0
        while r.successful() and r.t < (t[-1]) :
            i += 1
            r.integrate(t[i])
            S[i,:] = r.y
        


        # calculate the fidelity
        S_end = S[-1,:]

        # plot results
        if params['run_params'].show_plot:
            subplot(2,2,1)
            pulse.plotIntensity()
            #pulse.plotefield()
            subplot(2,2,2)
#            plot(convert(t, ARU_TO_FEMTO), S[:,0], 'r-')
#            plot(convert(t, ARU_TO_FEMTO), S[:,1], 'b-')
            plot(convert(t, ARU_TO_FEMTO), S[:,2], 'g-')
            xlabel('time (fs)')
            ylabel('occupation')
            grid(True)
            axis('tight')
            ylim(-1.0, 1.0)
            subplot(2,2,3)
            #pulse.plotEfieldAmp()
            #pulse.plotefield()
            pulse.plotinstfreq()
            #xlim(1.0, 1.15)
            subplot(2,2,4)
            pulse.plotIntEfield()
            #xlim(1.0, 1.15)
            #ylim(-10.0, 10.0)
            show(block=False)

        # write pulse data to file
        pulse.intensityAC()
        pulse.writePulse()

    # calculate the fidelity
    final_state = blochvector_to_state(S_end)
    final_dmatrix = state_to_dmatrix(final_state)
    print("Final state", final_dmatrix)
    print '-'*10
    final_dmatrix = 1/2.0 * np.array([[1-S_end[2], S_end[0]-1j*S_end[1]] ,[S_end[0] + 1j*S_end[1], 1+S_end[2]]], dtype=complex)
    print("Final state", final_dmatrix)
    fidelity = trace(dot(final_dmatrix, desired_dmatrix)).real
    print '-'*10
    print("desired state", desired_dmatrix)
    if optimize:
        return 1.0 - fidelity
    else:
        data = {"fidelity": fidelity, "t": convert(t, ARU_TO_FEMTO), "S": S, "pulse_data": pulse.data()}
        return data












