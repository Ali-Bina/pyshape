#! /home/ajan/anaconda/bin/python
from pulse import *
from bloch import *
from matmanip import *
from scipy.integrate import ode

def obj_threedot(x, params, ampmaskfunc, phasemaskfunc, optimize):

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

    # apply mask to pulse
    pulse.applyMask(ampmaskfunc, phasemaskfunc, params, **kwargs)
    pulse.ifft()

    # define integration parameters
    rtol = params['rkdp_params'].rkdp_rtol
    atol = params['rkdp_params'].rkdp_atol
    t_start = params['run_params'].run_t_start
    t_end = params['run_params'].run_t_end
    dt = params['rkdp_params'].rkdp_tdep_step
    NSTEPS = int(ceil((t_end - t_start)/dt))
    dt = (t_end - t_start)/NSTEPS
    if optimize:
        nsteps = params['rkdp_params'].rkdp_nsteps_opt
    else:
        nsteps = params['rkdp_params'].rkdp_nsteps_tdep

    #  ----  FIRST DOT  ----  #
    # define initial conditions
    initial_dmatrix = state_to_dmatrix(params['qdot_params'].qdot_initial_state)
    rho0_1 = dmatrix_to_array(initial_dmatrix)

    # define desired final state
    qdot1_desired_dmatrix = state_to_dmatrix(params['qdot_params'].qdot_desired_state)
    rhod_1 = dmatrix_to_array(qdot1_desired_dmatrix)

    # integrate bloch equations
    bloch_1 = BlochEqns(params, pulse.efield, pulse.t, DOT1)                    # define bloch equations with shaped electric field
    r_1 = ode(bloch_1.operator)
    r_1.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)
    r_1.set_initial_value(rho0_1, t_start)

    #  ----  SECOND DOT  ----  #
    # define initial conditions
    qdot2_initial_dmatrix = state_to_dmatrix(params['qdot2_params'].qdot_initial_state)
    rho0_2 = dmatrix_to_array(qdot2_initial_dmatrix)

    # define desired final state
    qdot2_desired_dmatrix = state_to_dmatrix(params['qdot2_params'].qdot_desired_state)
    rhod_2 = dmatrix_to_array(qdot2_desired_dmatrix)

    # integrate bloch equations
    bloch_2 = BlochEqns(params, pulse.efield, pulse.t, DOT2)                    # define bloch equations with shaped electric field
    r_2 = ode(bloch_2.operator)
    r_2.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)
    r_2.set_initial_value(rho0_2, t_start)

    #  ----  THIRD DOT  ----  #
    # define initial conditions
    qdot3_initial_dmatrix = state_to_dmatrix(params['qdot3_params'].qdot_initial_state)
    rho0_3 = dmatrix_to_array(qdot3_initial_dmatrix)

    # define desired final state
    qdot3_desired_dmatrix = state_to_dmatrix(params['qdot3_params'].qdot_desired_state)
    rhod_3 = dmatrix_to_array(qdot3_desired_dmatrix)

    # integrate bloch equations
    bloch_3 = BlochEqns(params, pulse.efield, pulse.t, DOT3)                    # define bloch equations with shaped electric field
    r_3 = ode(bloch_3.operator)
    r_3.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)
    r_3.set_initial_value(rho0_3, t_start)

    if optimize:
        qdot1_rho_end = r_1.integrate(t_end)
        qdot2_rho_end = r_2.integrate(t_end)
        qdot3_rho_end = r_3.integrate(t_end)

    else:
        t = zeros([NSTEPS])
        rho_1 = zeros([NSTEPS, 16])
        # integrate the Bloch equations
        i = 0
        while r_1.successful() and r_1.t < t_end:
            r_1.integrate(r_1.t + dt)
            rho_1[i,:] = r_1.y
            t[i] = r_1.t
            i =  i + 1
        qdot1_rho_end = rho_1[-1,:]

        rho_2 = zeros([NSTEPS, 16])
        # integrate the Bloch equations
        i = 0
        while r_2.successful() and r_2.t < t_end:
            r_2.integrate(r_2.t + dt)
            rho_2[i,:] = r_2.y
            i =  i + 1
        qdot2_rho_end = rho_2[-1,:]

        rho_3 = zeros([NSTEPS, 16])
        # integrate the Bloch equations
        i = 0
        while r_3.successful() and r_3.t < t_end:
            r_3.integrate(r_3.t + dt)
            rho_3[i,:] = r_3.y
            i =  i + 1
        qdot3_rho_end = rho_3[-1,:]

        # plot the results
        if params['run_params'].show_plot:
            subplot(2,2,1)
            pulse.plotIntensity()
            subplot(2,2,2)
            plot(convert(t, ARU_TO_FEMTO), rho_1[:,0], 'b-')
            plot(convert(t, ARU_TO_FEMTO), rho_1[:,1], 'g-')
            plot(convert(t, ARU_TO_FEMTO), rho_1[:,2], 'r-')
            plot(convert(t, ARU_TO_FEMTO), rho_1[:,3], 'k-')
            plot(convert(t, ARU_TO_FEMTO), rho_2[:,0], 'b--')
            plot(convert(t, ARU_TO_FEMTO), rho_2[:,1], 'g--')
            plot(convert(t, ARU_TO_FEMTO), rho_2[:,2], 'r--')
            plot(convert(t, ARU_TO_FEMTO), rho_2[:,3], 'k--')
            plot(convert(t, ARU_TO_FEMTO), rho_3[:,0], 'b:')
            plot(convert(t, ARU_TO_FEMTO), rho_3[:,1], 'g:')
            plot(convert(t, ARU_TO_FEMTO), rho_3[:,2], 'r:')
            plot(convert(t, ARU_TO_FEMTO), rho_3[:,3], 'k:')
            xlabel('time (fs)')
            ylabel('occupation')
            grid(True)
            axis('tight')
            subplot(2,2,3)
            pulse.plotEfieldAmp()
            xlim(1.0, 1.15)
            subplot(2,2,4)
            pulse.plotEfieldPhase()
            xlim(1.0, 1.15)
            show()

        # write pulse data to file
        pulse.intensityAC()
        pulse.writePulse()

    # calculate the fidelity for dot 1
    qdot1_final_dmatrix = array_to_dmatrix(qdot1_rho_end)
    fidelity_1 = trace(dot(qdot1_final_dmatrix, qdot1_desired_dmatrix)).real

    # calculate the fidelity for dot 2
    qdot2_final_dmatrix = array_to_dmatrix(qdot2_rho_end)
    fidelity_2 = trace(dot(qdot2_final_dmatrix, qdot2_desired_dmatrix)).real

    # calculate the fidelity for dot 3
    qdot3_final_dmatrix = array_to_dmatrix(qdot3_rho_end)
    fidelity_3 = trace(dot(qdot3_final_dmatrix, qdot3_desired_dmatrix)).real

    # calculate the overall fidelity for the gate
    dp_state_physical = kron(kron(qdot1_final_dmatrix, qdot2_final_dmatrix), qdot3_final_dmatrix)
    dp_state_ideal = kron(kron(qdot1_desired_dmatrix, qdot2_desired_dmatrix), qdot3_desired_dmatrix)
    fidelity = trace(dot(dp_state_physical, dp_state_desired)).real

    if optimize:
        # prevent pulse area from going negative (not sure why this happens)
        if params['pulse_params'].pulse_EO < params['mask_params'].x_bounds[3][0] or params['pulse_params'].pulse_EO > params['mask_params'].x_bounds[3][1]:
            fidelity = 0.0
        return 1.0 - fidelity
    else:
        data = {"fidelity": fidelity, "fidelity dot1": fidelity_1, "fidelity dot2": fidelity_2, "fidelity dot3": fidelity_3, "t": convert(t, ARU_TO_FEMTO), "rho_1": rho_1, "rho_2": rho_2, "rho_3": rho_3, "pulse_data": pulse.data(), "rho0_1": rho0_1, "rhod_1": rhod_1, "rho0_2": rho0_2, "rhod_2": rhod_2, "rho0_3": rho0_3, "rhod_3": rhod_3}
        return data









