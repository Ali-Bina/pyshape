from pulse import *
from bloch_two_level_dm_sx import *
from matmanip import *
from scipy.integrate import ode
import matplotlib.pyplot as plt

def obj_twodot_two_level_dm(x, params, ampmaskfunc, phasemaskfunc, optimize):

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
    hstart = params['rkdp_params'].rkdp_hstart
    hmax = params['rkdp_params'].rkdp_hmax

    t_start = params['run_params'].run_t_start
    t_end = params['run_params'].run_t_end
    t_step = params['rkdp_params'].rkdp_tdep_step

    array_len = int((convert(t_end, ARU_TO_FEMTO) - convert(t_start, ARU_TO_FEMTO))/convert(t_step, ARU_TO_FEMTO))

    if optimize:
        nsteps = params['rkdp_params'].rkdp_nsteps_opt
    else:
        nsteps = params['rkdp_params'].rkdp_nsteps_tdep

    #  ----  FIRST DOT  ----  #
    # define initial conditions
    initial_dmatrix = state_to_dmatrix(params['qdot_params'].qdot_initial_state)
    rho0_1 = dmatrix_to_array_two_level_dm(initial_dmatrix)

    # define desired final state
    qdot1_desired_dmatrix = state_to_dmatrix(params['qdot_params'].qdot_desired_state)
    rhod_1 = dmatrix_to_array_two_level_dm(qdot1_desired_dmatrix)

    # integrate bloch equations
    bloch_1 = BlochEqns(params, pulse.efield, pulse.t, DOT1)                    # define bloch equations with shaped electric field
    r_1 = ode(bloch_1.operator)
    r_1.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, first_step=hstart)
    r_1.set_initial_value(rho0_1, t_start)

    #  ----  SECOND DOT  ----  #
    # define initial conditions
    qdot2_initial_dmatrix = state_to_dmatrix(params['qdot2_params'].qdot_initial_state)
    rho0_2 = dmatrix_to_array_two_level_dm(qdot2_initial_dmatrix)

    # define desired final state
    qdot2_desired_dmatrix = state_to_dmatrix(params['qdot2_params'].qdot_desired_state)
    rhod_2 = dmatrix_to_array_two_level_dm(qdot2_desired_dmatrix)

    # integrate bloch equations
    bloch_2 = BlochEqns(params, pulse.efield, pulse.t, DOT2)                    # define bloch equations with shaped electric field
    r_2 = ode(bloch_2.operator)
    r_2.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, first_step=hstart)
    r_2.set_initial_value(rho0_2, t_start)

    if optimize:
        qdot1_rho_end = r_1.integrate(t_end - t_step)
        qdot2_rho_end = r_2.integrate(t_end - t_step)

    else:
        t = numpy.linspace(t_start, t_end, array_len, endpoint=False)
        rho_1 = numpy.zeros([array_len, 3])
        rho_2 = numpy.zeros([array_len, 3])

        # integrate the Bloch equations
        i = 0
        while r_1.successful() and r_1.t < t[-1]:
            i += 1
            r_1.integrate(t[i])
            rho_1[i,:] = r_1.y
        qdot1_rho_end = rho_1[-1,:]

        # integrate the Bloch equations
        i = 0
        while r_2.successful() and r_2.t < t[-1]:
            i += 1
            r_2.integrate(t[i])
            rho_2[i,:] = r_2.y
        qdot2_rho_end = rho_2[-1,:]

        # plot the results
        if params['run_params'].show_plot:
            plt.subplot(2,2,1)
            pulse.plotIntensity()
            plt.subplot(2,2,2)
            plt.plot(convert(t, ARU_TO_FEMTO), 1.0 - rho_1[:,0], 'b-')
            plt.plot(convert(t, ARU_TO_FEMTO), rho_1[:,0], 'g-')
            plt.plot(convert(t, ARU_TO_FEMTO), rho_1[:,1], 'r-')
            plt.plot(convert(t, ARU_TO_FEMTO), rho_1[:,2], 'k-')
            plt.plot(convert(t, ARU_TO_FEMTO), 1.0 - rho_2[:,0], 'b--')
            plt.plot(convert(t, ARU_TO_FEMTO), rho_2[:,0], 'g--')
            plt.plot(convert(t, ARU_TO_FEMTO), rho_2[:,1], 'r--')
            plt.plot(convert(t, ARU_TO_FEMTO), rho_2[:,2], 'k--')
            plt.xlabel('time (fs)')
            plt.ylabel('occupation')
            plt.grid(True)
            plt.axis('tight')
            plt.subplot(2,2,3)
            pulse.plotEfieldPhase()
            plt.subplot(2,2,4)
            pulse.plotrealimagefield()
            plt.show()

        # write pulse data to file
        pulse.intensityAC()
        pulse.writePulse()

    # calculate the fidelity for dot 1
    qdot1_final_dmatrix = array_to_dmatrix_two_level_dm(qdot1_rho_end)
    fidelity_1 = numpy.trace(numpy.dot(qdot1_final_dmatrix, qdot1_desired_dmatrix)).real

    # calculate the fidelity for dot 2
    qdot2_final_dmatrix = array_to_dmatrix_two_level_dm(qdot2_rho_end)
    fidelity_2 = numpy.trace(numpy.dot(qdot2_final_dmatrix, qdot2_desired_dmatrix)).real

    # calculate the overall fidelity for the gate

    dp_state_physical = numpy.kron(qdot1_final_dmatrix, qdot2_final_dmatrix)
    dp_state_desired = numpy.kron(qdot1_desired_dmatrix, qdot2_desired_dmatrix)
    fidelity = numpy.trace(numpy.dot(dp_state_physical, dp_state_desired)).real

    if optimize:
        return 1.0 - fidelity
    else:
        data = {"fidelity": fidelity, "fidelity dot1": fidelity_1, "fidelity dot2": fidelity_2, "t": convert(t, ARU_TO_FEMTO), "rho_1": rho_1, "rho_2": rho_2, "pulse_data": pulse.data(), "rho0_1": rho0_1, "rhod_1": rhod_1, "rho0_2": rho0_2, "rhod_2": rhod_2}
        return data









