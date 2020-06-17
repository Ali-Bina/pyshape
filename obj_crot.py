#! /home/ajan/anaconda/bin/python
from pulse import *
from bloch import *
from scipy.integrate import ode

def obj_crot(x, params, ampmaskfunc, phasemaskfunc, optimize):

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

    # define bloch equations with shaped electric field
    bloch = BlochEqns(params, pulse.efield, pulse.t, DOT1)

    # define integration parameters
    rtol = params['rkdp_params'].rkdp_rtol
    atol = params['rkdp_params'].rkdp_atol
    nsteps = params['rkdp_params'].rkdp_nsteps
    t_start = params['run_params'].run_t_start
    t_end = params['run_params'].run_t_end

    initial_index = array([0, 1, 2, 3])
    final_index = array([0, 1, 3, 2])
    NTRIALS = len(initial_index)
    fidelity_results = zeros(NTRIALS)

    # integrate bloch equations
    r = ode(bloch.operator)
    r.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)


    dt = params['rkdp_params'].rkdp_tdep_step
    NSTEPS = int((t_end - t_start)/dt)
    dt = (t_end - t_start)/NSTEPS

    # integrate Bloch equations for the four initial conditions
    for i in range(NTRIALS):
        # define initial conditions
        rho0 = zeros(16)
        rho0[initial_index[i]] = 1.0

        # integrate Bloch equations
        r = ode(bloch.operator)
        r.set_integrator('dop853', method='adams', rtol=rtol, atol=atol)
        r.set_initial_value(rho0, t_start)

        if optimize:

            rho = r.integrate(t_end)
            fidelity_results[i] = rho[final_index[i]]

        else:
            t = zeros([NSTEPS])
            rho = zeros([NTRIALS, NSTEPS, 16])

            # integrate the Bloch equations
            j = 0
            while r.successful() and r.t < (t_end - dt):
                r.integrate(r.t + dt)
                rho[i, j,:] = r.y
                t[j] = r.t
                j =  j + 1

            fidelity_results[i] = rho[i, NSTEPS - 1, final_index[i]]

    # average the results for the four initial states
    fidelity = average(fidelity_results)
    print fidelity_results
    print fidelity

    if optimize is False:
        # plot the results
        if params['run_params'].show_plot:
            subplot(2,2,1)
            pulse.plotIntensity()
            subplot(2,2,2)
            plot(convert(t, ARU_TO_FEMTO), rho[2,:,0])
            plot(convert(t, ARU_TO_FEMTO), rho[2,:,1])
            plot(convert(t, ARU_TO_FEMTO), rho[2,:,2])
            plot(convert(t, ARU_TO_FEMTO), rho[2,:,3])
            xlabel('time (fs)')
            ylabel('occupation')
            grid(True)
            axis('tight')
            subplot(2,2,3)
            pulse.plotEfieldAmp()
            subplot(2,2,4)
            pulse.plotEfieldPhase()
            show()


    if optimize:
        return 1.0 - fidelity
    else:
        data = {"fidelity": fidelity, "t": convert(t, ARU_TO_FEMTO), "rho_1": rho, "pulse_data": pulse.data()}
        return data












