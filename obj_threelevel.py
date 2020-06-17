#! /home/ajan/anaconda/bin/python
from pulse import *
from bloch_three_level import *
from matmanip import *
from scipy.integrate import ode

def obj_threelevel(x, params, ampmaskfunc, phasemaskfunc, optimize):

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
    pulse.rabifreq(params['qdot_params'].qdot_d_yo)

    # define bloch equations with shaped electric field
    bloch = BlochEqns(params, pulse.efield, pulse.t, DOT1)

    # define integration parameters
    rtol = params['rkdp_params'].rkdp_rtol
    atol = params['rkdp_params'].rkdp_atol
    nsteps = params['rkdp_params'].rkdp_nsteps
    t_start = params['run_params'].run_t_start
    t_end = params['run_params'].run_t_end

    # define initial conditions
    #rho0_1 = initial_rho(params['qdot_params'].qdot_initial_state, params)
    rho0_1 = zeros(5)
    rho0_1[0] = 1.0


    # integrate bloch equations
    r = ode(bloch.operator)
    r.set_integrator('dop853', method='adams', rtol=rtol, atol=atol, nsteps=nsteps)
    r.set_initial_value(rho0_1, t_start)

    dt = params['rkdp_params'].rkdp_tdep_step
    NSTEPS = int(ceil((t_end - t_start)/dt))
    dt = (t_end - t_start)/NSTEPS

    if optimize:
        S = r.integrate(t_end)
        fidelity = (S[2] + 1.0)/2.0
    else:
        t = zeros([NSTEPS])
        rho = zeros([NSTEPS, 5])
        i = 0

        # integrate the Bloch equations
        while r.successful() and r.t < t_end:
            r.integrate(r.t + dt)
            rho[i,:] = r.y
            t[i] = r.t
            i =  i + 1

        # calculate the fidelity
        fstate = rho[-1,:]
        fidelity = fstate[1] + fstate[2]

        # plot the results
        if params['run_params'].show_plot:
            subplot(2,2,1)
            pulse.plotIntensity()
            #pulse.plotinstfreq()
            subplot(2,2,2)
            plot(convert(t, ARU_TO_FEMTO), rho[:,0], 'r')
            plot(convert(t, ARU_TO_FEMTO), rho[:,1], 'b')
            plot(convert(t, ARU_TO_FEMTO), rho[:,2], 'g')
            xlabel('time (fs)')
            ylabel('occupation')
            grid(True)
            axis('tight')
            subplot(2,2,3)
            pulse.plotEfieldAmp()
            #pulse.plotrabifreq()
            #xlim(1.0, 1.15)
            subplot(2,2,4)
            pulse.plotEfieldPhase()
            #xlim(1.0, 1.15)
            #ylim(-10.0, 10.0)
            show()

        # write pulse data to file
        pulse.intensityAC()
        pulse.writePulse()

    if optimize:
        return 1.0 - fidelity
    else:
        data = {"fidelity": fidelity, "t": convert(t, ARU_TO_FEMTO), "rho_1": rho, "pulse_data": pulse.data()}
        return data










