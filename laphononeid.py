#! /home/ajan/anaconda/bin/python
# Filename: pulse.py

import matplotlib.pyplot as plt
import numpy 
from scipy import integrate
from constants import *
from convert import *


def eidkernel(params):
    if (params['laphnoneid_params'].read_kernel):
        print '\nReading LA Phonon EID Kernel from file...'
        Omega_eid = numpy.loadtxt('data/eid/Omega_eid.txt', delimiter=',', unpack=True, skiprows=1)
        K_real_eid = numpy.loadtxt('data/eid/K_real.txt', delimiter=',', unpack=True, skiprows=1)
        K_imag_eid = numpy.loadtxt('data/eid/K_imag.txt', delimiter=',', unpack=True, skiprows=1)
        Omega_eid = convert(Omega_eid, EV_TO_ARU)/1000.0                   # convert from meV to ARU
        K_real_eid = convert(K_real_eid, INV_PICO_TO_ARU)                  # convert from ps-1 to ARU
        K_imag_eid = convert(K_imag_eid, EV_TO_ARU)/1000.0                 # convert from meV to ARU

        '''
        subplot(2,1,1)
        plot(1000.0*convert(Omega_eid, ARU_TO_EV), convert(K_real_eid, ARU_TO_INV_PICO), '-', label='K real')
        xlabel('Omega (meV)')
        ylabel('Re[$K(\hbar\Omega)$] (ps$^-1$)')
        subplot(2,1,2)
        plot(1000.0*convert(Omega_eid, ARU_TO_EV), 1000.0*convert(K_imag_eid, ARU_TO_EV), '-', label='K imag')
        xlabel('Rabi Energy $\hbar\Omega$ (meV)')
        ylabel('Im[$\hbar K(\hbar\Omega)$] (meV)')
        show()
        '''
        params['laphnoneid_params'].Omega_eid = Omega_eid
        params['laphnoneid_params'].K_real_eid = K_real_eid
        params['laphnoneid_params'].K_imag_eid = K_imag_eid
    else:
        print '\nCalculating LA Phonon EID Kernel...'

        # eid parameters from config file
        omega_c_eid = params['laphnoneid_params'].omega_c_eid
        alpha_eid = params['laphnoneid_params'].alpha_eid
        Omega_start_eid = params['laphnoneid_params'].Omega_start_eid
        Omega_end_eid = params['laphnoneid_params'].Omega_end_eid
        Omega_step_eid = params['laphnoneid_params'].Omega_step_eid
        T_eid = params['laphnoneid_params'].T_eid

        t_start_eid = convert(0.0, FEMTO_TO_ARU)
        t_end_eid = convert(5000.0, FEMTO_TO_ARU)
        t_step_eid = convert(0.1, FEMTO_TO_ARU)
        t_array_len_eid = int((t_end_eid - t_start_eid)/t_step_eid)
        t_step_eid = (t_end_eid - t_start_eid)/t_array_len_eid
        t_eid = numpy.linspace(t_start_eid, t_end_eid, t_array_len_eid)

        def Kt_integrand(omega, t, T):
            return alpha_eid*omega**3*numpy.exp(-(omega/omega_c_eid)**2)*numpy.cos(omega*t)/numpy.tanh(H_BAR_ARU*omega/(2.0*KB_ARU*T_eid))
            #return alpha_eid*omega**3*numpy.exp(-(omega/omega_c_eid)**2)*(1-(2/3.0)*pow(omega/omega_c_eid, 2.0)+(1/5.0)*pow(omega/omega_c_eid, 4.0))*numpy.cos(omega*t)/numpy.tanh(H_BAR_ARU*omega/(2.0*KB_ARU*T_eid))
            #return alpha_eid*pow(omega, 3.0)*(numpy.exp(-pow(omega/omega_c_eid, 2.0))+0.5*numpy.exp(-pow(omega/(2.5*omega_c_eid), 2.0))-numpy.exp(-pow(omega/(omega_c_eid),2.0)*(1.313)))*numpy.cos(omega*t)/numpy.tanh(H_BAR_ARU*omega/(2.0*KB_ARU*T_eid))

        def Kt(t, T):
            return integrate.quad(Kt_integrand, 0.0, 4.0*omega_c_eid, args=(t, T_eid))[0]

        Kt = numpy.vectorize(Kt)
        Kt_t_T = Kt(t_eid, T_eid)
        array_len_eid = int((Omega_end_eid - Omega_start_eid)/Omega_step_eid)
        Omega_step_eid = (Omega_end_eid - Omega_start_eid)/array_len_eid
        Omega_eid = numpy.linspace(Omega_start_eid, Omega_end_eid, array_len_eid)


        def K_integrand_real(t, Omega, T):
            Kt = numpy.interp(t, t_eid, Kt_t_T)
            return numpy.cos(Omega*t)*Kt

        def K_integrand_imag(t, Omega, T):
            Kt = numpy.interp(t, t_eid, Kt_t_T)
            return numpy.sin(Omega*t)*Kt

        def K_real(Omega, T):
            #print convert(Omega, ARU_TO_EV)
            return integrate.quad(K_integrand_real, 0.0, t_end_eid, args=(Omega, T))[0]

        def K_imag(Omega, T):
            #print convert(Omega, ARU_TO_EV)
            return integrate.quad(K_integrand_imag, 0.0, t_end_eid, args=(Omega, T))[0]

        K_real = numpy.vectorize(K_real)
        K_imag = numpy.vectorize(K_imag)

        K_real_eid = K_real(Omega_eid, T_eid)
        K_imag_eid = K_imag(Omega_eid, T_eid)

        #alpha = convert(0.0272, PICO2_TO_ARU)
        #K_real_eid = (PI/2.0)*alpha*pow(Omega_eid, 3.0)*numpy.exp(-pow(Omega_eid/omega_c_eid, 2.0))/numpy.tanh(Omega_eid/(2.0*KB_ARU*T_eid))

        print '...calculation complete.'
        numpy.savetxt('data/eid/Omega_eid.txt', numpy.transpose(1000.0*convert(Omega_eid, ARU_TO_EV)), delimiter=',', header='hbar Omega (meV)')
        numpy.savetxt('data/eid/K_real.txt', numpy.transpose(convert(K_real_eid, ARU_TO_INV_PICO)), delimiter=',', header='Re[K(hbar Omega)] (ps-1)')
        numpy.savetxt('data/eid/K_imag.txt', numpy.transpose(1000.0*convert(K_imag_eid, ARU_TO_EV)), delimiter=',', header='Im[hbar K(hbar Omega)] (meV)')
       # plt.subplot(2,1,1)
       # plt.plot(1000.0*convert(Omega_eid, ARU_TO_EV), convert(K_real_eid, ARU_TO_INV_PICO), '-', label='K real')
       # plt.xlabel('Omega (meV)')
       # plt.ylabel('Re[$K(\hbar\Omega)$] (ps$^-1$)')
       # plt.subplot(2,1,2)
       # plt.plot(1000.0*convert(Omega_eid, ARU_TO_EV), 1000.0*convert(K_imag_eid, ARU_TO_EV), '-', label='K imag')
       # plt.xlabel('Rabi Energy $\hbar\Omega$ (meV)')
       # plt.ylabel('Im[$\hbar K(\hbar\Omega)$] (meV)')
#        plt.show(block=False)
       # plt.close()
        params['laphnoneid_params'].Omega_eid = Omega_eid
        params['laphnoneid_params'].K_real_eid = K_real_eid
        params['laphnoneid_params'].K_imag_eid = K_imag_eid











