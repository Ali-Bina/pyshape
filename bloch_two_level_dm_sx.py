#!/usr/bin/python
# Filename: pulse.py

import numpy

#from numpy import *
from pulse import *
import scipy.interpolate

from pulse import *
from convert import *
from read_config import *
from phasemasks import *
from constants import *

from time import sleep
import sys

class BlochEqns(object):

    def __init__(self, params, efield, elec_time, dot_number):

        if dot_number==DOT1:
            dot_param_object = 'qdot_params'
        elif dot_number==DOT2:
            dot_param_object = 'qdot2_params'
        elif dot_number==DOT3:
            dot_param_object = 'qdot3_params'

        self.omega_yo = params[dot_param_object].qdot_omega_yo
        self.d_yo = params[dot_param_object].qdot_d_yo

        self.rwa = params['run_params'].run_rwa

        if (params['run_params'].run_dephasing == True):
            self.gamma_oo = 1.0/params[dot_param_object].qdot_T_oo
            self.gamma_yy = 1.0/params[dot_param_object].qdot_T_yy
            self.gamma_yo = 1.0/params[dot_param_object].qdot_T_yo
            
        else:
            self.gamma_oo = 0.0
            self.gamma_yy = 0.0
            self.gamma_yo = 0.0

        self.pulse_omega = params['pulse_params'].pulse_omega_o
        self.pulse_width = params['pulse_params'].pulse_width
        self.xcomp = params['pulse_params'].pulse_xcomp
        self.ycomp = params['pulse_params'].pulse_ycomp
        self.pulse_delay = params['pulse_params'].pulse_delay
        self.pulse_phase = params['pulse_params'].pulse_phase
        self.elec_time = elec_time
        self.dt = self.elec_time[1]-self.elec_time[0]

        if (self.rwa == True):
            efield = efield*numpy.exp(-1j*self.pulse_omega*(self.elec_time - self.pulse_delay))

        self.elec_field_real = efield.real
        self.elec_field_imag = efield.imag
        self.run_wl_eid = params['run_params'].run_wl_eid
        self.run_phonon_eid = params['run_params'].run_phonon_eid

        if self.run_phonon_eid:
            # eid parameters from config file
            self.omega_c_eid = params['laphnoneid_params'].omega_c_eid
            self.alpha_eid = params['laphnoneid_params'].alpha_eid
            self.Omega_start_eid = params['laphnoneid_params'].Omega_start_eid
            self.Omega_end_eid = params['laphnoneid_params'].Omega_end_eid
            self.Omega_step_eid = params['laphnoneid_params'].Omega_step_eid
            self.T_eid = params['laphnoneid_params'].T_eid

            self.Omega_eid = params['laphnoneid_params'].Omega_eid
            self.K_real_eid = params['laphnoneid_params'].K_real_eid
            self.K_imag_eid = params['laphnoneid_params'].K_imag_eid

        self.drhodt = numpy.zeros(3)

    def operator(self, t, rho):
        '''Returns derivatives of rho at time t'''
        # linear interpolation to find electric field at time t
        i = numpy.floor((t - self.elec_time[0])/self.dt)
        efield_real_interp = 0.0
        efield_imag_interp = 0.0
        if t < self.elec_time[-1]:
            efield_real_interp = self.elec_field_real[i] + (self.elec_field_real[i+1] - self.elec_field_real[i])*(t-self.elec_time[i])/self.dt
            efield_imag_interp = self.elec_field_imag[i] + (self.elec_field_imag[i+1] - self.elec_field_imag[i])*(t-self.elec_time[i])/self.dt

        efield_interp = complex(efield_real_interp, efield_imag_interp)
        efield_interp_x = self.xcomp*efield_interp
        efield_interp_y = self.ycomp*efield_interp

        # determine instantaneous Rabi frequencies
        mu_yo = -self.d_yo * efield_interp_y
        mu_oy = numpy.conj(mu_yo)

        # assign density matrix elements
      	rho_oo = 1.0 - rho[0]
        rho_yy = rho[0]
        rho_yo = complex(rho[1], rho[2])
        rho_oy = numpy.conj(rho_yo)


        ''' functions to extract la phonon dephasing params from kernel '''
        def Kr(Gamma):
            # real part of kernel
            return (PI/2.0)*self.alpha_eid*pow(Gamma, 3.0)*numpy.exp(-pow(Gamma/self.omega_c_eid, 2.0))/numpy.tanh(Gamma/(2.0*KB_ARU*self.T_eid))

        def J(Gamma):
            # la-phonon - electron interaction spectrum
            return self.alpha_eid*pow(Gamma, 3.0)*numpy.exp(-pow(Gamma/self.omega_c_eid, 2.0))

        Omega = mu_yo
        Delta = self.omega_yo - self.pulse_omega
        Gamma = pow(abs(Omega)**2 + Delta**2, 0.5)


        if (self.rwa is True):
            # calculate derivatives in rotating wave approximation
            # phonon mechanism
            sx = 2*rho_yo.real
            sy = -2*rho_yo.imag
            sz = 2*rho_yy - 1

            if self.run_phonon_eid:
                dsxdt = -Delta*sy - Omega.imag*sz + (abs(Omega)/Gamma)**2*Kr(Gamma)*sx + Delta*abs(Omega)*Kr(Gamma)*sz/Gamma**2 + PI*abs(Omega)*J(Gamma)/(2.0*Gamma)
                dsydt = Delta*sx - Omega.real*sz + (abs(Omega)/Gamma)**2*Kr(Gamma)*sy
                #dsydt = Delta*sx - Omega.real*sz - (abs(Omega)*Ki(Gamma)/Gamma)*sz + (abs(Omega)/Gamma)**2*Kr(Gamma)*sy
                dszdt = Omega.imag*sx + Omega.real*sy
            else:
                dsxdt = -Delta*sy - Omega.imag*sz
                dsydt = Delta*sx - Omega.real*sz
                dszdt = Omega.imag*sx + Omega.real*sy
            drho_yy = complex(dszdt/2)
            drho_yo = complex(dsxdt/2, -1*dsydt/2)
        else:
            # full differential equations
            drho_yy = -(1j/2.0)*(mu_yo*rho_oy - mu_oy*rho_yo) - self.gamma_yy*rho_yy
            drho_yo = -(1j/2.0)*(2.0*self.omega_yo*rho_yo + mu_yo*(rho_oo - rho_yy)) - self.gamma_yo*rho_yo


        # assign derivatives
        self.drhodt[0] = drho_yy.real
        self.drhodt[1] = drho_yo.real
        self.drhodt[2] = drho_yo.imag


        return self.drhodt

