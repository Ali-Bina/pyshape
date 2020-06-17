#! /home/ajan/anaconda/bin/python
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
	self.phi_2 = params['mask_params'].x[0]
	self.alpha = 2.0*self.phi_2/(pow(self.pulse_width, 4.0)/pow(2.0*numpy.log(2.0), 2.0) + pow(2.0*self.phi_2, 2.0))

	if self.run_phonon_eid == True:
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


        if self.run_wl_eid == True:
            def find_nearest_value(array,value):
                idx = (abs(array-value)).argmin()
                return array[idx]
            def find_nearest_index(array,value):
                idx = (abs(array-value)).argmin()
                return idx

            # find full-width half max of chirped pulse
            I = pow(abs(efield).real, 2.0)
            max_value = max(I)
            t1_index = find_nearest_index(I, 0.5*max_value)
            t1 = elec_time[t1_index]
            tmax = elec_time[-1]
            self.chirped_fwhm = 2.0*(abs(tmax/2.0 - t1))
            self.tl_fwhm = self.pulse_width

            # determine the pulse area
            dipole = convert(self.d_yo, ARU_TO_DEBYE)                       # has to be in debye for conversion
            width = convert(self.pulse_width, ARU_TO_FEMTO)                 # had to be in femtoseconds for conversion
            shape = params['pulse_params'].pulse_shape
            kwargs = {'dipole_moment': dipole, 'pulse_width': width, 'pulse_shape': shape}
            pulse_area = convert(params['pulse_params'].pulse_EO, ARU_TO_AREA, **kwargs)

            # wl dephasing constant
            self.c = params[dot_param_object].qdot_eid_c
            self.b = params[dot_param_object].qdot_eid_b

            print 'the pulse area is: {0:.1f} PI radians'.format(pulse_area)
            print "the tl pulse width is: {0:.3f} fs".format(convert(self.tl_fwhm, ARU_TO_FEMTO))
            print "the chirped pulse width is: {0:.3f} fs".format(convert(self.chirped_fwhm, ARU_TO_FEMTO))
            # note: pulse area is in units of PI radians
            self.gamma_wl_eid = self.c*pow(pulse_area, 2.0)/(self.chirped_fwhm*self.chirped_fwhm) + self.b*pulse_area/self.chirped_fwhm
        else:
            self.gamma_wl_eid = 0.0
        self.idk = 0
        self.t_old = 0.0
    def operator(self, t, rho):
	
	''' functions to extract la phonon dephasing params from kernel '''
        def Kr(Gamma):
            # real part of kernel
            return (PI/2.0)*self.alpha_eid*pow(Gamma, 3.0)*numpy.exp(-pow(Gamma/self.omega_c_eid, 2.0))/numpy.tanh(Gamma/(2.0*KB_ARU*self.T_eid))

        def Ki(Gamma):
            # imaginary part of kernel
            if Gamma < max(self.Omega_eid):
                return 0.0 #interp(Gamma, self.Omega_eid, self.K_imag_eid)
            else:
                #print 'Gamma out of bounds: {0:.3f} meV'.format(1000.0*convert(Gamma, ARU_TO_EV))
                return 0.0

        def J(Gamma):
            # la-phonon - electron interaction spectrum
            return self.alpha_eid*pow(Gamma, 3.0)*numpy.exp(-pow(Gamma/self.omega_c_eid, 2.0))

        '''Returns derivatives of rho at time t'''
        # linear interpolation to find electric field at time t
        jk = 1
        if jk == 1:
            efield_real_interp = numpy.interp(t, self.elec_time, self.elec_field_real)
            efield_imag_interp = numpy.interp(t, self.elec_time, self.elec_field_imag)
        else:
            i = numpy.floor((t - self.elec_time[0])/self.dt)
            efield_real_interp = 0.0
            efield_imag_interp = 0.0
            if t < self.elec_time[-1]:
                efield_real_interp = self.elec_field_real[i] + (self.elec_field_real[i+1] - self.elec_field_real[i])*(t-self.elec_time[i])/self.dt
                efield_imag_interp = self.elec_field_imag[i] + (self.elec_field_imag[i+1] - self.elec_field_imag[i])*(t-self.elec_time[i])/self.dt

        self.idk = self.idk + 1.0
        #print self.idk, convert(t - self.t_old, ARU_TO_FEMTO), convert(t, ARU_TO_FEMTO)
        self.t_old = t
        #print self.idk, convert(t, ARU_TO_FEMTO), t/self.elec_time[-1], numpy.interp(t, self.elec_time, self.elec_field_real), efield_real_interp, 100*(numpy.interp(t, self.elec_time, self.elec_field_real) - efield_real_interp)/(efield_real_interp + sys.float_info.epsilon)
        #sleep(0.01)

        efield_interp = complex(efield_real_interp, efield_imag_interp)
        efield_interp_x = self.xcomp*efield_interp
        efield_interp_y = self.ycomp*efield_interp

        # determine instantaneous Rabi frequencies
        Omega = self.d_yo * abs(efield_interp)
        mu_yo = -self.d_yo * efield_interp_y
        mu_oy = numpy.conj(mu_yo)
        Delta = self.omega_yo - (self.pulse_omega + 2.0*self.alpha*(t-self.pulse_delay))
        Gamma = pow(Omega**2 + Delta**2, 0.5)


        # assign density matrix elements
        rho_oo = rho[0]
        rho_xx = rho[1]
        rho_yy = rho[2]
        rho_bb = rho[3]

        rho_xo = complex(rho[4], rho[5])
        rho_yo = complex(rho[6], rho[7])
        rho_yx = complex(rho[8], rho[9])
        rho_bx = complex(rho[10], rho[11])
        rho_by = complex(rho[12], rho[13])
        rho_bo = complex(rho[14], rho[15])
        rho_ox = numpy.conj(rho_xo)
        rho_oy = numpy.conj(rho_yo)
        rho_xy = numpy.conj(rho_yx)
        rho_xb = numpy.conj(rho_bx)
        rho_yb = numpy.conj(rho_by)
        rho_ob = numpy.conj(rho_bo)

        if (self.rwa is True) and not(self.run_phonon_eid) :
            # calculate derivatives in rotating wave approximation
            drho_yy = -(1j/2.0)*(mu_yo*rho_oy- mu_oy*rho_yo) - self.gamma_yy*rho_yy
            drho_oo = -(1j/2.0)*(-mu_yo*rho_oy+ mu_oy*rho_yo) + self.gamma_yy*rho_yy
            drho_yo = -(1j/2.0)*(2.0*(self.omega_yo - self.pulse_omega)*rho_yo + mu_yo*(rho_oo - rho_yy)) - (self.gamma_yo + self.gamma_wl_eid)*rho_yo
        elif (self.rwa is False) and not(self.run_phonon_eid) :
            # full differential equations
            drho_yy = -(1j/2.0)*(mu_yo*(rho_oy-rho_yo) - mu_oy*(rho_oy-rho_yo)) - self.gamma_yy*rho_yy
            drho_oo = -(1j/2.0)*(-mu_yo*(rho_oy-rho_yo) + mu_oy*(rho_oy-rho_yo)) + self.gamma_yy*rho_yy
            drho_yo = -(1j/2.0)*(2.0*self.omega_yo*rho_yo + mu_yo*(rho_oo - rho_yy)+mu_oy*(rho_oo - rho_yy)) - self.gamma_yo*rho_yo
        elif (self.rwa is True) and (self.run_phonon_eid) :
            drho_yy = -(1j/2.0)*(mu_yo*rho_oy- mu_oy*rho_yo) - self.gamma_yy*rho_yy
            drho_oo = -(1j/2.0)*(-mu_yo*rho_oy+ mu_oy*rho_yo) + self.gamma_yy*rho_yy
            drho_yo = -(1j/2.0)*(2.0*(self.omega_yo - self.pulse_omega)*rho_yo + (1+(Kr(Gamma)/Gamma))*mu_yo*(rho_oo - rho_yy)) + numpy.pi*Omega*J(Gamma)/(2.0*Gamma)+ Delta*Omega*Kr(Gamma)*(rho_oo - rho_yy)/Gamma**2-2*(Omega/Gamma)**2*Kr(Gamma)*rho_yo

#        else:
#            # full differential equations
#            drho_yy = -(1j/2.0)*(mu_yo*rho_oy - mu_oy*rho_yo) - self.gamma_yy*rho_yy
#            drho_oo = -(1j/2.0)*(-mu_yo*rho_oy + mu_oy*rho_yo) + self.gamma_yy*rho_yy
#            drho_yo = -(1j/2.0)*(2.0*self.omega_yo*rho_yo + mu_yo*(rho_oo - rho_yy)) - self.gamma_yo*rho_yo


        drhodt = numpy.zeros(16)

        # assign derivatives
        drhodt[0] = drho_oo.real
        drhodt[1] = 0.0
        drhodt[2] = drho_yy.real
        drhodt[3] = 0.0

        drhodt[4] = 0.0
        drhodt[5] = 0.0
        drhodt[6] = drho_yo.real
        drhodt[7] = drho_yo.imag
        drhodt[8] = 0.0
        drhodt[9] = 0.0
        drhodt[10] = 0.0
        drhodt[11] = 0.0
        drhodt[12] = 0.0
        drhodt[13] = 0.0
        drhodt[14] = 0.0
        drhodt[15] = 0.0

        return drhodt

