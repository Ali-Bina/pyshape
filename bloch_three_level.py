#! /home/ajan/anaconda/bin/python
# Filename: pulse.py

from matplotlib.pyplot import *
from numpy import *
from pulse import *
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import simps
import sys
from scipy import interp
import datetime


from pulse import *
from convert import *
from read_config import *
from phasemasks import *
from constants import *

class BlochEqns(object):

    def __init__(self, params, efield, elec_time, dot_number):

        if dot_number==DOT1:
            dot_param_object = 'qdot_params'
        elif dot_number==DOT2:
            dot_param_object = 'qdot2_params'
        elif dot_number==DOT3:
            dot_param_object = 'qdot3_params'

        self.omega_yo = params[dot_param_object].qdot_omega_yo
        self.omega_yx = params[dot_param_object].qdot_omega_xyfine
        self.omega_bind = params[dot_param_object].qdot_omega_bind
        self.omega_xo = self.omega_yo - self.omega_yx
        self.omega_bo = self.omega_xo + self.omega_yo - self.omega_bind
        self.omega_bx = self.omega_yo - self.omega_bind
        self.omega_by = self.omega_xo - self.omega_bind

        self.d_xo = params[dot_param_object].qdot_d_xo
        self.d_yo = params[dot_param_object].qdot_d_yo
        self.d_xy = params[dot_param_object].qdot_d_xy
        self.d_bo = params[dot_param_object].qdot_d_bo
        self.d_bx = params[dot_param_object].qdot_d_bx
        self.d_by = params[dot_param_object].qdot_d_by

        self.rwa = params['run_params'].run_rwa

        if (params['run_params'].run_dephasing == True):
            self.gamma_oo = 1.0/params[dot_param_object].qdot_T_oo
            self.gamma_xx = 1.0/params[dot_param_object].qdot_T_xx
            self.gamma_yy = 1.0/params[dot_param_object].qdot_T_yy
            self.gamma_bb = 1.0/params[dot_param_object].qdot_T_bb
            self.gamma_xo = 1.0/params[dot_param_object].qdot_T_xo
            self.gamma_yo = 1.0/params[dot_param_object].qdot_T_yo
            self.gamma_bo = 1.0/params[dot_param_object].qdot_T_bo
            self.gamma_xy = 1.0/params[dot_param_object].qdot_T_xy
            self.gamma_bx = 1.0/params[dot_param_object].qdot_T_bx
            self.gamma_by = 1.0/params[dot_param_object].qdot_T_by
        else:
            self.gamma_oo = 0.0
            self.gamma_xx = 0.0
            self.gamma_yy = 0.0
            self.gamma_bb = 0.0
            self.gamma_xo = 0.0
            self.gamma_yo = 0.0
            self.gamma_bo = 0.0
            self.gamma_xy = 0.0
            self.gamma_bx = 0.0
            self.gamma_by = 0.0

        self.pulse_omega = params['pulse_params'].pulse_omega_o
        self.pulse_width = params['pulse_params'].pulse_width
        self.xcomp = params['pulse_params'].pulse_xcomp
        self.ycomp = params['pulse_params'].pulse_ycomp
        self.pulse_delay = params['pulse_params'].pulse_delay
        self.pulse_phase = params['pulse_params'].pulse_phase
        self.elec_time = elec_time
        if (self.rwa == True):
            efield = efield*exp(-1j*self.pulse_omega*(self.elec_time - self.pulse_delay))

        self.elec_field_real = efield.real
        self.elec_field_imag = efield.imag
        self.run_wl_eid = params['run_params'].run_wl_eid


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
        mu_oy = conj(mu_yo)


        # assign density matrix elements
        rho_oo = rho[0]
        rho_xx = rho[1]
        rho_yy = rho[2]

        rho_yo = complex(rho[3], rho[4])

        rho_oy = conj(rho_yo)

        if (self.rwa is True):
            # calculate derivatives in rotating wave approximation
            drho_xx = self.gamma_yy*rho_yy
            drho_yy = (1j/2.0)*(mu_oy*rho_yo - mu_yo*rho_oy) - self.gamma_yy*rho_yy
            drho_oo = -(1j/2.0)*(mu_oy*rho_yo - mu_yo*rho_oy)
            drho_yo = -(1j/2.0)*(2.0*(self.omega_yo - self.pulse_omega)*rho_yo + mu_yo*(rho_oo - rho_yy)) - (self.gamma_yy/2.0 + self.gamma_wl_eid)*rho_yo
        else:
            # full differential equations
            drho_xx = -(1j/2.0)*(mu_xo*rho_ox - mu_ox*rho_xo) - self.gamma_xx*rho_xx + (self.gamma_bb/2)*rho_bb
            drho_yy = -(1j/2.0)*(mu_yo*rho_oy - mu_oy*rho_yo) - self.gamma_yy*rho_yy + (self.gamma_bb/2)*rho_bb
            drho_oo = -(1j/2.0)*(mu_ox*rho_xo - mu_xo*rho_ox - mu_yo*rho_oy + mu_oy*rho_yo) + self.gamma_xx*rho_xx + self.gamma_yy*rho_yy
            drho_yo = -(1j/2.0)*(2.0*self.omega_yo*rho_yo - mu_xo*rho_yx + mu_yo*(rho_oo - rho_yy)) - self.gamma_yo*rho_yo


        drhodt = zeros(5)

        # assign derivatives
        drhodt[0] = drho_oo.real
        drhodt[1] = drho_xx.real
        drhodt[2] = drho_yy.real

        drhodt[3] = drho_yo.real
        drhodt[4] = drho_yo.imag




        return drhodt

