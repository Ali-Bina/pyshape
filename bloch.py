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
        self.dt = self.elec_time[1]-self.elec_time[0]

        if (self.rwa == True):
            efield = efield*numpy.exp(-1j*self.pulse_omega*(self.elec_time - self.pulse_delay))

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

        self.drhodt = numpy.zeros(16)

    def operator(self, t, rho):
        '''Returns derivatives of rho at time t'''
        # linear interpolation to find electric field at time t
        i = int(numpy.floor((t - self.elec_time[0])/self.dt))
        efield_real_interp = 0.0
        efield_imag_interp = 0.0
        if t < self.elec_time[-1]:
            efield_real_interp = self.elec_field_real[i] + (self.elec_field_real[i+1] - self.elec_field_real[i])*(t-self.elec_time[i])/self.dt
            efield_imag_interp = self.elec_field_imag[i] + (self.elec_field_imag[i+1] - self.elec_field_imag[i])*(t-self.elec_time[i])/self.dt

        efield_interp = complex(efield_real_interp, efield_imag_interp)
        efield_interp_x = self.xcomp*efield_interp
        efield_interp_y = self.ycomp*efield_interp

        # determine instantaneous Rabi frequencies
        mu_xo = -self.d_xo * efield_interp_x
        mu_yo = -self.d_yo * efield_interp_y
        mu_bo = -self.d_bo * efield_interp
        mu_xy = -self.d_xy * efield_interp
        mu_bx = -self.d_bx * efield_interp_x
        mu_by = -self.d_by * efield_interp_y
        mu_ox = numpy.conj(mu_xo)
        mu_oy = numpy.conj(mu_yo)
        mu_ob = numpy.conj(mu_bo)
        mu_yx = numpy.conj(mu_xy)
        mu_xb = numpy.conj(mu_bx)
        mu_yb = numpy.conj(mu_by)

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

        if (self.rwa is True):
            # calculate derivatives in rotating wave approximation
            drho_xx = -(1j/2.0)*(mu_xo*rho_ox - mu_ox*rho_xo - mu_bx*rho_xb + mu_xb*rho_bx) - self.gamma_xx*rho_xx + (self.gamma_bb/2.0)*rho_bb
            drho_yy = -(1j/2.0)*(mu_yo*rho_oy - mu_oy*rho_yo - mu_by*rho_yb + mu_yb*rho_by) - self.gamma_yy*rho_yy + (self.gamma_bb/2.0)*rho_bb
            drho_bb = -(1j/2.0)*(mu_bx*rho_xb - mu_xb*rho_bx - mu_yb*rho_by + mu_by*rho_yb) - self.gamma_bb*rho_bb
            drho_oo = -(1j/2.0)*(mu_ox*rho_xo - mu_xo*rho_ox - mu_yo*rho_oy + mu_oy*rho_yo) + self.gamma_xx*rho_xx + self.gamma_yy*rho_yy
            drho_xo = -(1j/2.0)*(2.0*(self.omega_xo - self.pulse_omega)*rho_xo + mu_xb*rho_bo - mu_yo*rho_xy + mu_xo*(rho_oo - rho_xx)) - (self.gamma_xo + self.gamma_wl_eid)*rho_xo
            drho_yo = -(1j/2.0)*(2.0*(self.omega_yo - self.pulse_omega)*rho_yo + mu_yb*rho_bo - mu_xo*rho_yx + mu_yo*(rho_oo - rho_yy)) - (self.gamma_yo + self.gamma_wl_eid)*rho_yo
            drho_yx = -(1j/2.0)*(2.0*self.omega_yx*rho_yx + mu_yo*rho_ox - mu_ox*rho_yo + mu_yb*rho_bx - mu_bx*rho_yb) - (self.gamma_xy + self.gamma_wl_eid)*rho_xy
            drho_bo = -(1j/2.0)*(2.0*(self.omega_bo - 2.0*self.pulse_omega)*rho_bo + mu_bx*rho_xo - mu_xo*rho_bx + mu_by*rho_yo - mu_yo*rho_by) - (self.gamma_bo + self.gamma_wl_eid)*rho_bo
            drho_bx = -(1j/2.0)*(2.0*(self.omega_bx - self.pulse_omega)*rho_bx + mu_by*rho_yx - mu_ox*rho_bo + mu_bx*(rho_xx - rho_bb)) - (self.gamma_bx + self.gamma_wl_eid)*rho_bx
            drho_by = -(1j/2.0)*(2.0*(self.omega_by - self.pulse_omega)*rho_by + mu_bx*rho_xy - mu_oy*rho_bo + mu_by*(rho_yy - rho_bb)) - (self.gamma_by + self.gamma_wl_eid)*rho_by
        else:
            # full differential equations
            drho_xx = -(1j/2.0)*(mu_xo*rho_ox - mu_ox*rho_xo - mu_bx*rho_xb + mu_xb*rho_bx) - self.gamma_xx*rho_xx + (self.gamma_bb/2.0)*rho_bb
            drho_yy = -(1j/2.0)*(mu_yo*rho_oy - mu_oy*rho_yo - mu_by*rho_yb + mu_yb*rho_by) - self.gamma_yy*rho_yy + (self.gamma_bb/2.0)*rho_bb
            drho_bb = -(1j/2.0)*(mu_bx*rho_xb - mu_xb*rho_bx - mu_yb*rho_by + mu_by*rho_yb) - self.gamma_bb*rho_bb
            drho_oo = -(1j/2.0)*(mu_ox*rho_xo - mu_xo*rho_ox - mu_yo*rho_oy + mu_oy*rho_yo) + self.gamma_xx*rho_xx + self.gamma_yy*rho_yy
            drho_xo = -(1j/2.0)*(2.0*self.omega_xo*rho_xo + mu_xb*rho_bo - mu_yo*rho_xy + mu_xo*(rho_oo - rho_xx)) - self.gamma_xo*rho_xo
            drho_yo = -(1j/2.0)*(2.0*self.omega_yo*rho_yo + mu_yb*rho_bo - mu_xo*rho_yx + mu_yo*(rho_oo - rho_yy)) - self.gamma_yo*rho_yo
            drho_yx = -(1j/2.0)*(2.0*self.omega_yx*rho_yx + mu_yo*rho_ox - mu_ox*rho_yo + mu_yb*rho_bx - mu_bx*rho_yb) - self.gamma_xy*rho_xy
            drho_bo = -(1j/2.0)*(2.0*self.omega_bo*rho_bo + mu_bx*rho_xo - mu_xo*rho_bx + mu_by*rho_yo - mu_yo*rho_by) - self.gamma_bo*rho_bo
            drho_bx = -(1j/2.0)*(2.0*self.omega_bx*rho_bx + mu_by*rho_yx - mu_ox*rho_bo + mu_bx*(rho_xx - rho_bb)) - self.gamma_bx*rho_bx
            drho_by = -(1j/2.0)*(2.0*self.omega_by*rho_by + mu_bx*rho_xy - mu_oy*rho_bo + mu_by*(rho_yy - rho_bb)) - self.gamma_by*rho_by





        # assign derivatives
        self.drhodt[0] = drho_oo.real
        self.drhodt[1] = drho_xx.real
        self.drhodt[2] = drho_yy.real
        self.drhodt[3] = drho_bb.real

        self.drhodt[4] = drho_xo.real
        self.drhodt[5] = drho_xo.imag
        self.drhodt[6] = drho_yo.real
        self.drhodt[7] = drho_yo.imag
        self.drhodt[8] = drho_yx.real
        self.drhodt[9] = drho_yx.imag
        self.drhodt[10] = drho_bx.real
        self.drhodt[11] = drho_bx.imag
        self.drhodt[12] = drho_by.real
        self.drhodt[13] = drho_by.imag
        self.drhodt[14] = drho_bo.real
        self.drhodt[15] = drho_bo.imag

        return self.drhodt

