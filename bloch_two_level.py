#! /home/ajan/anaconda/bin/python
# Filename: pulse.py

from matplotlib.pyplot import *
from numpy import *
from pulse import *
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import quad
from scipy import special
import sys
from scipy import interp
import datetime
import ctypes as ct


from pulse import *
from convert import *
from read_config import *
from phasemasks import *

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

        self.d_yo = params[dot_param_object].qdot_d_yo

        self.rwa = params['run_params'].run_rwa

        if (params['run_params'].run_dephasing == True):
            self.gamma_yy = 1.0/params[dot_param_object].qdot_T_yy
            self.gamma_yo = 1.0/params[dot_param_object].qdot_T_yo
        else:
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

        self.array_len = len(elec_time)
        self.efield_phase = angle(efield)
        self.efield_inst_freq = zeros(self.array_len)
        #self.efield_inst_freq[0:self.array_len-1] = diff(unwrap(self.efield_phase))/diff(elec_time)
        #self.efield_inst_freq[self.array_len-1] = self.efield_inst_freq[self.array_len-2]
        self.efield_inst_freq[0:self.array_len-1] = diff(unwrap(self.efield_phase))/self.dt
        self.efield_inst_freq[self.array_len-1] = self.efield_inst_freq[self.array_len-2]
        #print len(self.efield_inst_freq)
        #print len(self.elec_time)
        #print len(efield)
        #print len(self.efield_phase)

        if (self.rwa == True):
            efield = efield*exp(-1j*self.pulse_omega*(self.elec_time - self.pulse_delay))

        self.elec_field_real = efield.real
        self.elec_field_imag = efield.imag

        self.phi_2 = params['mask_params'].x[0]
        self.alpha = 2.0*self.phi_2/(pow(self.pulse_width, 4.0)/pow(2.0*log(2.0), 2.0) + pow(2.0*self.phi_2, 2.0))

        self.run_dephasing = params['run_params'].run_dephasing
        self.run_wl_eid = params['run_params'].run_wl_eid
        self.run_phonon_eid = params['run_params'].run_phonon_eid
        self.run_nonmarkovian_eid = params['run_params'].run_nonmarkovian_eid

        if self.run_phonon_eid == True:
            # eid parameters from config file
            self.dot_shape = params['laphnoneid_params'].dot_shape
            self.omega_c_eid = params['laphnoneid_params'].omega_c_eid
            self.omega_z_eid = params['laphnoneid_params'].omega_z_eid
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
            peak_intensity = max(I)
            t1_index = find_nearest_index(I, 0.5*peak_intensity)
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
            #self.gamma_wl_eid = self.c*pow(pulse_area, 2.0)/(self.tl_fwhm*self.chirped_fwhm)
            self.gamma_wl_eid = self.c*pow(pulse_area/self.tl_fwhm, 2.0)*(1/self.chirped_fwhm)
        else:
            self.gamma_wl_eid = 0.0

        if self.run_nonmarkovian_eid:
            self.kappa_0 = params['nonmarkovian_params'].kappa_0
            self.kappa_1 = params['nonmarkovian_params'].kappa_1
            self.kappa_2 = params['nonmarkovian_params'].kappa_2
        else:
            self.kappa_0 = 0.0
            self.kappa_1 = 0.0
            self.kappa_2 = 0.0
        '''
        # print type of dephasing
        if self.run_phonon_eid and not(self.run_wl_eid or self.run_dephasing):
            print "dephasing model: LA phonons only"
        elif self.run_wl_eid and not(self.run_phonon_eid or self.run_dephasing):
            print "dephasing model: wetting layer only"
        elif self.run_phonon_eid and self.run_dephasing and not(self.run_wl_eid):
            print "dephasing model: LA phonons and constant dephasing"
        elif self.run_wl_eid and self.run_dephasing and not(self.run_phonon_eid):
            print "dephasing model: wetting layer and constant dephasing"
        elif self.run_dephasing and not(self.run_wl_eid or self.run_phonon_eid):
            print "dephasing model: constant dephasing only"
        elif not(self.run_wl_eid or self.run_phonon_eid or self.run_dephasing or self.run_nonmarkovian_eid):
            print "dephasing model: no dephasing"
        elif self.run_nonmarkovian_eid:
            print "dephasing model: non markovian"
        else:
            print "dephasing model: INCOMPATABILE"
        '''
    def operator(self, t, s):

        ''' functions to extract la phonon dephasing params from kernel '''
        def Kr(Gamma):
#            # real part of kernel
            if self.dot_shape == 'SPHERICAL':  
                return (PI/2.0)*self.alpha_eid*Gamma*Gamma*Gamma*exp(-(Gamma*Gamma/(self.omega_c_eid*self.omega_c_eid)))/tanh(Gamma/(2.0*KB_ARU*self.T_eid)) #harmonic oscillator ground state
            elif self.dot_shape == 'PANCAKE':
                def F(theta,omega):
                    return exp(-1*pow(omega,2)*(pow(sin(theta)/self.omega_c_eid,2)+pow(cos(theta)/(self.omega_z_eid),2)))*sin(theta)  
                return (PI/2.0)*self.alpha_eid*pow(Gamma, 3.0)*quad(F,0,PI,args=(Gamma))[0]/tanh(Gamma/(2.0*KB_ARU*self.T_eid))
            elif self.dot_shape == 'PANCAKE_C':
                lib = ct.CDLL("./lib.so")
                N = ct.c_int(1000)
                omega_c = ct.c_double(self.omega_c_eid)
                omega_l = ct.c_double(self.omega_z_eid)
                lib.integrate.argtypes = [ct.c_int, ct.c_double, ct.c_double, ct.c_double]
                lib.integrate.restype = ct.c_double
                return 0.5*(PI/2.0)*self.alpha_eid*lib.integrate(N, ct.c_double(Gamma), omega_c, omega_l)/tanh(Gamma/(2.0*KB_ARU*self.T_eid))
            elif self.dot_shape == 'P_shell0':
                return (PI/2.0)*self.alpha_eid*pow(Gamma, 3.0)*exp(-pow(Gamma/self.omega_c_eid, 2.0))*(1+(-2/3.0)*pow(Gamma/self.omega_c_eid, 2.0)+(1/5.0)*pow(Gamma/self.omega_c_eid, 4.0))/tanh(Gamma/(2.0*KB_ARU*self.T_eid))
            elif self.dot_shape == 'P_shell1':
                return (PI/2.0)*self.alpha_eid*pow(Gamma, 3.0)*exp(-pow(Gamma/self.omega_c_eid, 2.0))*(1+(-2/3.0)*pow(Gamma/self.omega_c_eid, 2.0)+(2/15.0)*pow(Gamma/self.omega_c_eid, 4.0))/tanh(Gamma/(2.0*KB_ARU*self.T_eid))
            elif self.dot_shape == 'PANCAKE_P_shell0':
                def F(theta,omega):
                    st=sin(theta)
                    ct=cos(theta)
                    return exp(-1*pow(omega,2)*(pow(st/self.omega_c_eid,2)+pow(ct/(self.omega_z_eid),2)))*st*(1-(2*pow(omega/self.omega_z_eid,2)*pow(ct,2))+pow(omega/self.omega_z_eid,4)*pow(ct,4))  
                return (PI/2.0)*self.alpha_eid*pow(Gamma, 3.0)*quad(F,0,PI,args=(Gamma))[0]/tanh(Gamma/(2.0*KB_ARU*self.T_eid))
            elif self.dot_shape == 'PANCAKE_P_shell0_C':
                lib = ct.CDLL("./lib.so")
                N = ct.c_int(1000)
                omega_c = ct.c_double(self.omega_c_eid)
                omega_l = ct.c_double(self.omega_z_eid)
                lib.integratepshell.argtypes = [ct.c_int, ct.c_double, ct.c_double, ct.c_double]
                lib.integratepshell.restype = ct.c_double
                return 0.5*(PI/2.0)*self.alpha_eid*lib.integratepshell(N, ct.c_double(Gamma), omega_c, omega_l)/tanh(Gamma/(2.0*KB_ARU*self.T_eid))
                
            
            
        def Ki(Gamma):
            # imaginary part of kernel
            if Gamma < max(self.Omega_eid):
                return 0.0 #interp(Gamma, self.Omega_eid, self.K_imag_eid)
            else:
                #print 'Gamma out of bounds: {0:.3f} meV'.format(1000.0*convert(Gamma, ARU_TO_EV))
                return 0.0

        def J(Gamma):
            # la-phonon - electron interaction spectrum
            if self.dot_shape == 'SPHERICAL':
                return self.alpha_eid*Gamma*Gamma*Gamma*exp(-(Gamma*Gamma/(self.omega_c_eid*self.omega_c_eid))) #harmonic oscillator ground state
            elif self.dot_shape == 'PANCAKE':
                def F(theta,omega):
                    return exp(-1*pow(omega,2)*(pow(sin(theta)/self.omega_c_eid,2)+pow(cos(theta)/(self.omega_z_eid),2)))*sin(theta) 
                return self.alpha_eid*pow(Gamma, 3.0)*quad(F,0,PI,args=(Gamma))[0]
            elif self.dot_shape == 'PANCAKE_C':
                lib = ct.CDLL("./lib.so")
                N = ct.c_int(1000)
                omega_c = ct.c_double(self.omega_c_eid)
                omega_l = ct.c_double(self.omega_z_eid)
                lib.integrate.argtypes = [ct.c_int, ct.c_double, ct.c_double, ct.c_double]
                lib.integrate.restype = ct.c_double
                return 0.5*self.alpha_eid*lib.integrate(N, ct.c_double(Gamma), omega_c, omega_l)
            elif self.dot_shape == 'P_shell0':
                return self.alpha_eid*pow(Gamma, 3.0)*exp(-pow(Gamma/self.omega_c_eid, 2.0))*(1+(-2/3.0)*pow(Gamma/self.omega_c_eid, 2.0)+(1/5.0)*pow(Gamma/self.omega_c_eid, 4.0))
            elif self.dot_shape == 'P_shell1':
                return self.alpha_eid*pow(Gamma, 3.0)*exp(-pow(Gamma/self.omega_c_eid, 2.0))*(1+(-2/3.0)*pow(Gamma/self.omega_c_eid, 2.0)+(2/15.0)*pow(Gamma/self.omega_c_eid, 4.0))
            elif self.dot_shape == 'PANCAKE_P_shell0':
                def F(theta,omega):
                    st=sin(theta)
                    ct=cos(theta)
                    return exp(-1*pow(omega,2)*(pow(st/self.omega_c_eid,2)+pow(ct/(self.omega_z_eid),2)))*st*(1-(2*pow(omega/self.omega_z_eid,2)*pow(ct,2))+pow(omega/self.omega_z_eid,4)*pow(ct,4))
                return self.alpha_eid*pow(Gamma, 3.0)*quad(F,0,PI,args=(Gamma))[0]
            elif self.dot_shape == 'PANCAKE_P_shell0_C':
                lib = ct.CDLL("./lib.so")
                N = ct.c_int(1000)
                omega_c = ct.c_double(self.omega_c_eid)
                omega_l = ct.c_double(self.omega_z_eid)
                lib.integratepshell.argtypes = [ct.c_int, ct.c_double, ct.c_double, ct.c_double]
                lib.integratepshell.restype = ct.c_double
                return 0.5*self.alpha_eid*lib.integratepshell(N, ct.c_double(Gamma), omega_c, omega_l)
#            
             
        '''Returns derivatives of rho at time t'''
        # linear interpolation to find electric field at time t
        #efield_real_interp = interp(t, self.elec_time, self.elec_field_real)
        #efield_imag_interp = interp(t, self.elec_time, self.elec_field_imag)
        #inst_freq = interp(t, self.elec_time, self.efield_inst_freq)

        # replaces interpolation functions above - results in 1.88X speed-up of code
        i = int(numpy.floor((t - self.elec_time[0])/self.dt))
        efield_real_interp = 0.0
        efield_imag_interp = 0.0
        if t < self.elec_time[-1]:
            efield_real_interp = self.elec_field_real[i] + (self.elec_field_real[i+1] - self.elec_field_real[i])*(t-self.elec_time[i])/self.dt
            efield_imag_interp = self.elec_field_imag[i] + (self.elec_field_imag[i+1] - self.elec_field_imag[i])*(t-self.elec_time[i])/self.dt
            inst_freq_interp = self.efield_inst_freq[i] + (self.efield_inst_freq[i+1] - self.efield_inst_freq[i])*(t-self.elec_time[i])/self.dt

        efield_interp = complex(efield_real_interp, efield_imag_interp)

        # determine instantaneous Rabi frequencies
        Omega = self.d_yo * abs(efield_interp)
        if Omega == 0.0:
            Omega = sys.float_info.epsilon
        Delta = self.omega_yo - (self.pulse_omega + 2.0*self.alpha*(t-self.pulse_delay))
        #Delta1 = Delta
        #Delta = self.omega_yo - (inst_freq_interp)
        #print 1000*convert(Delta1, ARU_TO_EV), 1000*convert(Delta, ARU_TO_EV)
        Gamma = pow(Omega**2 + Delta**2, 0.5)

        # assign density matrix elements
        sx = s[0]
        sy = s[1]
        sz = s[2]

        # calculate derivatives
        if self.run_phonon_eid and not(self.run_wl_eid or self.run_dephasing):
            # phonon mechanism
            dsxdt = Delta*sy - (Omega/Gamma)**2*Kr(Gamma)*sx - Delta*Omega*Kr(Gamma)*sz/Gamma**2 - pi*Omega*J(Gamma)/(2.0*Gamma)
            dsydt = -Delta*sx + Omega*sz + (Omega*Ki(Gamma)/Gamma)*sz - (Omega/Gamma)**2*Kr(Gamma)*sy
            dszdt = -Omega*sy
        elif self.run_wl_eid and not(self.run_phonon_eid or self.run_dephasing):
            # wetting layer mechanism
            dsxdt = Delta*sy - self.gamma_wl_eid*sx
            dsydt = -Delta*sx + Omega*sz - self.gamma_wl_eid*sy
            dszdt = -Omega*sy
        elif self.run_phonon_eid and self.run_dephasing and not(self.run_wl_eid):
            # phonon mechanism and pure dephasing
            dsxdt = Delta*sy - (Omega/Gamma)**2*Kr(Gamma)*sx - Delta*Omega*Kr(Gamma)*sz/Gamma**2 - pi*Omega*J(Gamma)/(2.0*Gamma) - self.gamma_yo*sx
            dsydt = -Delta*sx + Omega*sz + (Omega*Ki(Gamma)/Gamma)*sz - (Omega/Gamma)**2*Kr(Gamma)*sy - self.gamma_yo*sy
            dszdt = -Omega*sy - self.gamma_yy*sz
        elif self.run_wl_eid and self.run_dephasing and not(self.run_phonon_eid):
            # wetting layer mechanism and pure dephasing
            dsxdt = Delta*sy - (self.gamma_yo + self.gamma_wl_eid)*sx
            dsydt = -Delta*sx + Omega*sz - (self.gamma_yo + self.gamma_wl_eid)*sy
            dszdt = -Omega*sy - self.gamma_yy*sz
        elif self.run_dephasing and not(self.run_wl_eid or self.run_phonon_eid):
            # pure dephasing
            dsxdt = Delta*sy - self.gamma_yo*sx
            dsydt = -Delta*sx + Omega*sz - self.gamma_yo*sy
            dszdt = -Omega*sy - self.gamma_yy*sz
        elif not(self.run_wl_eid or self.run_phonon_eid or self.run_dephasing or self.run_nonmarkovian_eid):
            # no dephasing
            dsxdt = Delta*sy
            dsydt = -Delta*sx + Omega*sz
            dszdt = -Omega*sy
        elif self.run_nonmarkovian_eid:
            dsxdt = Delta*sy - (self.kappa_0 + self.kappa_1*Omega + self.kappa_2*pow(Omega, 2.0))*sx
            dsydt = -Delta*sx + Omega*sz - (self.kappa_0 + self.kappa_1*Omega + self.kappa_2*pow(Omega, 2.0))*sy
            dszdt = -Omega*sy
	elif self.run_phonon_eid and self.run_wl_eid and not (self.run_dephasing):
            # phonon mechanism and wetting layer mechanism
            dsxdt = Delta*sy - (Omega/Gamma)**2*Kr(Gamma)*sx - Delta*Omega*Kr(Gamma)*sz/Gamma**2 - pi*Omega*J(Gamma)/(2.0*Gamma) - self.gamma_wl_eid*sx
            dsydt = -Delta*sx + Omega*sz + (Omega*Ki(Gamma)/Gamma)*sz - (Omega/Gamma)**2*Kr(Gamma)*sy - self.gamma_wl_eid*sy
            dszdt = -Omega*sy
        else:
            print "CHECK DEPHASING SETTINGS - BLOCH_TWO_LEVEL.PY"
            dsxdt = Delta*sy
            dsydt = -Delta*sx + Omega*sz
            dszdt = -Omega*sy
        dsdt = zeros(3)
        # assign derivatives
        dsdt[0] = dsxdt
        dsdt[1] = dsydt
        dsdt[2] = dszdt

        return dsdt
