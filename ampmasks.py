#! /home/ajan/anaconda/bin/python
# Filename: phasemasks.py
import numpy
from constants import *
from convert import *


def noAmpMask(omega, omega_o, amp_mask, params, **kwargs):
    return amp_mask

# only useful for the CROT gate
def chenAmpMask(omega, omega_o, amp_mask, params, **kwargs):
    # if optimizing, use kwargs values from optimization routine, else read from config file
    if params['run_params'].run_mask == 'ampmask_chen':
        p1width = kwargs['p1width']
        p2width = kwargs['p2width']
        aratio = kwargs['Aratio']
    else:
        p1width = params['ampmask_params'].x[params['ampmask_params'].param_list.index('p1width')]
        p2width = params['ampmask_params'].x[params['ampmask_params'].param_list.index('p2width')]
        aratio = params['ampmask_params'].x[params['ampmask_params'].param_list.index('Aratio')]

    # define transition energies
    omega_yo =  params['qdot_params'].qdot_omega_yo
    omega_xy =  params['qdot_params'].qdot_omega_xyfine
    omega_bind = params['qdot_params'].qdot_omega_bind
    omega_xo = omega_yo - omega_xy
    omega_by = omega_xo - omega_bind

    p1omega = omega_by
    p2omega = omega_yo

    # calculate amplitude mask
    num1 = pow(omega - p1omega, 2.0)
    den1 = 8.0*numpy.log(2.0)/pow(p1width, 2.0)
    num2 = pow(omega - p2omega, 2.0)
    den2 = 8.0*numpy.log(2.0)/pow(p2width, 2.0)
    shape1 = numpy.exp(-num1/den1)
    shape2 = numpy.exp(-num2/den2)
    amp_mask = abs(shape1 - aratio*shape2)

    return amp_mask

def slmccAmpMask(omega, omega_o, phase_mask, params, **kwargs):
    # if optimizing, use kwargs values from optimization routine, else read from config file
    if params['run_params'].run_mask == 'mask_slmcc':
        g0102 = kwargs['g0102']
        T2 = kwargs['T2']
        phi_2 = kwargs['phi_2']
    else:
        g0102 = params['ampmask_params'].x[params['ampmask_params'].param_list.index('g0102')]
        T2 = params['ampmask_params'].x[params['ampmask_params'].param_list.index('T2')]
        phi_2 = params['ampmask_params'].x[params['ampmask_params'].param_list.index('phi_2')]

    phi = T2*(abs(omega) - omega_o) - phi_2
    # calculate phase mask
    amp_mask = pow(1.0 + pow(g0102, 2.0) + 2*g0102*numpy.cos(phi), 0.5) / max(pow(1.0 + pow(g0102, 2.0) + 2*g0102*numpy.cos(phi), 0.5))
    amp_mask = pow(1.0 + pow(g0102, 2.0) + 2*g0102*numpy.cos(phi), 0.5) / pow(1.0 + pow(g0102, 2.0) + 2*g0102, 0.5)

    return amp_mask

def dichromaticMask(omega,omega_o,amp_mask, params,**kwargs):
    #test implimentation of dichromatic pulses
    if params['run_params'].run_mask == 'dichrome':
        delta=kwargs['delta']
        
    else:
        delta = params['ampmask_params'].x[params['ampmask_params'].param_list.index('delta')]
    #setting position of spectral hole to dot transition
    omega_o = params['qdot_params'].qdot_omega_yo
    #setting the depth of hole to pulse amplitude
    amplitude = params['ampmask_params'].x[params['ampmask_params'].param_list.index('pulse_area')]
    
    amp_mask= 1-(numpy.exp(-2.0*numpy.log(2.0)*pow((omega-omega_o), 2.0 )/pow(delta, 2.0)))
    
    return amp_mask

