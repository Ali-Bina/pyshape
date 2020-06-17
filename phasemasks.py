#! /home/ajan/anaconda/bin/python
# Filename: phasemasks.py
import numpy
from constants import *
from convert import *

def noPhaseMask(omega, omega_o, phase_mask, params, **kwargs):
    return phase_mask


def cosinePhaseMask(omega, omega_o, phase_mask, params, **kwargs):
    # if optimizing, use kwargs values from optimization routine, else read from config file
    if params['run_params'].run_mask == 'phasemask_cos':
        phi_alpha = kwargs['phi_alpha']
        phi_gamma = kwargs['phi_gamma']
        phi_delta = kwargs['phi_delta']
    else:
        phi_alpha = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_alpha')]
        phi_gamma = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_gamma')]
        phi_delta = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_delta')]

    # calculate phase mask
    phase_mask = phi_alpha*numpy.cos(phi_gamma*(abs(omega) - omega_o) - phi_delta)

    return phase_mask

def polyPhaseMask(omega, omega_o, phase_mask, params, **kwargs):
    # if optimizing, use kwargs values from optimization routine, else read from config file
    if params['run_params'].run_mask == 'phasemask_poly':
        phi_2 = kwargs['phi_2']
        phi_3 = kwargs['phi_3']
        phi_4 = kwargs['phi_4']
    else:
        phi_2 = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_2')]
        phi_3 = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_3')]
        phi_4 = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_4')]


    # calculate phase mask
    phase_mask =phi_2*pow((abs(omega) - omega_o), 2.0)/2.0 + phi_3*pow((abs(omega) - omega_o), 3.0)/6.0 + phi_4*pow((abs(omega) - omega_o), 4.0)/24.0

    return phase_mask

def slmccPhaseMask(omega, omega_o, phase_mask, params, **kwargs):
    # if optimizing, use kwargs values from optimization routine, else read from config file
    if params['run_params'].run_mask == 'mask_slmcc':
        g0102 = kwargs['g0102']
        T2 = kwargs['T2']
        phi_2 = kwargs['phi_2']
    else:
        g0102 = params['phasemask_params'].x[params['phasemask_params'].param_list.index('g0102')]
        T2 = params['phasemask_params'].x[params['phasemask_params'].param_list.index('T2')]
        phi_2 = params['phasemask_params'].x[params['phasemask_params'].param_list.index('phi_2')]

    phi = T2*(abs(omega) - omega_o) - phi_2
    # calculate phase mask
    phase_mask = numpy.arctan( g0102*numpy.sin(phi)/(1.0 + g0102*numpy.cos(phi)) )

    return phase_mask
