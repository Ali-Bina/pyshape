#! /home/ajan/anaconda/bin/python
# Filename: phasemasks.py
import numpy

def array_to_dmatrix(rho):
    rho_oo = rho[0]
    rho_xx = rho[1]
    rho_yy = rho[2]
    rho_bb = rho[3]
    rho_xo = complex(rho[4], rho[5])
    rho_ox = numpy.conjugate(rho_xo)
    rho_yo = complex(rho[6], rho[7])
    rho_oy = numpy.conjugate(rho_yo)
    rho_yx = complex(rho[8], rho[9])
    rho_xy = numpy.conjugate(rho_yx)
    rho_bx = complex(rho[10], rho[11])
    rho_xb = numpy.conjugate(rho_bx)
    rho_by = complex(rho[12], rho[13])
    rho_yb = numpy.conjugate(rho_by)
    rho_bo = complex(rho[14], rho[15])
    rho_ob = numpy.conjugate(rho_bo)

    dmatrix = numpy.zeros((4,4), dtype=complex)

    dmatrix[0,0] = rho_oo
    dmatrix[1,0] = rho_xo
    dmatrix[2,0] = rho_yo
    dmatrix[3,0] = rho_bo

    dmatrix[0,1] = rho_ox
    dmatrix[1,1] = rho_xx
    dmatrix[2,1] = rho_yx
    dmatrix[3,1] = rho_bx

    dmatrix[0,2] = rho_oy
    dmatrix[1,2] = rho_xy
    dmatrix[2,2] = rho_yy
    dmatrix[3,2] = rho_by

    dmatrix[0,3] = rho_ob
    dmatrix[1,3] = rho_xb
    dmatrix[2,3] = rho_yb
    dmatrix[3,3] = rho_bb

    return dmatrix

def dmatrix_to_array(dmatrix):
    rho_oo = dmatrix[0,0]
    rho_xo = dmatrix[1,0]
    rho_yo = dmatrix[2,0]
    rho_bo = dmatrix[3,0]

    rho_ox = dmatrix[0,1]
    rho_xx = dmatrix[1,1]
    rho_yx = dmatrix[2,1]
    rho_bx = dmatrix[3,1]

    rho_oy = dmatrix[0,2]
    rho_xy = dmatrix[1,2]
    rho_yy = dmatrix[2,2]
    rho_by = dmatrix[3,2]

    rho_ob = dmatrix[0,3]
    rho_xb = dmatrix[1,3]
    rho_yb = dmatrix[2,3]
    rho_bb = dmatrix[3,3]

    rho = numpy.zeros(16)

    rho[0] = rho_oo.real
    rho[1] = rho_xx.real
    rho[2] = rho_yy.real
    rho[3] = rho_bb.real
    rho[4] = rho_xo.real
    rho[5] = rho_xo.imag
    rho[6] = rho_yo.real
    rho[7] = rho_yo.imag
    rho[8] = rho_yx.real
    rho[9] = rho_yx.imag
    rho[10] = rho_bx.real
    rho[11] = rho_bx.imag
    rho[12] = rho_by.real
    rho[13] = rho_by.imag
    rho[14] = rho_bo.real
    rho[15] = rho_bo.imag

    return rho

def state_to_dmatrix(state):
    '''convert row vector with amplitude coefficients to density matrix'''
    dmatrix = numpy.outer(state, numpy.conj(state))
    return dmatrix

def state_to_blochvector(state):
    ''' convert row vector with amplitude coefficients to bloch vector'''
    co = state[0]             # ground state
    cy = state[2]             # Y state (transition energies in config file are defined with respect to this state)

    s1 = (co*numpy.conj(cy) + numpy.conj(co)*cy).real
    s2 = (-1j*(co*numpy.conj(cy) - numpy.conj(co)*cy)).real
    s3 = (pow(abs(cy), 2.0) - pow(abs(co), 2.0)).real

    bloch_vector  = numpy.array([s1, s2, s3])

    return bloch_vector


def blochvector_to_state(S):
    '''convert bloch vector to row vector with amplitude coefficients'''
    state = numpy.zeros((2,), dtype=complex)

    # express in terms of {coC1*, |c1|^2} or {|co|^2, co*C1} and renormalize
    # the if-else statement is necessary so that we don't get a null vector
    inversion = (S[2] + 1.0)/2.0
    if inversion != 0.0:
        state[0] = (S[0] + 1j*S[1])/2.0
        state[1] = (S[2] + 1.0)/2.0
    else:
        state[0] = (1.0 - S[2])/2.0
        state[1] = (S[0] - 1j*S[1])/2.0

    state = normalize_state(state)
    return state

def normalize_state(state):
    '''normalize any given state'''
    N = len(state)
    state_norm = 0.0
    for i in range(N):
        state_norm = state_norm + pow(abs(state[i]), 2.0)
    state_norm = pow(state_norm, 0.5)
    state = state/state_norm
    return state
