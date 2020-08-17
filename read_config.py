from constants import *
from convert import *
from matmanip import *
import numpy
from configobj import ConfigObj
from validate import Validator
import sys
import cmath

# class to hold run parameters
class RunParams(object):
    def __init__(self, run_t_start=None, run_t_end=None, run_t_step=None, run_dephasing=None, run_phonon_eid = None, run_wl_eid = None, run_nonmarkovian_eid=None, run_rwa=None, run_mask=None, run_ampmask=None, run_phasemask=None, NITER=None, gate=None, optimize=None, sobol_seed=None, show_plot=None):
        self.run_t_start = run_t_start
        self.run_t_end = run_t_end
        self.run_t_step = run_t_step
        self.run_dephasing = run_dephasing
        self.run_phonon_eid = run_phonon_eid
        self.run_wl_eid = run_wl_eid
        self.run_nonmarkovian_eid = run_nonmarkovian_eid
        self.run_rwa = run_rwa
        self.run_mask = run_mask
        self.run_ampmask = run_ampmask
        self.run_phasemask = run_phasemask
        self.NITER = NITER
        self.gate = gate
        self.optimize = optimize
        self.sobol_seed = sobol_seed
        self.show_plot = show_plot

# class to hold runge-kutta integration parameters
class RkdpParams(object):
    def __init__(self, rkdp_atol=None, rkdp_rtol=None, rkdp_nsteps_opt=None, rkdp_nsteps_tdep=None, rkdp_hstart=None, rkdp_hmax=None, rkdp_tdep_step=None):
        self.rkdp_atol = rkdp_atol
        self.rkdp_rtol = rkdp_rtol
        self.rkdp_nsteps_opt = rkdp_nsteps_opt
        self.rkdp_nsteps_tdep = rkdp_nsteps_tdep
        self.rkdp_hstart = rkdp_hstart
        self.rkdp_hmax = rkdp_hmax
        self.rkdp_tdep_step = rkdp_tdep_step

# class to hold mask parameters
class MaskParams(object):
    def __init__(self, param_list=None, x=None, x_bounds=None, x_convert=None, x_units=None, NOPTS=None):
        self.param_list = param_list
        self.x = x
        self.x_bounds = x_bounds
        self.x_convert = x_convert
        self.x_units = x_units
        self.NOPTS = NOPTS

# class to hold mask parameters
class AmpMaskParams(object):
    def __init__(self, param_list=None, x=None, x_bounds=None, x_convert=None, x_units=None, NOPTS=None):
        self.param_list = param_list
        self.x = x
        self.x_bounds = x_bounds
        self.x_convert = x_convert
        self.x_units = x_units
        self.NOPTS = NOPTS

# class to hold mask parameters
class PhaseMaskParams(object):
    def __init__(self, param_list=None, x=None, x_bounds=None, x_convert=None, x_units=None, NOPTS=None):
        self.param_list = param_list
        self.x = x
        self.x_bounds = x_bounds
        self.x_convert = x_convert
        self.x_units = x_units
        self.NOPTS = NOPTS

# class to hold pulse parameters
class PulseParams(object):
    def __init__(self, pulse_omega_o=None, pulse_detun=None, pulse_width=None, pulse_delay=None, pulse_chirp=None, pulse_shape=None, pulse_pol=None, pulse_phase=None, pulse_area=None, pulse_EO=None, pulse_dipole=None, pulse_xcomp=None, pulse_ycomp=None, pulse_hole_width=None):
        self.pulse_omega_o = pulse_omega_o
        self.pulse_detun = pulse_detun
        self.pulse_width = pulse_width
        self.pulse_delay = pulse_delay
        self.pulse_chirp = pulse_chirp
        self.pulse_shape = pulse_shape
        self.pulse_pol = pulse_pol
        self.pulse_phase = pulse_phase
        self.pulse_area = pulse_area
        self.pulse_EO = pulse_EO
        self.pulse_dipole = pulse_dipole
        self.pulse_xcomp = pulse_xcomp
        self.pulse_ycomp = pulse_ycomp
        self.pulse_hole_width = pulse_hole_width

# class to hold quantum dot parameters
class QDotParams(object):
    def __init__(self, qdot_omega_yo=None, qdot_omega_xyfine=None, qdot_omega_bind=None, qdot_d_xo=None, qdot_d_yo=None, qdot_d_xy=None, qdot_d_bo=None, qdot_d_bx=None, qdot_d_by=None, qdot_T_oo=None, qdot_T_xx=None, qdot_T_yy=None, qdot_T_bb=None, qdot_T_xo=None, qdot_T_yo=None, qdot_T_bo=None, qdot_T_xy=None, qdot_T_bx=None, qdot_T_by=None, qdot_eid_c=None, qdot_eid_b=None, qdot_initial_state=None, qdot_desired_state=None):
        self.qdot_omega_yo = qdot_omega_yo
        self.qdot_omega_xyfine = qdot_omega_xyfine
        self.qdot_omega_bind = qdot_omega_bind
        self.qdot_d_xo = qdot_d_xo
        self.qdot_d_yo = qdot_d_yo
        self.qdot_d_xy = qdot_d_xy
        self.qdot_d_bo = qdot_d_bo
        self.qdot_d_bx = qdot_d_bx
        self.qdot_d_by = qdot_d_by
        self.qdot_T_oo = qdot_T_oo
        self.qdot_T_xx = qdot_T_xx
        self.qdot_T_yy = qdot_T_yy
        self.qdot_T_bb = qdot_T_bb
        self.qdot_T_xo = qdot_T_xo
        self.qdot_T_yo = qdot_T_yo
        self.qdot_T_bo = qdot_T_bo
        self.qdot_T_xy = qdot_T_xy
        self.qdot_T_bx = qdot_T_bx
        self.qdot_T_by = qdot_T_by
        self.qdot_eid_c = qdot_eid_c
        self.qdot_eid_b = qdot_eid_b
        self.qdot_initial_state = qdot_initial_state
        self.qdot_desired_state = qdot_desired_state

# class to hold la-phonon eid parameters
class LAPhononEIDParams(object):
    def __init__(self, omega_c_eid=None, alpha_eid=None, Omega_start_eid=None, Omega_end_eid=None, Omega_step_eid=None, T_eid=None, read_kernel=None, Omega_eid=None, K_real_eid=None, K_imag_eid=None):
        self.omega_c_eid = omega_c_eid
        self.alpha_eid = alpha_eid
        self.Omega_start_eid = Omega_start_eid
        self.Omega_end_eid = Omega_end_eid
        self.Omega_step_eid = Omega_step_eid
        self.T_eid = T_eid
        self.read_kernel = read_kernel
        if Omega_eid is None:
            self.Omega_eid = []
        else:
            self.Omega_eid = Omega_eid
        if K_real_eid is None:
            self.K_real_eid = []
        else:
            self.K_real_eid = K_real_eid
        if K_imag_eid is None:
            self.K_imag_eid = []
        else:
            self.K_imag_eid = K_imag_eid

# class to hold non-markovian eid parameters
class NonMarkovianParams(object):
    def __init__(self, kappa_0=None, kappa_1=None, kappa_2=None):
        self.kappa_0 = kappa_0
        self.kappa_1 = kappa_1
        self.kappa_2 = kappa_2

def read_config (config_file, configspec_file):

    # read in configuration file and validate against configspec
    config = ConfigObj(config_file, configspec='configspec.ini')
    validator = Validator()
    result = config.validate(validator)

    if result != True:
        print "Error reading config file"
        sys.exit(1)

    #------------------------------------------------------------#


    # read run parameters
    dump_to_screen = config['run']['dump']
    run_t_start = config['run']['t_start']
    run_t_end = config['run']['t_end']
    run_t_step = config['run']['t_step']
    run_dephasing = config['run']['dephasing']
    run_phonon_eid = config['run']['phonon_eid']
    run_wl_eid = config['run']['wl_eid']
    run_nonmarkovian_eid = config['run']['non_markovian_eid']
    run_rwa = config['run']['rwa']
    run_mask = config['run']['mask']
    run_ampmask = config['run']['ampmask']
    run_phasemask = config['run']['phasemask']
    NITER = config['run']['NITER']
    gate = config['run']['gate']
    optimize = config['run']['optimize']
    sobol_seed = config['run']['sobol_seed']
    show_plot = config['run']['show_plot']
    two_level_bool = config['run']['two_level']

    if optimize is True and run_mask != run_ampmask and run_mask != run_phasemask:
        sys.exit("The mask does not match the ampmask or the phasemask - fix config file.\nProgram exiting now!")
    if (gate == 'twodottwoleveldm') and (two_level_bool != True):
        print "Change two_level boolean to True if using twodottwoleveldm mask"
        quit()
    elif (gate != 'twodottwoleveldm') and (two_level_bool == True):
        print "Change two_level boolean to False if not using twodottwoleveldm mask"
        quit()
        
    # write parameters to screen if requested
    if dump_to_screen == True:
        print "reading run parameters..."
        print "run start time:", run_t_start, "fs"
        print "run end time:", run_t_end, "fs"
        print "run time step:", run_t_step, "fs"
        print "dephasing:", run_dephasing
        print "phonon eid:", run_phonon_eid
        print "wl eid:", run_wl_eid
        print "non markovian eid:", run_nonmarkovian_eid
        print "rwa:", run_rwa
        print "mask:", run_mask
        print "NITER:", NITER
        print "gate", gate, "\n"

    # convert run pararmeters to Atomic Rydberg Units (ARU)
    run_t_start = convert(run_t_start, FEMTO_TO_ARU)    # convert run start time from fs to ARU
    run_t_end = convert(run_t_end, FEMTO_TO_ARU)        # convert run start time from fs to ARU
    run_t_step = convert(run_t_step, FEMTO_TO_ARU)        # convert run time step from fs to ARU

    run_params = RunParams(run_t_start, run_t_end, run_t_step, run_dephasing, run_phonon_eid, run_wl_eid, run_nonmarkovian_eid, run_rwa, run_mask, run_ampmask, run_phasemask, NITER, gate, optimize, sobol_seed, show_plot)

    #------------------------------------------------------------#

    # read rkdp parameters
    rkdp_atol = config['rkdp']['atol']
    rkdp_rtol = config['rkdp']['rtol']
    rkdp_nsteps_opt = config['rkdp']['nsteps_opt']
    rkdp_nsteps_tdep = config['rkdp']['nsteps_tdep']
    rkdp_hstart = config['rkdp']['h_start']
    rkdp_hmax = config['rkdp']['h_max']
    rkdp_tdep_step = config['rkdp']['tdep_step']

    # write parameters to screen if requested
    if dump_to_screen == True:
        print 'reading rkdp parameters...'
        print "rkdp absolute tolerance:", rkdp_atol
        print "rkdp relative tolerance:", rkdp_rtol
        print "rkdp max steps opt:", rkdp_nsteps_opt
        print "rkdp max steps tdep:", rkdp_nsteps_tdep
        print "rkdp initial step size:", rkdp_hstart
        print "rkdp max steps size:", rkdp_hmax
        print "tdep step:", rkdp_tdep_step, "\n"

    # convert time dependent step from fs to ARU
    rkdp_hstart = convert(rkdp_hstart, FEMTO_TO_ARU)
    rkdp_hmax = convert(rkdp_hmax, FEMTO_TO_ARU)
    rkdp_tdep_step = convert(rkdp_tdep_step, FEMTO_TO_ARU)
    rkdp_params = RkdpParams(rkdp_atol, rkdp_rtol, rkdp_nsteps_opt, rkdp_nsteps_tdep, rkdp_hstart, rkdp_hmax, rkdp_tdep_step)

    #------------------------------------------------------------#
    # read in parameters for chosen mask
    param_list = config['masks'][run_mask]['param_list']
    NOPTS = len(param_list)
    # read in parameters and create tuple with bounds for optimization
    # note parameters are in the same as param_list
    x = numpy.zeros(NOPTS)        # values of parameters
    x_lower = numpy.zeros(NOPTS)        # values of parameters
    x_upper = numpy.zeros(NOPTS)        # values of parameters
    x_bounds = ()                        # tuple with bounds for optimization
    x_convert = []
    x_units = []

    # read in values and convert to ARU
    for i in range(NOPTS):
        # read in parameters
        x[i] = config['masks'][run_mask][param_list[i]]
        param_ARU = param_list[i] + '_ARU'
        x_units.append(config['masks'][run_mask][param_ARU][2])

        # read in bounds
        param_range = param_list[i] + '_range'
        x_lower[i] = config['masks'][run_mask][param_range][0]
        x_upper[i] = config['masks'][run_mask][param_range][1]

    # write parameters to screen if requested
    if dump_to_screen == True:
        print 'reading ' + run_mask + ' parameters...'
        for i in range(NOPTS):
            print param_list[i], ': Value:', x[i], ', Bounds:', x_lower[i], ',', x_upper[i], x_units[i]
        print  "\n"

    # convert parameters according to conversion factor
    for i in range(NOPTS):
        # read in conversion factor
        param_ARU = param_list[i] + '_ARU'
        x_convert.append(config['masks'][run_mask][param_ARU][0])

        # if conversion is from AREA (in fraction of PI) to ARU, pass extra kwargs
        if (x_convert[i] == AREA_TO_ARU):
            dipole = config['pulse']['dipole']
            width = config['pulse']['width']
            shape = config['pulse']['shape']
            kwargs = {'dipole_moment': dipole, 'pulse_width': width, 'pulse_shape': shape}
        else:
            kwargs = {}

        x[i] = convert(x[i], x_convert[i], **kwargs)
        x_lower[i] = convert(x_lower[i], x_convert[i], **kwargs)
        x_upper[i] = convert(x_upper[i], x_convert[i], **kwargs)

        # add bounds for i-th variable to tuple
        x_bounds = x_bounds + (numpy.array([x_lower[i], x_upper[i]]), )

    mask_params = MaskParams(param_list, x, x_bounds, x_convert, x_units, NOPTS)

    #------------------------------------------------------------#
    # read in parameters for ampmask
    param_list = config['masks'][run_ampmask]['param_list']
    NOPTS = len(param_list)
    # read in parameters and create tuple with bounds for optimization
    # note parameters are in the same as param_list
    x = numpy.zeros(NOPTS)        # values of parameters
    x_lower = numpy.zeros(NOPTS)        # values of parameters
    x_upper = numpy.zeros(NOPTS)        # values of parameters
    x_bounds = ()                        # tuple with bounds for optimization
    x_convert = []
    x_units = []

    # read in values and convert to ARU
    for i in range(NOPTS):
        # read in parameters
        x[i] = config['masks'][run_ampmask][param_list[i]]
        param_ARU = param_list[i] + '_ARU'
        x_units.append(config['masks'][run_ampmask][param_ARU][2])

        # read in bounds
        param_range = param_list[i] + '_range'
        x_lower[i] = config['masks'][run_ampmask][param_range][0]
        x_upper[i] = config['masks'][run_ampmask][param_range][1]

    # write parameters to screen if requested
    if dump_to_screen == True:
        print 'reading ' + run_ampmask + ' parameters...'
        for i in range(NOPTS):
            print param_list[i], ': Value:', x[i], ', Bounds:', x_lower[i], ',', x_upper[i], x_units[i]
        print  "\n"

    # convert parameters according to conversion factor
    for i in range(NOPTS):
        # read in conversion factor
        param_ARU = param_list[i] + '_ARU'
        x_convert.append(config['masks'][run_ampmask][param_ARU][0])

        # if conversion is from AREA (in fraction of PI) to ARU, pass extra kwargs
        if (x_convert[i] == AREA_TO_ARU):
            dipole = config['pulse']['dipole']
            width = config['pulse']['width']
            shape = config['pulse']['shape']
            kwargs = {'dipole_moment': dipole, 'pulse_width': width, 'pulse_shape': shape}
        else:
            kwargs = {}

        x[i] = convert(x[i], x_convert[i], **kwargs)
        x_lower[i] = convert(x_lower[i], x_convert[i], **kwargs)
        x_upper[i] = convert(x_upper[i], x_convert[i], **kwargs)

        # add bounds for i-th variable to tuple
        x_bounds = x_bounds + (numpy.array([x_lower[i], x_upper[i]]), )

    ampmask_params = AmpMaskParams(param_list, x, x_bounds, x_convert, x_units, NOPTS)

    #------------------------------------------------------------#
    # read in parameters for phasemask
    param_list = config['masks'][run_phasemask]['param_list']
    NOPTS = len(param_list)
    # read in parameters and create tuple with bounds for optimization
    # note parameters are in the same as param_list
    x = numpy.zeros(NOPTS)        # values of parameters
    x_lower = numpy.zeros(NOPTS)        # values of parameters
    x_upper = numpy.zeros(NOPTS)        # values of parameters
    x_bounds = ()                        # tuple with bounds for optimization
    x_convert = []
    x_units = []

    # read in values and convert to ARU
    for i in range(NOPTS):
        # read in parameters
        x[i] = config['masks'][run_phasemask][param_list[i]]
        param_ARU = param_list[i] + '_ARU'
        x_units.append(config['masks'][run_phasemask][param_ARU][2])

        # read in bounds
        param_range = param_list[i] + '_range'
        x_lower[i] = config['masks'][run_phasemask][param_range][0]
        x_upper[i] = config['masks'][run_phasemask][param_range][1]

    # write parameters to screen if requested
    if dump_to_screen == True:
        print 'reading ' + run_phasemask + ' parameters...'
        for i in range(NOPTS):
            print param_list[i], ': Value:', x[i], ', Bounds:', x_lower[i], ',', x_upper[i], x_units[i]
        print  "\n"

    # convert parameters according to conversion factor
    for i in range(NOPTS):
        # read in conversion factor
        param_ARU = param_list[i] + '_ARU'
        x_convert.append(config['masks'][run_phasemask][param_ARU][0])

        # if conversion is from AREA (in fraction of PI) to ARU, pass extra kwargs
        if (x_convert[i] == AREA_TO_ARU):
            dipole = config['pulse']['dipole']
            width = config['pulse']['width']
            shape = config['pulse']['shape']
            kwargs = {'dipole_moment': dipole, 'pulse_width': width, 'pulse_shape': shape}
        else:
            kwargs = {}

        x[i] = convert(x[i], x_convert[i], **kwargs)
        x_lower[i] = convert(x_lower[i], x_convert[i], **kwargs)
        x_upper[i] = convert(x_upper[i], x_convert[i], **kwargs)

        # add bounds for i-th variable to tuple
        x_bounds = x_bounds + (numpy.array([x_lower[i], x_upper[i]]), )

    phasemask_params = PhaseMaskParams(param_list, x, x_bounds, x_convert, x_units, NOPTS)
    #------------------------------------------------------------#

    # read pulse parameters
    pulse_omega_o = config['pulse']['omega_o']
    pulse_detun = config['pulse']['detun']
    pulse_width = config['pulse']['width']
    pulse_delay = config['pulse']['delay']
    pulse_chirp = config['pulse']['chirp']
    pulse_shape = config['pulse']['shape']
    pulse_pol = config['pulse']['pol']
    pulse_phase = config['pulse']['phase']
    pulse_area = config['pulse']['area']
    pulse_dipole = config['pulse']['dipole']
    pulse_hole_width = config['pulse']['hole_width']

    # write parameters to screen if requested
    if dump_to_screen == True:
        print "reading pulse parameters..."
        print "pulse parameters:"
        print "pulse center freq is:", pulse_omega_o, "eV"
        print "pulse shape is:", pulse_shape
        print "pulse phase is:", pulse_phase, "PI radians"
        print "pulse area is:", pulse_area, "PI radians"
        print "pulse polarization is:", pulse_pol
        print "pulse coupling strength is:", pulse_dipole, "Debye\n"
        if pulse_shape == DICHROMATIC:
            print "pulse hole width is:", pulse_hole_width, "eV"

    # determine unit vector of pulse polarization
    if pulse_pol == POL_H:
        pulse_xcomp = complex(1.0, 0.0)
        pulse_ycomp = complex(0.0, 0.0)
    elif pulse_pol == POL_V:
        pulse_xcomp = complex(0.0, 0.0)
        pulse_ycomp = complex(1.0, 0.0)
    elif pulse_pol == POL_DIA:
        pulse_xcomp = (1.0/SQRT2)*complex(1.0, 0.0)
        pulse_ycomp = (1.0/SQRT2)*complex(1.0, 0.0)
    elif pulse_pol == POL_ADIA:
        pulse_xcomp = (1.0/SQRT2)*complex(1.0, 0.0)
        pulse_ycomp = (1.0/SQRT2)*complex(-1.0, 0.0)
    elif pulse_pol == POL_LCP:
        pulse_xcomp = (1.0/SQRT2)*complex(1.0, 0.0)
        pulse_ycomp = (1.0/SQRT2)*complex(0.0, 1.0)
    elif pulse_pol == POL_RCP:
        pulse_xcomp = (1.0/SQRT2)*complex(1.0, 0.0)
        pulse_ycomp = (1.0/SQRT2)*complex(0.0, -1.0)
    else:
        print "Invalid Pulse Polarization!"

    # convert from fraction of PI to radians
    pulse_phase = pulse_phase*PI
    pulse_area = pulse_area*PI

    # determine peak electric field for given pulse shape, pulse width, pulse area, and dipole moment
    if pulse_shape == GAUSSIAN:
        pulse_EO = H_BAR*pulse_area*pow(GAUSSIAN_CONST/PI, 0.5)/(pulse_dipole*DEBYE_TO_CM*pulse_width*1.0e-15);
    elif pulse_shape == SECH:
        pulse_EO = SECH_CONST*H_BAR*pulse_area/(pulse_dipole*DEBYE_TO_CM*PI*pulse_width*1.0e-15)

    elif pulse_shape == SQUARE:
        pulse_EO = H_BAR*pulse_area/(pulse_dipole*DEBYE_TO_CM*pulse_width*1.0e-15)
    elif pulse_shape == LORENTZIAN: 
        pulse_EO = H_BAR*pulse_area/(pulse_dipole*DEBYE_TO_CM*pulse_width*1.0e-15)
    elif pulse_shape == DICHROMATIC: 
        pulse_EO = H_BAR*pulse_area*pow(GAUSSIAN_CONST/PI, 0.5)/(pulse_dipole*DEBYE_TO_CM*pulse_width*1.0e-15)
    else:
        print "Invalid Pulse Shape!"


    # convert pulse pararmeters to Atomic Rydberg Units (ARU)
    pulse_omega_o = convert(pulse_omega_o, EV_TO_ARU)    # convert from eV to ARU (angular frequency)
    pulse_detun = convert(pulse_detun, EV_TO_ARU)        # convert from eV to ARU (angular frequency)
    pulse_width = convert(pulse_width, FEMTO_TO_ARU)    # convert pulse width from fs to ARU
    pulse_delay = convert(pulse_delay, FEMTO_TO_ARU)    # convert pulse width from fs to ARU
    pulse_EO = convert(pulse_EO, ELEC_TO_ARU)            # convert pulse electric field to ARU
    pulse_chirp = convert(pulse_chirp, FEMTO2_TO_ARU)    # convert pulse chirp from fs^2 to ARU
    pulse_hole_width = convert(pulse_hole_width, EV_TO_ARU) # convert from eV to ARU (spectral hole for dichromatic pulse)
    
    # write data to PulseParams data structure

    pulse_params = PulseParams(pulse_omega_o, pulse_detun, pulse_width, pulse_delay, pulse_chirp, pulse_shape, pulse_pol, pulse_phase, pulse_area, pulse_EO, pulse_dipole, pulse_xcomp, pulse_ycomp,pulse_hole_width)
    


    #------------------------------------------------------------#

    # read quantum dot parameters

    # write parameters to screen if requested
    if dump_to_screen == True:
        print "reading qdot parameters..."
    qdot_omega_yo = config['qdot']['omega_yo']
    qdot_omega_xyfine = config['qdot']['omega_xyfine']
    qdot_omega_bind = config['qdot']['omega_bind']
    qdot_d_xo = config['qdot']['d_xo']
    qdot_d_yo = config['qdot']['d_yo']
    qdot_d_xy = config['qdot']['d_xy']
    qdot_d_bo = config['qdot']['d_bo']
    qdot_d_bx = config['qdot']['d_bx']
    qdot_d_by = config['qdot']['d_by']
    qdot_T_oo = config['qdot']['T_oo']
    qdot_T_xx = config['qdot']['T_xx']
    qdot_T_yy = config['qdot']['T_yy']
    qdot_T_bb = config['qdot']['T_bb']
    qdot_T_xo = config['qdot']['T_xo']
    qdot_T_yo = config['qdot']['T_yo']
    qdot_T_bo = config['qdot']['T_bo']
    qdot_T_xy = config['qdot']['T_xy']
    qdot_T_bx = config['qdot']['T_bx']
    qdot_T_by = config['qdot']['T_by']
    qdot_eid_c = config['qdot']['eid_c']
    qdot_eid_b = config['qdot']['eid_b']

    initial_state_co_r = config['qdot']['initial_state_co_r']
    initial_state_co_phi = config['qdot']['initial_state_co_phi']
    initial_state_cx_r = config['qdot']['initial_state_cx_r']
    initial_state_cx_phi = config['qdot']['initial_state_cx_phi']
    initial_state_cy_r = config['qdot']['initial_state_cy_r']
    initial_state_cy_phi = config['qdot']['initial_state_cy_phi']
    initial_state_cb_r = config['qdot']['initial_state_cb_r']
    initial_state_cb_phi = config['qdot']['initial_state_cb_phi']

    desired_state_co_r = config['qdot']['desired_state_co_r']
    desired_state_co_phi = config['qdot']['desired_state_co_phi']
    desired_state_cx_r = config['qdot']['desired_state_cx_r']
    desired_state_cx_phi = config['qdot']['desired_state_cx_phi']
    desired_state_cy_r = config['qdot']['desired_state_cy_r']
    desired_state_cy_phi = config['qdot']['desired_state_cy_phi']
    desired_state_cb_r = config['qdot']['desired_state_cb_r']
    desired_state_cb_phi = config['qdot']['desired_state_cb_phi']

    # read in initial state and normalize
    qdot_initial_co = cmath.rect(initial_state_co_r, initial_state_co_phi*PI)
    qdot_initial_cx = cmath.rect(initial_state_cx_r, initial_state_cx_phi*PI)
    qdot_initial_cy = cmath.rect(initial_state_cy_r, initial_state_cy_phi*PI)
    qdot_initial_cb = cmath.rect(initial_state_cb_r, initial_state_cb_phi*PI)
    if config['run']['two_level']:
        qdot_initial_state = numpy.array([qdot_initial_co, qdot_initial_cy], dtype=complex)
    else:
        qdot_initial_state = numpy.array([qdot_initial_co, qdot_initial_cx, qdot_initial_cy, qdot_initial_cb], dtype=complex)
    qdot_initial_state = normalize_state(qdot_initial_state)


    # read in desired state and normalize
    qdot_desired_co = cmath.rect(desired_state_co_r, desired_state_co_phi*PI)
    qdot_desired_cx = cmath.rect(desired_state_cx_r, desired_state_cx_phi*PI)
    qdot_desired_cy = cmath.rect(desired_state_cy_r, desired_state_cy_phi*PI)
    qdot_desired_cb = cmath.rect(desired_state_cb_r, desired_state_cb_phi*PI)

    if config['run']['two_level']:
        qdot_desired_state = numpy.array([qdot_desired_co, qdot_desired_cy], dtype=complex)
    else:
        qdot_desired_state = numpy.array([qdot_desired_co, qdot_desired_cx, qdot_desired_cy, qdot_desired_cb], dtype=complex)
    qdot_desired_state = normalize_state(qdot_desired_state)

    if dump_to_screen == True:
        print "Initial State for QDOT1 is: {0}".format(qdot_initial_state)
        print "Desired State for QDOT1 is: {0}".format(qdot_desired_state)

    # convert quantum dot parameters to Atomic Rydberg Units (ARU)

    # convert quantum dot energies from eV to ARU (angular frequency)
    qdot_omega_yo = convert(qdot_omega_yo, EV_TO_ARU)
    qdot_omega_xyfine = convert(qdot_omega_xyfine, EV_TO_ARU)
    qdot_omega_bind = convert(qdot_omega_bind, EV_TO_ARU)

    # quantum dot dipole moments from Debye to ARU
    qdot_d_xo = convert(qdot_d_xo, DEBYE_TO_ARU)
    qdot_d_yo = convert(qdot_d_yo, DEBYE_TO_ARU)
    qdot_d_xy = convert(qdot_d_xy, DEBYE_TO_ARU)
    qdot_d_bo = convert(qdot_d_bo, DEBYE_TO_ARU)
    qdot_d_bx = convert(qdot_d_bx, DEBYE_TO_ARU)
    qdot_d_by = convert(qdot_d_by, DEBYE_TO_ARU)

    # quantum dot decay constants from ps to ARU
    qdot_T_xx = convert(qdot_T_xx, PICO_TO_ARU)
    qdot_T_yy = convert(qdot_T_yy, PICO_TO_ARU)
    qdot_T_bb = convert(qdot_T_bb, PICO_TO_ARU)
    qdot_T_xo = convert(qdot_T_xo, PICO_TO_ARU)
    qdot_T_yo = convert(qdot_T_yo, PICO_TO_ARU)
    qdot_T_bo = convert(qdot_T_bo, PICO_TO_ARU)
    qdot_T_xy = convert(qdot_T_xy, PICO_TO_ARU)
    qdot_T_bx = convert(qdot_T_bx, PICO_TO_ARU)
    qdot_T_by = convert(qdot_T_by, PICO_TO_ARU)

    qdot_eid_c = convert(qdot_eid_c, PICO_TO_ARU)



    # write data to QDotParams data structure
    qdot_params = QDotParams(qdot_omega_yo, qdot_omega_xyfine, qdot_omega_bind, qdot_d_xo, qdot_d_yo, qdot_d_xy, qdot_d_bo, qdot_d_bx, qdot_d_by, qdot_T_oo, qdot_T_xx, qdot_T_yy, qdot_T_bb, qdot_T_xo, qdot_T_yo, qdot_T_bo, qdot_T_xy, qdot_T_bx, qdot_T_by, qdot_eid_c, qdot_eid_b, qdot_initial_state, qdot_desired_state)

    # read and write data for the second quantum dot
    qdot2_omega_yo = config['qdot2']['omega_yo']
    qdot2_omega_xyfine = config['qdot2']['omega_xyfine']
    qdot2_omega_bind = config['qdot2']['omega_bind']
    qdot2_d_xo = config['qdot2']['d_xo']
    qdot2_d_yo = config['qdot2']['d_yo']
    qdot2_d_xy = config['qdot2']['d_xy']
    qdot2_d_bo = config['qdot2']['d_bo']
    qdot2_d_bx = config['qdot2']['d_bx']
    qdot2_d_by = config['qdot2']['d_by']
    qdot2_T_oo = config['qdot2']['T_oo']
    qdot2_T_xx = config['qdot2']['T_xx']
    qdot2_T_yy = config['qdot2']['T_yy']
    qdot2_T_bb = config['qdot2']['T_bb']
    qdot2_T_xo = config['qdot2']['T_xo']
    qdot2_T_yo = config['qdot2']['T_yo']
    qdot2_T_bo = config['qdot2']['T_bo']
    qdot2_T_xy = config['qdot2']['T_xy']
    qdot2_T_bx = config['qdot2']['T_bx']
    qdot2_T_by = config['qdot2']['T_by']
    qdot2_eid_c = config['qdot2']['eid_c']
    qdot2_eid_b = config['qdot2']['eid_b']

    initial_state_co_r = config['qdot2']['initial_state_co_r']
    initial_state_co_phi = config['qdot2']['initial_state_co_phi']
    initial_state_cx_r = config['qdot2']['initial_state_cx_r']
    initial_state_cx_phi = config['qdot2']['initial_state_cx_phi']
    initial_state_cy_r = config['qdot2']['initial_state_cy_r']
    initial_state_cy_phi = config['qdot2']['initial_state_cy_phi']
    initial_state_cb_r = config['qdot2']['initial_state_cb_r']
    initial_state_cb_phi = config['qdot2']['initial_state_cb_phi']

    desired_state_co_r = config['qdot2']['desired_state_co_r']
    desired_state_co_phi = config['qdot2']['desired_state_co_phi']
    desired_state_cx_r = config['qdot2']['desired_state_cx_r']
    desired_state_cx_phi = config['qdot2']['desired_state_cx_phi']
    desired_state_cy_r = config['qdot2']['desired_state_cy_r']
    desired_state_cy_phi = config['qdot2']['desired_state_cy_phi']
    desired_state_cb_r = config['qdot2']['desired_state_cb_r']
    desired_state_cb_phi = config['qdot2']['desired_state_cb_phi']

    # read in initial state and normalize it if necessary
    qdot2_initial_co = cmath.rect(initial_state_co_r, initial_state_co_phi*PI)
    qdot2_initial_cx = cmath.rect(initial_state_cx_r, initial_state_cx_phi*PI)
    qdot2_initial_cy = cmath.rect(initial_state_cy_r, initial_state_cy_phi*PI)
    qdot2_initial_cb = cmath.rect(initial_state_cb_r, initial_state_cb_phi*PI)

    if config['run']['two_level']:
        qdot2_initial_state = numpy.array([qdot2_initial_co, qdot2_initial_cy], dtype=complex)
    else:
        qdot2_initial_state = numpy.array([qdot2_initial_co, qdot2_initial_cx, qdot2_initial_cy, qdot2_initial_cb], dtype=complex)
    qdot2_initial_state = normalize_state(qdot2_initial_state)


    # read in desired state and normalize it if necessary
    qdot2_desired_co = cmath.rect(desired_state_co_r, desired_state_co_phi*PI)
    qdot2_desired_cx = cmath.rect(desired_state_cx_r, desired_state_cx_phi*PI)
    qdot2_desired_cy = cmath.rect(desired_state_cy_r, desired_state_cy_phi*PI)
    qdot2_desired_cb = cmath.rect(desired_state_cb_r, desired_state_cb_phi*PI)

    if config['run']['two_level']:
        qdot2_desired_state = numpy.array([qdot2_desired_co, qdot2_desired_cy], dtype=complex)
    else:
        qdot2_desired_state = numpy.array([qdot2_desired_co, qdot2_desired_cx, qdot2_desired_cy, qdot2_desired_cb], dtype=complex)
    qdot2_desired_state = normalize_state(qdot2_desired_state)
    
    if dump_to_screen == True:
        print "Initial State for QDOT2 is: {0}".format(qdot2_initial_state)
        print "Desired State for QDOT2 is: {0}".format(qdot2_desired_state)

    # convert quantum dot parameters to Atomic Rydberg Units (ARU)

    # convert quantum dot energies from eV to ARU (angular frequency)
    qdot2_omega_yo = convert(qdot2_omega_yo, EV_TO_ARU)
    qdot2_omega_xyfine = convert(qdot2_omega_xyfine, EV_TO_ARU)
    qdot2_omega_bind = convert(qdot2_omega_bind, EV_TO_ARU)

    # quantum dot dipole moments from Debye to ARU
    qdot2_d_xo = convert(qdot2_d_xo, DEBYE_TO_ARU)
    qdot2_d_yo = convert(qdot2_d_yo, DEBYE_TO_ARU)
    qdot2_d_xy = convert(qdot2_d_xy, DEBYE_TO_ARU)
    qdot2_d_bo = convert(qdot2_d_bo, DEBYE_TO_ARU)
    qdot2_d_bx = convert(qdot2_d_bx, DEBYE_TO_ARU)
    qdot2_d_by = convert(qdot2_d_by, DEBYE_TO_ARU)

    # quantum dot decay constants from ps to ARU
    qdot2_T_xx = convert(qdot2_T_xx, PICO_TO_ARU)
    qdot2_T_yy = convert(qdot2_T_yy, PICO_TO_ARU)
    qdot2_T_bb = convert(qdot2_T_bb, PICO_TO_ARU)
    qdot2_T_xo = convert(qdot2_T_xo, PICO_TO_ARU)
    qdot2_T_yo = convert(qdot2_T_yo, PICO_TO_ARU)
    qdot2_T_bo = convert(qdot2_T_bo, PICO_TO_ARU)
    qdot2_T_xy = convert(qdot2_T_xy, PICO_TO_ARU)
    qdot2_T_bx = convert(qdot2_T_bx, PICO_TO_ARU)
    qdot2_T_by = convert(qdot2_T_by, PICO_TO_ARU)

    qdot2_eid_c = convert(qdot2_eid_c, PICO_TO_ARU)

    qdot2_params = QDotParams(qdot2_omega_yo, qdot2_omega_xyfine, qdot2_omega_bind, qdot2_d_xo, qdot2_d_yo, qdot2_d_xy, qdot2_d_bo, qdot2_d_bx, qdot2_d_by, qdot2_T_oo, qdot2_T_xx, qdot2_T_yy, qdot2_T_bb, qdot2_T_xo, qdot2_T_yo, qdot2_T_bo, qdot2_T_xy, qdot2_T_bx, qdot2_T_by, qdot2_eid_c, qdot2_eid_b, qdot2_initial_state, qdot2_desired_state)


    # read and write data for the second quantum dot
    qdot3_omega_yo = config['qdot3']['omega_yo']
    qdot3_omega_xyfine = config['qdot3']['omega_xyfine']
    qdot3_omega_bind = config['qdot3']['omega_bind']
    qdot3_d_xo = config['qdot3']['d_xo']
    qdot3_d_yo = config['qdot3']['d_yo']
    qdot3_d_xy = config['qdot3']['d_xy']
    qdot3_d_bo = config['qdot3']['d_bo']
    qdot3_d_bx = config['qdot3']['d_bx']
    qdot3_d_by = config['qdot3']['d_by']
    qdot3_T_oo = config['qdot3']['T_oo']
    qdot3_T_xx = config['qdot3']['T_xx']
    qdot3_T_yy = config['qdot3']['T_yy']
    qdot3_T_bb = config['qdot3']['T_bb']
    qdot3_T_xo = config['qdot3']['T_xo']
    qdot3_T_yo = config['qdot3']['T_yo']
    qdot3_T_bo = config['qdot3']['T_bo']
    qdot3_T_xy = config['qdot3']['T_xy']
    qdot3_T_bx = config['qdot3']['T_bx']
    qdot3_T_by = config['qdot3']['T_by']
    qdot3_eid_c = config['qdot3']['eid_c']
    qdot3_eid_b = config['qdot3']['eid_b']

    initial_state_co_r = config['qdot3']['initial_state_co_r']
    initial_state_co_phi = config['qdot3']['initial_state_co_phi']
    initial_state_cx_r = config['qdot3']['initial_state_cx_r']
    initial_state_cx_phi = config['qdot3']['initial_state_cx_phi']
    initial_state_cy_r = config['qdot3']['initial_state_cy_r']
    initial_state_cy_phi = config['qdot3']['initial_state_cy_phi']
    initial_state_cb_r = config['qdot3']['initial_state_cb_r']
    initial_state_cb_phi = config['qdot3']['initial_state_cb_phi']

    desired_state_co_r = config['qdot3']['desired_state_co_r']
    desired_state_co_phi = config['qdot3']['desired_state_co_phi']
    desired_state_cx_r = config['qdot3']['desired_state_cx_r']
    desired_state_cx_phi = config['qdot3']['desired_state_cx_phi']
    desired_state_cy_r = config['qdot3']['desired_state_cy_r']
    desired_state_cy_phi = config['qdot3']['desired_state_cy_phi']
    desired_state_cb_r = config['qdot3']['desired_state_cb_r']
    desired_state_cb_phi = config['qdot3']['desired_state_cb_phi']

    # read in initial state and normalize it if necessary
    qdot3_initial_co = cmath.rect(initial_state_co_r, initial_state_co_phi*PI)
    qdot3_initial_cx = cmath.rect(initial_state_cx_r, initial_state_cx_phi*PI)
    qdot3_initial_cy = cmath.rect(initial_state_cy_r, initial_state_cy_phi*PI)
    qdot3_initial_cb = cmath.rect(initial_state_cb_r, initial_state_cb_phi*PI)

    if config['run']['two_level']:
        qdot3_initial_state = numpy.array([qdot3_initial_co, qdot3_initial_cy], dtype=complex)
    else:
        qdot3_initial_state = numpy.array([qdot3_initial_co, qdot3_initial_cx, qdot3_initial_cy, qdot3_initial_cb], dtype=complex)
    qdot3_initial_state = normalize_state(qdot3_initial_state)


    # read in desired state and normalize it if necessary
    qdot3_desired_co = cmath.rect(desired_state_co_r, desired_state_co_phi*PI)
    qdot3_desired_cx = cmath.rect(desired_state_cx_r, desired_state_cx_phi*PI)
    qdot3_desired_cy = cmath.rect(desired_state_cy_r, desired_state_cy_phi*PI)
    qdot3_desired_cb = cmath.rect(desired_state_cb_r, desired_state_cb_phi*PI)

    if config['run']['two_level']:
        qdot3_desired_state = numpy.array([qdot3_desired_co, qdot3_desired_cy], dtype=complex)
    else:
        qdot3_desired_state = numpy.array([qdot3_desired_co, qdot3_desired_cx, qdot3_desired_cy, qdot3_desired_cb], dtype=complex)
    qdot3_desired_state = normalize_state(qdot3_desired_state)

    if dump_to_screen == True:
        print "Initial State for QDOT3 is: {0}".format(qdot3_initial_state)
        print "Desired State for QDOT3 is: {0}".format(qdot3_desired_state)

    # convert quantum dot parameters to Atomic Rydberg Units (ARU)

    # convert quantum dot energies from eV to ARU (angular frequency)
    qdot3_omega_yo = convert(qdot3_omega_yo, EV_TO_ARU)
    qdot3_omega_xyfine = convert(qdot3_omega_xyfine, EV_TO_ARU)
    qdot3_omega_bind = convert(qdot3_omega_bind, EV_TO_ARU)

    # quantum dot dipole moments from Debye to ARU
    qdot3_d_xo = convert(qdot3_d_xo, DEBYE_TO_ARU)
    qdot3_d_yo = convert(qdot3_d_yo, DEBYE_TO_ARU)
    qdot3_d_xy = convert(qdot3_d_xy, DEBYE_TO_ARU)
    qdot3_d_bo = convert(qdot3_d_bo, DEBYE_TO_ARU)
    qdot3_d_bx = convert(qdot3_d_bx, DEBYE_TO_ARU)
    qdot3_d_by = convert(qdot3_d_by, DEBYE_TO_ARU)

    # quantum dot decay constants from ps to ARU
    qdot3_T_xx = convert(qdot3_T_xx, PICO_TO_ARU)
    qdot3_T_yy = convert(qdot3_T_yy, PICO_TO_ARU)
    qdot3_T_bb = convert(qdot3_T_bb, PICO_TO_ARU)
    qdot3_T_xo = convert(qdot3_T_xo, PICO_TO_ARU)
    qdot3_T_yo = convert(qdot3_T_yo, PICO_TO_ARU)
    qdot3_T_bo = convert(qdot3_T_bo, PICO_TO_ARU)
    qdot3_T_xy = convert(qdot3_T_xy, PICO_TO_ARU)
    qdot3_T_bx = convert(qdot3_T_bx, PICO_TO_ARU)
    qdot3_T_by = convert(qdot3_T_by, PICO_TO_ARU)

    qdot3_eid_c = convert(qdot3_eid_c, PICO_TO_ARU)

    qdot3_params = QDotParams(qdot3_omega_yo, qdot3_omega_xyfine, qdot3_omega_bind, qdot3_d_xo, qdot3_d_yo, qdot3_d_xy, qdot3_d_bo, qdot3_d_bx, qdot3_d_by, qdot3_T_oo, qdot3_T_xx, qdot3_T_yy, qdot3_T_bb, qdot3_T_xo, qdot3_T_yo, qdot3_T_bo, qdot3_T_xy, qdot3_T_bx, qdot3_T_by, qdot3_eid_c, qdot3_eid_b, qdot3_initial_state, qdot3_desired_state)
    #------------------------------------------------------------#


    # read la-phonon eid parameters
    omega_c_eid = config['laphononeid']['omega_c_eid']
    alpha_eid = config['laphononeid']['alpha_eid']
    Omega_start_eid = config['laphononeid']['Omega_start_eid']
    Omega_end_eid = config['laphononeid']['Omega_end_eid']
    Omega_step_eid = config['laphononeid']['Omega_step_eid']
    omega_end_eid = config['laphononeid']['omega_end_eid']
    T_eid = config['laphononeid']['T_eid']
    read_kernel = config['laphononeid']['read_kernel']

    # write parameters to screen if requested
    if dump_to_screen and run_phonon_eid:
        print "LA phonon EID parameters..."
        print "omega_c_eid:", omega_c_eid, "eV"
        print "alpha_eid:", alpha_eid, "ps^2"
        print "Omega_start_eid: ", Omega_start_eid, "eV"
        print "Omega_end_eid: ", Omega_end_eid, "eV"
        print "Omega_step_eid: ", Omega_step_eid, "eV"
        print "T_eid: ", T_eid, "K"

    # convert eid pararmeters to Atomic Rydberg Units (ARU)
    omega_c_eid = convert(omega_c_eid, EV_TO_ARU)            # convert from eV to ARU
    alpha_eid = convert(alpha_eid, PICO2_TO_ARU)                        # convert from ps^2 to ARU
    Omega_start_eid = convert(Omega_start_eid, EV_TO_ARU)                # convert from fs to ARU
    Omega_end_eid = convert(Omega_end_eid, EV_TO_ARU)                # convert from fs to ARU
    Omega_step_eid = convert(Omega_step_eid, EV_TO_ARU)                # convert from fs to ARU

    # write data to EIDParams data structure
    laphononeid_params = LAPhononEIDParams(omega_c_eid, alpha_eid, Omega_start_eid, Omega_end_eid, Omega_step_eid, T_eid, read_kernel, None , None)


    #------------------------------------------------------------#

    # read non-markovian eid parameters
    kappa_0 = config['nonmarkovian']['kappa_0']
    kappa_1 = config['nonmarkovian']['kappa_1']
    kappa_2 = config['nonmarkovian']['kappa_2']

    # write parameters to screen if requested
    if dump_to_screen and run_phonon_eid:
        print "Non-Markovian reservoir EID parameters..."
        print "kappa_0:", kappa_0, "ps-1"
        print "kappa_1:", kappa_1, "ps-1"
        print "kappa_2: ", kappa_2, "ps-1"


    # convert eid pararmeters to Atomic Rydberg Units (ARU)
    kappa_0 = convert(kappa_0, INV_PICO_TO_ARU)                        # convert from ps^-1 to ARU
    kappa_1 = convert(kappa_1, INV_PICO_TO_ARU)                        # convert from ps^-1 to ARU
    kappa_2 = convert(kappa_2, INV_PICO_TO_ARU)                        # convert from ps^-1 to ARU


    # write data to EIDParams data structure
    nonmarkovian_params = NonMarkovianParams(kappa_0, kappa_1, kappa_2)


    #------------------------------------------------------------#

    return {'run_params': run_params, 'rkdp_params': rkdp_params, 'mask_params': mask_params, 'ampmask_params': ampmask_params, 'phasemask_params': phasemask_params, 'pulse_params': pulse_params, 'qdot_params': qdot_params, 'qdot2_params': qdot2_params, 'qdot3_params': qdot3_params, 'laphnoneid_params': laphononeid_params, 'nonmarkovian_params': nonmarkovian_params}

def write_params (x, config_file):
    config = ConfigObj(config_file, configspec='configspec.ini')
    validator = Validator()
    result = config.validate(validator)

    # determine mask that is being used
    run_mask = config['run']['mask']

    # read in list of parameters that are being optimized
    param_list = config['masks'][run_mask]['param_list']
    NOPTS = len(param_list)

    # write each parameters to file
    for i in range(NOPTS):
        config['masks'][run_mask][param_list[i]] = x[i]

    config.write()

def read_params (config_file):
    config = ConfigObj(config_file, configspec='configspec.ini')
    validator = Validator()
    result = config.validate(validator)

    # determine mask that is being used
    run_mask = config['run']['mask']

    # read in list of parameters that are being optimized
    param_list = config['masks'][run_mask]['param_list']
    NOPTS = len(param_list)
    x = numpy.zeros(NOPTS)
    # write each parameters to file
    for i in range(NOPTS):
        x[i] = config['masks'][run_mask][param_list[i]]

    return x


# convert parameters from/to ARU units
def convert_aru(x, config_file, conversion):
    config = ConfigObj(config_file, configspec='configspec.ini')
    validator = Validator()
    result = config.validate(validator)

    # determine mask that is being used
    run_mask = config['run']['mask']
    param_list = config['masks'][run_mask]['param_list']
    NOPTS = len(param_list)

    x_convert = []
    # convert parameters according to conversion factor
    for i in range(NOPTS):
        # read in conversion factor
        param_ARU = param_list[i] + '_ARU'
        if (conversion == TO_ARU):
            x_convert.append(config['masks'][run_mask][param_ARU][0])
        elif (conversion == FROM_ARU):
            x_convert.append(config['masks'][run_mask][param_ARU][1])
        else:
            "do nothing"
        # if conversion is from AREA (in fraction of PI) to ARU, pass extra kwargs
        if (x_convert[i] == AREA_TO_ARU or x_convert[i] == ARU_TO_AREA):
            dipole = config['pulse']['dipole']
            width = config['pulse']['width']
            shape = config['pulse']['shape']
            kwargs = {'dipole_moment': dipole, 'pulse_width': width, 'pulse_shape': shape}
        else:
            kwargs = {}
        x[i] = convert(x[i], x_convert[i], **kwargs)

    return x



