#! /home/ajan/anaconda/bin/python
from constants import *

def convert(param, unit, **kwargs):
    if (unit == FEMTO_TO_ARU):             # convert from fs to ARU
        return (param * FEMTO) / ARU_TIME
    elif (unit == ARU_TO_FEMTO):
        return (param * ARU_TIME) / FEMTO
    elif (unit == FEMTO2_TO_ARU):       # convert from fs^2 to ARU
        return (param * pow((FEMTO / ARU_TIME), 2))
    elif (unit == ARU_TO_FEMTO2):
        return (param * pow((ARU_TIME / FEMTO), 2))
    elif (unit == FEMTO3_TO_ARU):       # convert from fs^3 to ARU
        return (param * pow((FEMTO / ARU_TIME), 3))
    elif (unit == ARU_TO_FEMTO3):
        return (param * pow((ARU_TIME / FEMTO), 3))
    elif (unit == FEMTO4_TO_ARU):       # convert from fs^4 to ARU
        return (param * pow((FEMTO / ARU_TIME), 4))
    elif (unit == ARU_TO_FEMTO4):
        return (param * pow((ARU_TIME / FEMTO), 4))
    elif (unit == PICO_TO_ARU):   # convert from ps to ARU
        return (param * PICO) / ARU_TIME
    elif (unit == ARU_TO_PICO):
        return (param * ARU_TIME) / PICO
    elif (unit == PICO2_TO_ARU):       # convert from ps^2 to ARU
        return (param * pow((PICO / ARU_TIME), 2))
    elif (unit == INV_PICO_TO_ARU):    # convert from ps^-1 to ARU
        return (param / PICO) * ARU_TIME
    elif (unit == ARU_TO_INV_PICO):
        return (param / ARU_TIME) * PICO
    elif (unit == ELEC_TO_ARU):   # convert from V/m to ARU
        return param / ARU_ELECTRIC_FIELD
    elif (unit == ARU_TO_ELEC):
        return param * ARU_ELECTRIC_FIELD
    elif (unit == DEBYE_TO_ARU):  # convert from debye to ARU
        return param / ARU_DEBYE
    elif (unit == ARU_TO_DEBYE):
        return param * ARU_DEBYE
    elif (unit == EV_TO_ARU):     # convert from eV to ARU (angular frequency)
        return (param * E_CHARGE / H_BAR) / ARU_FREQUENCY
    elif (unit == ARU_TO_EV):
        return (param * ARU_FREQUENCY * H_BAR) / E_CHARGE
    elif (unit == ARU_TO_FBWIDTH):     # convert pulse width from femto to freq bandwidth in meV (hbar*DeltaOmega)
        return 4.0e18 * H_BAR * log(2.0)/(E_CHARGE * param * ARU_TIME / FEMTO)
    elif (unit ==  FRAC_TO_ANGLE):
        return param*PI
    elif (unit == ANGLE_TO_FRAC):
        return param/PI
    elif (unit == AREA_TO_ARU):            # convert pulse area from fraction of PI to ELEC (V/m) to ARU
        d = kwargs['dipole_moment']
        pulse_width = kwargs['pulse_width']
        shape = kwargs['pulse_shape']
        ampmask=kwargs['ampmask']
        pulse_hole_width_mask=kwargs['holeWidthMask']
        if (shape == GAUSSIAN):
            if ampmask == 'dichrome':
                
                def iFTA2(t, delta, E0, tau):
                    import mpmath as mp
                    cg = delta
                    cg1 = E0
                    cg5 = tau
                    cg7 = t
                    return -cg1 * mp.exp(-32 * mp.log(2) ** 3 * cg7 ** 2 / (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2) / cg5 ** 2) * (cg * 4 ** (-cg7 ** 2 * (cg ** 2 * cg5 ** 2 - 16 * mp.log(2) ** 2) / cg5 ** 2 / (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2)) * cg5 - 4 ** (-cg7 ** 2 * cg ** 2 / (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2)) * mp.sqrt(cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2)) * (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2) ** (-0.1e1 / 0.2e1)
                def integrate(delta, EO, Tau):
                    from scipy.integrate import quad
                    import mpmath as mp
                    import numpy as np

                    g=lambda t: np.abs(iFTA2(t, delta, EO, Tau))
                    return quad(g, -np.inf, np.inf)[0]
                pulse_hole_width_mask  *= 1.602e-19 / H_BAR / (2 * PI) *1e15 #hole width in fs
                A = integrate(pulse_hole_width_mask, 1, pulse_width) / 1e15
                print "I ran (to EO)", A
                return H_BAR*(param*PI)/(d*DEBYE_TO_CM)/A

            return (param*PI)*(H_BAR*pow(GAUSSIAN_CONST/PI, 0.5)/(d*DEBYE_TO_CM*pulse_width*1.0e-15)) / ARU_ELECTRIC_FIELD
        elif (shape == SECH):
            return (param*PI)*(H_BAR*SECH_CONST/(d*DEBYE_TO_CM*PI*pulse_width*1.0e-15)) / ARU_ELECTRIC_FIELD
        elif (shape == SQUARE):
            return (param*PI)*(H_BAR/(d*DEBYE_TO_CM*pulse_width*1.0e-15)) / ARU_ELECTRIC_FIELD
        elif (shape == LORENTZIAN):
            return (param*PI)*(H_BAR/(d*DEBYE_TO_CM*PI*1.0e-15))/ARU_ELECTRIC_FIELD
        elif (shape == DICHROMATIC):
            return (param*PI)*(H_BAR*pow(GAUSSIAN_CONST/PI, 0.5)/(d*DEBYE_TO_CM*pulse_width*1.0e-15)) / ARU_ELECTRIC_FIELD
    elif (unit == ARU_TO_AREA):            # convert from ARU to ELEC (V/m) to pulse area in fraction of PI
        d = kwargs['dipole_moment']                # dipole moment in Debye
        pulse_width = kwargs['pulse_width']              # pulse width in femtoseconds
        shape = kwargs['pulse_shape']
        ampmask=kwargs['ampmask']
        pulse_hole_width_mask=kwargs['holeWidthMask']
        if (shape == GAUSSIAN):
            if ampmask == 'dichrome':

                def iFTA2(t, delta, E0, tau):
                    import mpmath as mp
                    cg = delta
                    cg1 = E0
                    cg5 = tau
                    cg7 = t
                    return -cg1 * mp.exp(-32 * mp.log(2) ** 3 * cg7 ** 2 / (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2) / cg5 ** 2) * (cg * 4 ** (-cg7 ** 2 * (cg ** 2 * cg5 ** 2 - 16 * mp.log(2) ** 2) / cg5 ** 2 / (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2)) * cg5 - 4 ** (-cg7 ** 2 * cg ** 2 / (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2)) * mp.sqrt(cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2)) * (cg ** 2 * cg5 ** 2 + 16 * mp.log(2) ** 2) ** (-0.1e1 / 0.2e1)
                def integrate(delta, EO, Tau):
                    from scipy.integrate import quad
                    import mpmath as mp
                    import numpy as np

                    g=lambda t: np.abs(iFTA2(t, delta, EO, Tau))
                    return quad(g, -np.inf, np.inf)[0]
                pulse_hole_width_mask  *= 1.602e-19 / H_BAR / (2 * PI) *1e15 #hole width in fs
                A = integrate(pulse_hole_width_mask, 1, pulse_width) / 1e15
                return 1/H_BAR*(param*ARU_ELECTRIC_FIELD)*(d*DEBYE_TO_CM)*A/PI

            return param * ARU_ELECTRIC_FIELD * (1/PI) * (d*DEBYE_TO_CM*pulse_width*1.0e-15)/(H_BAR*pow(GAUSSIAN_CONST/PI, 0.5))
        elif (shape == SECH):
            return param * ARU_ELECTRIC_FIELD * (1/PI)*(d*DEBYE_TO_CM*PI*pulse_width*1.0e-15)/(SECH_CONST*H_BAR)
        elif (shape == SQUARE):
            return param * ARU_ELECTRIC_FIELD * (1/PI)*(d*DEBYE_TO_CM*pulse_width*1.0e-15)/(H_BAR)
        elif (shape == LORENTZIAN):
            return param * ARU_ELECTRIC_FIELD * (1/PI)*(d*DEBYE_TO_CM*pulse_width*1.0e-15)/(H_BAR)
        elif (shape == DICHROMATIC):
            return param * ARU_ELECTRIC_FIELD * (1/PI) * (d*DEBYE_TO_CM*pulse_width*1.0e-15)/(H_BAR*pow(GAUSSIAN_CONST/PI, 0.5))
            
        elif (shape != SECH) and (shape != GAUSSIAN) and (shape != SQUARE) and (shape != LORENTZIAN) and (shape != DICHROMATIC):
	    print "Invalid Pulse shape"
    elif (unit == NONE):
        return param
    else:
        print unit, "is an invalid unit conversion"


