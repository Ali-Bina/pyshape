#! /home/ajan/anaconda/bin/python
import scipy.constants
import numpy

PI = scipy.constants.pi
TWOPI  = 2*PI
SQRT2  = numpy.sqrt(2)

# physical constants in mksa units
E_CHARGE = scipy.constants.e        # electron charge (Coulomb)
H_BAR = scipy.constants.hbar        # hbar (Joule.second)
C = scipy.constants.c               # speed of light (m/s)
EPSO = scipy.constants.epsilon_0    # permittivity of vacuum (F/m)
KB = scipy.constants.k              # boltzmann constant (J/K)

DEBYE_TO_CM = 1E-21/C

# physical constants in ARU units
EPSO_ARU = 1/(4*PI)
H_BAR_ARU = 1.0
E_CHARGE_ARU = 1.0
KB_ARU = scipy.constants.k/scipy.constants.physical_constants["atomic unit of energy"][0]
C_ARU = scipy.constants.c/scipy.constants.physical_constants["atomic unit of velocity"][0]

# conversion factors to convert mksa units to atomic rydberg units (aru)
ARU_TIME = scipy.constants.physical_constants["atomic unit of time"][0]                       # (seconds)
ARU_FREQUENCY = 1.0/scipy.constants.physical_constants["atomic unit of time"][0]              # (Hertz)
ARU_ELECTRIC_FIELD = scipy.constants.physical_constants["atomic unit of electric field"][0]   # (volts/meter)
ARU_DEBYE = scipy.constants.physical_constants["atomic unit of electric dipole mom."][0]/DEBYE_TO_CM      # (Coulomb.meter)

PICO = scipy.constants.pico                      # conversion from pico to 1e0
FEMTO = scipy.constants.femto                     # conversion from femto to 1e0

# pulse constants
GAUSSIAN_CONST = 2*numpy.log(2)
SECH_CONST = numpy.log(3 + 2*numpy.sqrt(2))


STATE_XX = 0
STATE_YY = 1
STATE_BB = 2
STATE_OO = 3

DOT1 = "DOT1"
DOT2 = "DOT2"
DOT3 = "DOT3"


# string constants used to convert between units (convert.h)
FEMTO_TO_ARU = "FEMTO_TO_ARU"
ARU_TO_FEMTO = "ARU_TO_FEMTO"
FEMTO2_TO_ARU = "FEMTO2_TO_ARU"
ARU_TO_FEMTO2 = "ARU_TO_FEMTO2"
FEMTO3_TO_ARU = "FEMTO3_TO_ARU"
ARU_TO_FEMTO3 = "ARU_TO_FEMTO3"
FEMTO4_TO_ARU = "FEMTO4_TO_ARU"
ARU_TO_FEMTO4 = "ARU_TO_FEMTO4"
PICO_TO_ARU = "PICO_TO_ARU"
ARU_TO_PICO = "ARU_TO_PICO"
PICO2_TO_ARU = "PICO2_TO_ARU"
INV_PICO_TO_ARU = "INV_PICO_TO_ARU"
ARU_TO_INV_PICO = "ARU_TO_INV_PICO"
ELEC_TO_ARU = "ELEC_TO_ARU"
ARU_TO_ELEC = "ARU_TO_ELEC"
DEBYE_TO_ARU = "DEBYE_TO_ARU"
ARU_TO_DEBYE = "ARU_TO_DEBYE"
EV_TO_ARU = "EV_TO_ARU"
ARU_TO_EV = "ARU_TO_EV"
AREA_TO_ELEC = "AREA_TO_ELEC"
AREA_TO_ELECARU = "AREA_TO_ELECARU"
ARU_TO_FBWIDTH = "ARU_TO_FBWIDTH"
EO_TO_AREA = "EO_TO_AREA"
FRAC_TO_ANGLE = "FRAC_TO_ANGLE"
ANGLE_TO_FRAC = "ANGLE_TO_FRAC"
AREA_TO_ARU = "AREA_TO_ARU"
ARU_TO_AREA = "ARU_TO_AREA"
NONE = "NONE"

TO_ARU = "TO_ARU"
FROM_ARU ="FROM_ARU"

GAUSSIAN = "GAUSSIAN"
SECH = "SECH"
SQUARE="SQUARE"
LORENTZIAN="LORENTZIAN"
DICHROMATIC="DICHROMATIC"

PUMP = "PUMP"
PROBE = "PROBE"

AMPMASK_CHEN = "AMPMASK_CHEN"
PHASEMASK_COS = "PHASEMASK_COS"
PHASEMASK_CHIRP = "PHASEMASK_CHIRP"
UNSHAPED = "UNSHAPED"

ARP = "ARP"
CROT = "CROT"
TLS = "TLS"
THREETLS = "3TLS"

OPT = "OPT"
TDEP = "TDEP"
EXP = "EXP"
DTRANS = "DTRANS"

NONLINEAR = "NONLINEAR"
GENETIC = "GENETIC"

POL_H = "POL_H"
POL_V = "POL_V"
POL_DIA = "POL_DIA"
POL_ADIA = "POL_ADIA"
POL_LCP = "POL_LCP"
POL_RCP = "POL_RCP"
POL_LIN = "POL_LIN"

