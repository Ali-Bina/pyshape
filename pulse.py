#!/usr/bin/python
# Filename: pulse.py
# version before fixing rabi freq

import matplotlib.pyplot as plt
import numpy
import scipy.fftpack
from scipy import integrate
from constants import *
from convert import *
from phasemasks import *
from ampmasks import *
import sys

class Pulse(object):
    '''Represents a pulse'''

    def __init__(self, params):
        self.t_start = params['run_params'].run_t_start
        self.t_end = params['run_params'].run_t_end
        self.t_step = params['run_params'].run_t_step    # desired time step
        self.array_len = int((convert(self.t_end, ARU_TO_FEMTO) - convert(self.t_start, ARU_TO_FEMTO))/convert(self.t_step, ARU_TO_FEMTO))

        self.shape = params['pulse_params'].pulse_shape
        self.chirp = params['pulse_params'].pulse_chirp
        self.tau = params['pulse_params'].pulse_width
        self.t_o = params['pulse_params'].pulse_delay
        self.E_o = params['pulse_params'].pulse_EO
        self.omega_o = params['pulse_params'].pulse_omega_o
        self.dipole = params['pulse_params'].pulse_dipole
        self.phase = params['pulse_params'].pulse_phase
        self.hole_width = params['pulse_params'].pulse_hole_width
        self.t = numpy.linspace(self.t_start, self.t_end, self.array_len, endpoint=False)        # time array

        self.omega = 2.0*PI*scipy.fftpack.fftfreq(self.array_len, d=self.t_step)                                # (angular) frequency array
        self.efield = numpy.zeros((self.array_len), 'complex')    # time domain electric field
        self.Efield = numpy.zeros((self.array_len), 'complex')    # frequency domain electric field
        self.intensity = numpy.zeros(self.array_len)  # temporal intensity of pulse
        self.inst_freq = numpy.zeros(self.array_len)  # temporal intensity of pulse
        self.rabi_freq = numpy.zeros(self.array_len)  # temporal intensity of pulse
        self.actual_rabi_freq = numpy.zeros(self.array_len) # for writing rabi_freq.txt 
        self.E_freq = numpy.zeros((self.array_len), 'complex')    # frequency domain electric field
        self.efield_envelope = numpy.zeros(self.array_len)   # temporal envelope of electric field abs(efield)
        self.intAC = numpy.zeros(self.array_len)  # intensity autocorrelation

        # create amplitude mask array - must default to one
        self.amp_mask = numpy.ones(self.array_len)

        # create phase mask array - must default to zero
        self.phase_mask = numpy.zeros(self.array_len, dtype=complex)
        #self.beta = 1.76/convert(260.0, FEMTO_TO_ARU)
        self.beta = 1.76/self.tau
        self.mu = 0.0#2.0
        self.phase = self.mu*numpy.log(numpy.cosh(self.beta*(self.t - self.t_o)))

        self.rwa = params['run_params'].run_rwa

        self.limit = 2.0*sys.float_info.epsilon

    def setShape(self):
        '''Sets the shape of the pulse'''
        if self.shape == GAUSSIAN:
            self.efield_envelope = self.E_o*numpy.exp(-2.0*numpy.log(2.0)*pow((self.t-self.t_o), 2.0 )/pow(self.tau, 2.0))
            for i in range(self.array_len):
                if self.efield_envelope[i] < self.limit:
                    self.efield_envelope[i] = 0.0
            if (self.rwa == True):
                self.efield = self.efield_envelope*( numpy.exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase))
            else:
                self.efield = 0.5*self.efield_envelope*( numpy.exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) +  numpy.exp(-1j*self.omega_o*(self.t - self.t_o))*numpy.exp(-1j*self.phase) )


        elif self.shape == SECH:
            self.efield_envelope = self.E_o*(1/numpy.cosh(SECH_CONST*(self.t-self.t_o)/self.tau))
            if (self.rwa == True):
                self.efield = self.efield_envelope*(numpy. exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) )
            else:
                self.efield = 0.5*self.efield_envelope*( numpy.exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) + numpy.exp(-1j*self.omega_o*(self.t - self.t_o))*numpy.exp(-1j*self.phase) )
        
       
        elif self.shape == SQUARE:
            for i in range(self.array_len):
                if self.t[i]>=(self.t_o-(self.tau/2)) and self.t[i]<=(self.t_o+(self.tau/2)):
                    self.efield_envelope[i] =self.E_o
                else:
                    self.efield_envelope[i] =0.0
                    
            if (self.rwa == True):
                self.efield = self.efield_envelope*(numpy. exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) )
            else:
                self.efield = 0.5*self.efield_envelope*( numpy.exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) + numpy.exp(-1j*self.omega_o*(self.t - self.t_o))*numpy.exp(-1j*self.phase) )
            
        elif self.shape== LORENTZIAN:
            self.efield_envelope=self.E_o*(1/(1+pow((PI*(self.t-self.t_o)/(2*self.tau)),2)))
            for i in range(self.array_len):
                if self.efield_envelope[i] < self.limit:
                    self.efield_envelope[i] = 0.0
            if (self.rwa == True):
                self.efield = self.efield_envelope*(numpy. exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) )
            else:
                self.efield = 0.5*self.efield_envelope*( numpy.exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) + numpy.exp(-1j*self.omega_o*(self.t - self.t_o))*numpy.exp(-1j*self.phase) )
      
        elif self.shape== DICHROMATIC:
            self.efield_envelope = self.E_o*numpy.exp(-2.0*numpy.log(2.0)*pow((self.t-self.t_o), 2.0 )/pow(self.tau, 2.0))
            for i in range(self.array_len):
                if self.efield_envelope[i] < self.limit:
                    self.efield_envelope[i] = 0.0

            if (self.rwa == True):
                self.efield = self.efield_envelope*(numpy. exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase))*numpy.cos(self.hole_width*self.t)
            else:
                self.efield = 0.5*self.efield_envelope*( numpy.exp(1j*self.omega_o*(self.t - self.t_o))*numpy.exp(1j*self.phase) + numpy.exp(-1j*self.omega_o*(self.t - self.t_o))*numpy.exp(-1j*self.phase) )


        # calculate temporal intensity
        self.intensity = numpy.real(pow(abs(self.efield), 2.0))
        
    def rabifreq(self,dipole):
        self.efield_envelope = abs(self.efield)
        self.rabi_freq = dipole*self.efield_envelope       
    
    def actual_rabifreq(self): # for correct rabi output in rabi_freq.txt
        self.actual_rabi_freq = 0.5*self.dipole*self.efield*(numpy. exp(-1j*self.omega_o*(self.t - self.t_o)))
    
    def detuning(self, omega_o):
        
        self.detuning = self.omega-omega_o
    
#    def rabifreq(self, dipole,delta):
#        self.efield_envelope = abs(self.efield)
#        self.rabi_freq = pow((pow(dipole*self.efield_envelope,2)+pow(delta,2)),0.5)
        
    def fft(self):
        '''Forward FFT of pulse'''
        self.Efield = scipy.fftpack.fft(self.efield)
        for i in range(self.array_len):
            if abs(self.Efield[i]) < self.limit:
                self.Efield[i] = complex(0.0, 0.0)
        self.Efield_phase = numpy.angle(self.Efield)



    def ifft(self):
        '''Inverse FFT of pulse'''
        self.efield = scipy.fftpack.ifft(self.Efield)
        for i in range(self.array_len):
            if abs(self.Efield[i]) < self.limit:
                self.Efield[i] = complex(0.0, 0.0)
        self.intensity = numpy.real(pow(abs(self.efield), 2.0))  # update temporal intensity if efield changes
        self.efield_phase = numpy.angle(self.efield)
        self.inst_freq[0:self.array_len-1] = numpy.diff(numpy.unwrap(self.efield_phase))/numpy.diff(self.t)
        self.inst_freq[self.array_len-1] = self.inst_freq[self.array_len-2]
        self.Efield_phase = numpy.angle(self.Efield) - self.t_o*(abs(self.omega) - self.omega_o)
        self.intensity = numpy.real(pow(abs(self.efield), 2.0))

    def intensityAC(self):
        '''Intensity Autocorrelation'''
        self.intAC = numpy.correlate(self.intensity/self.intensity.max(), self.intensity/self.intensity.max(), 'same')


    def applyMask(self, ampmaskfunc, phasemaskfunc, params, **kwargs):
        '''Apply Mask'''

        # generate amplitude mask using requested mask function
        self.amp_mask = ampmaskfunc(self.omega, self.omega_o, self.amp_mask, params, **kwargs)

        # generate phase mask using requested phase funciton
        self.phase_mask = phasemaskfunc(self.omega, self.omega_o, self.phase_mask, params, **kwargs)

        # apply amplitude and/or phase masks to pulse
        self.Efield = self.Efield*self.amp_mask*numpy.exp(1j*self.phase_mask)
    '''
    def miips(self):
        #Calculates the MIIPS trace
        NPHASE = 100
        NOMEGA =
        self.miips_ref_phase = linspace(0, 2.0*PI, NPHASE)
        self.miips_ref_function = 2.0*PI*sin(gamma*omega - delta)

        self.miips_intensity = zeros(NPHASE,
    '''

    def getPulseArea(self):
        '''Calculates pulse area'''
        return self.dipole*integrate.trapz((self.efield_envelope), x=self.t)/H_BAR_ARU

    def plotefield(self):
        #plot(convert(self.t, ARU_TO_FEMTO), self.efield.real, convert(self.t, ARU_TO_FEMTO), self.efield.imag, convert(self.t, ARU_TO_FEMTO), abs(self.efield))
        plt.plot(convert(self.t, ARU_TO_FEMTO), abs(self.efield)/max(abs(self.efield)))

    def plotrabifreq(self):
        plt.plot(convert(self.t, ARU_TO_FEMTO), 1000.0*convert(self.rabi_freq, ARU_TO_EV))
        plt.xlabel('$t$ (fs)')
        plt.ylabel('$\hbar\Omega$ (meV)')
        plt.grid(True)

    def plotefieldphase(self):
        plt.plot(diff(numpy.unwrap(self.efield_phase)))

    def plotinstfreq(self):
        plt.plot(convert(self.t, ARU_TO_FEMTO), 1000*convert(self.inst_freq - self.omega_o, ARU_TO_EV))
        plt.xlabel('$t$ (fs)')
        plt.ylabel('$\hbar\omega_{inst}-\hbar\omega_{o}$ (meV)')
        plt.grid(True)

    def plotEfield(self):
        temp = scipy.fftpack.fftshift(abs(self.Efield))
        temp = temp/max(temp)
        plt.plot(convert(scipy.fftpack.fftshift(self.omega), ARU_TO_EV), temp)
        plt.xlabel('$\hbar\omega$ (eV)')
        plt.ylabel('$|E(\omega)|$')
        plt.grid(True)
        plt.xlim(convert(self.omega_o, ARU_TO_EV)-0.05,convert(self.omega_o, ARU_TO_EV)+0.05)

    def plotIntEfield(self):
        temp = pow(scipy.fftpack.fftshift(abs(self.Efield)), 2.0)
        temp = temp/max(temp)
        plt.plot(convert(scipy.fftpack.fftshift(self.omega), ARU_TO_EV), temp)
        plt.xlabel('$\hbar\omega$ (eV)')
        plt.ylabel('$|E(\omega)|^2$')
        plt.grid(True)
        plt.xlim(convert(self.omega_o, ARU_TO_EV)-0.05,convert(self.omega_o, ARU_TO_EV)+0.05)

    def plotEfieldPhase(self):
        plt.plot(convert(scipy.fftpack.fftshift(self.omega), ARU_TO_EV), scipy.fftpack.fftshift(numpy.unwrap(self.Efield_phase)))
        #plot(convert(scipy.fftpack.fftshift(self.omega), ARU_TO_EV), scipy.fftpack.fftshift(self.phase_mask))
        plt.xlabel('$\hbar\omega$ (eV)')
        plt.ylabel('$\Phi(\omega)$ (radians)')
        plt.grid(True)
        plt.axis('tight')
        plt.xlim(convert(self.omega_o, ARU_TO_EV)-0.05,convert(self.omega_o, ARU_TO_EV)+0.05)

    def plotIntensity(self):
        plt.plot(convert(self.t, ARU_TO_FEMTO), self.intensity/self.intensity.max())
        plt.xlabel('time (fs)')
        plt.ylabel('intensity')
        plt.grid(True)
        plt.axis('tight')

    def plotIntAC(self):
        plt.plot(convert(self.t, ARU_TO_FEMTO), self.intAC/self.intAC.max())
        plt.xlabel('$t$ (fs)')
        plt.ylabel('IAC')
        plt.grid(True)

    def writePulse(self):
        self.actual_rabifreq()
        '''Write data to files'''
        numpy.savetxt("data/pulse/time.txt", convert(self.t, ARU_TO_FEMTO))
        numpy.savetxt("data/pulse/freq.txt", convert(scipy.fftpack.fftshift(self.omega), ARU_TO_EV))
        numpy.savetxt("data/pulse/efield.txt", self.efield)
        numpy.savetxt("data/pulse/Efield.txt", scipy.fftpack.fftshift(self.Efield))
        numpy.savetxt("data/pulse/IntAC.txt", self.intAC/self.intAC.max())
        temp = pow(abs(self.Efield), 2.0)
        numpy.savetxt("data/pulse/Efield_intensity.txt", scipy.fftpack.fftshift(temp))
        numpy.savetxt("data/pulse/Efield_phase.txt", scipy.fftpack.fftshift(numpy.unwrap(self.Efield_phase)))
        numpy.savetxt("data/pulse/phase.txt", scipy.fftpack.fftshift(self.phase_mask))
        numpy.savetxt("data/pulse/intensity.txt", self.intensity/self.intensity.max())
        numpy.savetxt("data/pulse/rabi_freq.txt", 1000.0*convert(self.actual_rabi_freq, ARU_TO_EV), header='Time dependence of Rabi Frequency (meV)')
#        print self.detuning,"and",self.rabi_freq
        numpy.savetxt("data/pulse/detuning.txt", 1000*convert(self.inst_freq - self.omega_o, ARU_TO_EV), header='Time dependence of detuning (meV)')



    def data(self):
        pulse_data = {"time": convert(self.t, ARU_TO_FEMTO), \
                "efield": convert(self.efield, ARU_TO_ELEC), \
                "freq": convert(scipy.fftpack.fftshift(self.omega), ARU_TO_EV), \
                "Efield": scipy.fftpack.fftshift(self.Efield), \
                "intensity": self.intensity/self.intensity.max(), \
                "AmpMask": scipy.fftpack.fftshift(self.amp_mask), \
                "PhaseMask": scipy.fftpack.fftshift(self.phase_mask)}
        return pulse_data




