#! /home/ajan/anaconda/bin/python
from pylab import *

F = loadtxt("freq_1.txt")
P = loadtxt("phase_1.txt")

alpha = 0.996025000877*pi
gamma = 307.117675056*1E-15
delta = -0.57506777497*pi

E_o = 1.07175
hbar = 1.054571596E-34
eV = 1.602176462E-19 

omega_1 = 1.0*eV/hbar
omega_2 = 1.2*eV/hbar
omega_o = E_o*eV/hbar

omega = linspace(omega_1, omega_2, 2000)


y = alpha*cos(gamma*(omega-omega_o)-delta)


plot(omega*hbar/eV, y)

#print f_o
#print gamma
#y = alpha*cos(2*pi*gamma*(f-f_o)-delta)
#y = alpha*cos(2*pi*gamma*(freq-f_o)-delta)

#plot(y)
plot(F, P, 'r')
xlim(1.0, 1.2)
grid('on')
show()
