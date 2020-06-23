#! /home/ajan/anaconda/bin/python

from pylab import *

# load data from my calculations
pyshape_folder = './pyshape_data/'
Omega_eid = loadtxt(pyshape_folder + 'Omega_eid.txt', delimiter=',', unpack=True, skiprows=1)
K_real_eid_00 = loadtxt(pyshape_folder + 'K_real_eid_00K.txt', delimiter=',', unpack=True, skiprows=1)
K_imag_eid_00 = loadtxt(pyshape_folder + 'K_imag_eid_00K.txt', delimiter=',', unpack=True, skiprows=1)

K_real_eid_10 = loadtxt(pyshape_folder + 'K_real_eid_10K.txt', delimiter=',', unpack=True, skiprows=1)
K_imag_eid_10 = loadtxt(pyshape_folder + 'K_imag_eid_10K.txt', delimiter=',', unpack=True, skiprows=1)

K_real_eid_20 = loadtxt(pyshape_folder + 'K_real_eid_20K.txt', delimiter=',', unpack=True, skiprows=1)
K_imag_eid_20 = loadtxt(pyshape_folder + 'K_imag_eid_20K.txt', delimiter=',', unpack=True, skiprows=1)

K_real_eid_30 = loadtxt(pyshape_folder + 'K_real_eid_30K.txt', delimiter=',', unpack=True, skiprows=1)
K_imag_eid_30 = loadtxt(pyshape_folder + 'K_imag_eid_30K.txt', delimiter=',', unpack=True, skiprows=1)

K_real_eid_40 = loadtxt(pyshape_folder + 'K_real_eid_40K.txt', delimiter=',', unpack=True, skiprows=1)
K_imag_eid_40 = loadtxt(pyshape_folder + 'K_imag_eid_40K.txt', delimiter=',', unpack=True, skiprows=1)

# load data from Ramsay calculations
Ramsay_folder = './Ramsay_JAP109_102415_2011_data/'
ROmega_eid_00, RK_real_eid_00 = loadtxt(Ramsay_folder + 'K_real_00K.txt', delimiter='\t', unpack=True, skiprows=1)
ROmega_eid_10, RK_real_eid_10 = loadtxt(Ramsay_folder + 'K_real_10K.txt', delimiter='\t', unpack=True, skiprows=1)
ROmega_eid_20, RK_real_eid_20 = loadtxt(Ramsay_folder + 'K_real_20K.txt', delimiter='\t', unpack=True, skiprows=1)
ROmega_eid_30, RK_real_eid_30 = loadtxt(Ramsay_folder + 'K_real_30K.txt', delimiter='\t', unpack=True, skiprows=1)
ROmega_eid_40, RK_real_eid_40 = loadtxt(Ramsay_folder + 'K_real_40K.txt', delimiter='\t', unpack=True, skiprows=1)

IOmega_eid_00, RK_imag_eid_00 = loadtxt(Ramsay_folder + 'K_imag_00K.txt', delimiter='\t', unpack=True, skiprows=1)
IOmega_eid_10, RK_imag_eid_10 = loadtxt(Ramsay_folder + 'K_imag_10K.txt', delimiter='\t', unpack=True, skiprows=1)
IOmega_eid_20, RK_imag_eid_20 = loadtxt(Ramsay_folder + 'K_imag_20K.txt', delimiter='\t', unpack=True, skiprows=1)
IOmega_eid_30, RK_imag_eid_30 = loadtxt(Ramsay_folder + 'K_imag_30K.txt', delimiter='\t', unpack=True, skiprows=1)
IOmega_eid_40, RK_imag_eid_40 = loadtxt(Ramsay_folder + 'K_imag_40K.txt', delimiter='\t', unpack=True, skiprows=1)


subplot(2,1,1)
plot(Omega_eid, K_real_eid_00, 'k-', label='0 K')
plot(Omega_eid, K_real_eid_10, 'r-', label='10 K')
plot(Omega_eid, K_real_eid_20, 'g-', label='20 K')
plot(Omega_eid, K_real_eid_30, 'b-', label='30 K')
plot(Omega_eid, K_real_eid_40, 'm-', label='40 K')
'''
plot(ROmega_eid_00, RK_real_eid_00, 'ko')
plot(ROmega_eid_10, RK_real_eid_10, 'ro')
plot(ROmega_eid_20, RK_real_eid_20, 'go')
plot(ROmega_eid_30, RK_real_eid_30, 'bo')
plot(ROmega_eid_40, RK_real_eid_40, 'mo')
'''
xlim(0, 5.0)
ylim(0.0, 1.0)
grid(True)
xlabel(r'Rabi Energy $\hbar \Omega$ (meV)')
ylabel(r'Re$[\hbar K(\hbar \Omega)]$ (ps$^{-1}$)')
legend()

subplot(2,1,2)
plot(Omega_eid, K_imag_eid_00, 'k-', label='0 K')
plot(Omega_eid, K_imag_eid_10, 'r-', label='10 K')
plot(Omega_eid, K_imag_eid_20, 'g-', label='20 K')
plot(Omega_eid, K_imag_eid_30, 'b-', label='30 K')
plot(Omega_eid, K_imag_eid_40, 'm-', label='40 K')
''' 
plot(IOmega_eid_00, RK_imag_eid_00, 'ko')
plot(IOmega_eid_10, RK_imag_eid_10, 'ro')
const = 1.00

# the pyshape data matches the Ramasy data if the 
# pulse area is scaled by 1.44 for 20, 30, and 40 K
plot(const*IOmega_eid_20, RK_imag_eid_20, 'go')
plot(const*IOmega_eid_30, RK_imag_eid_30, 'bo')
plot(const*IOmega_eid_40, RK_imag_eid_40, 'mo')
'''
xlim(0.0, 5.0)
ylim(-0.3, 0.5)
grid(True)
xlabel(r'Rabi Energy $\hbar \Omega$ (meV)')
ylabel(r'Im$[\hbar K(\hbar \Omega)]$ (meV)')
legend()

show()
