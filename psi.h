#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _PHONON_SPECTRAL_FUNCTIONS_H
#define _PHONON_SPECTRAL_FUNCTIONS_H

// Some fun constants that are useful

const double e = 2.71828182846;
const double rtpi = 1.7724538509; 


// Some functions that are also useful
__inline double F(double theta, double omega, double omega_c, double omega_L);
__inline double Fpshell(double theta, double omega, double omega_c, double omega_L);
double fp(double base, unsigned int exp);
__inline double sq(double x);
__inline double cu(double x);
__inline double p4(double x);
double integrate(unsigned int N, double omega, double omega_c, double omega_L);



#endif
