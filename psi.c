#include "psi.h"

__inline double F(double theta, double omega, double omega_c, double omega_L) {
    return exp(-1*sq(omega)*(sq(sin(theta)/omega_c) + sq(cos(theta)/omega_L)))*sin(theta);
}

__inline double Fpshell(double theta, double omega, double omega_c, double omega_L) {
    double st=sin(theta);
    double cs=cos(theta);
    return exp(-1*sq(omega)*(sq(st/omega_c) + sq(cs/omega_L)))*st*(1-(2*sq(omega/omega_L)*sq(cs)) + p4(omega/omega_L)*p4(cs));
}

double fp(double base, unsigned int exp) {
    double nb = 1;
    while (exp > 0) {
        nb *= base;
        exp --;
    }
    return nb;
}

__inline double sq(double x) {return x*x;}
__inline double cu(double x) {return x*x*x;}
__inline double p4(double x) {return x*x*x*x;}

double integrate(register unsigned int N, double omega, double omega_c, double omega_L) {
    double sum = 0;
    const double pi = 3.1415926535897; 
    double step = (pi)/((double)N);

    for (register unsigned int i = 1; i < N; i ++) {
        sum += 0.5*(F( ((double)i-1)*step, omega, omega_c, omega_L ) 
                  + F( (double)i*step, omega, omega_c, omega_L) )*step;
    }

    return cu(omega)*sum;
}

double integratepshell(register unsigned int N, double omega, double omega_c, double omega_L) {
    double sum = 0;
    const double pi = 3.1415926535897; 
    double step = (pi)/((double)N);

    for (register unsigned int i = 1; i < N; i ++) {
        sum += 0.5*(Fpshell( ((double)i-1)*step, omega, omega_c, omega_L ) 
                  + Fpshell( (double)i*step, omega, omega_c, omega_L) )*step;
    }

    return cu(omega)*sum;
}
