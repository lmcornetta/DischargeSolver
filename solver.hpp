/***********************************************************************
 !    Solver: header file
 !       Declaration of constants and  signature/definition of
 !       functions dM (dM/dt), dQ (dQ/dt) and precision.
 !
 !       The integration is based on a RK4 solver with adaptative
 !       step and local extrapolation. The non-constant precision
 !       decreases ('better' precision) for black holes near 
 !       the extreme condition
 !
 !                  - Lucas Cornetta 09/OCT/2021
 **********************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>

//Global variables in geometric units and conversions
#define M0             1477.8                       // solar mass
#define Q0             1.7e5*M0                     // charge parameter (eq 8 from HW)
#define HBAR           2.61e-70                     // Planck's constant
#define QE             1.379e-36                    // electron charge
#define ME             6.76e-58                     // electron mass
#define A              pow(M_PI,2)/(15*pow(HBAR,3)) // 'a' parameter from eq 12 of HW (Hawking radiation)
#define ALPHA          2.0228                       // alpha parameter from HW. By default we are using n=3
#define timeConversion 9.4728e15                    // from geometric to yr

double dM(double M, double Q);
double dQ(double M, double Q);
double precision(double M, double Q);

bool extreme(double M, double Q);
double dMExtreme(double M, double Q);
double dQExtreme(double M, double Q);
double precisionExtreme(double M, double Q);

double dM(double M, double Q) {/* Derivative dM/dt */
    if (fabs(Q) >= M || M <= .0) { return .0; }
    double rp = M + sqrt(pow(M,2)-pow(Q,2));
    double k = sqrt(pow(M,2)-pow(Q,2))/pow(rp,2);
    double T = HBAR*k/(2*M_PI);
    double sigma0 = M_PI*pow((3*M + sqrt(9*pow(M,2)-8*pow(Q,2))),4)
                    /(8*(3*pow(M,2)-2*pow(Q,2)+M*sqrt(9*pow(M,2)-8*pow(Q,2))));
    double Mp = -pow(T,2)*A*ALPHA*sigma0*pow(T,2) + Q*dQ(M,Q)/rp;
    return Mp;
}

double dQ(double M, double Q) {/* Derivative dQ/dt */
    if (fabs(Q) >= M || M < .0 || Q <= .0) { return .0; }
    double rp = M + sqrt(pow(M,2)-pow(Q,2));
    double Qp = -pow(QE,4)*pow(Q,3)*exp(-pow(rp,2)/(Q*Q0))/
                (2*pow(M_PI,3)*HBAR*pow(ME,2)*pow(rp,3));
    return Qp;
}

double precision(double M, double Q) {
    return 1e-1;
}

bool extreme(double M, double Q) {
    double z = Q/M;
    if (fabs(1 - pow(z,2)) < 1e-5) return true;
    return false;
}

double dMExtreme(double M, double Q) {
    double z = Q/M;
    double T = 0.5*HBAR*(sqrt(1 - pow(z,2)) - (1 - pow(z,2)))/M_PI;
    double Mp = -pow(T,2)*A*ALPHA*8*M_PI*pow(M,2)*pow(T,2) + 0.5*z*dQExtreme(M,Q);
    return Mp;
}

double dQExtreme(double M, double Q) {
    double z = Q/M;
    double Qp = -pow(QE,4)*pow(z,3)*exp(-4*z*M/Q0)/(16*M_PI*HBAR*pow(ME,2));
    return Qp;
}

double precisionExtreme(double M, double Q) {
    double z = Q/M;
    return sqrt(1-pow(z,12)) + 1e-20;
}

#endif
