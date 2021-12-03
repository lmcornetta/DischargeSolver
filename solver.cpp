/***********************************************************************
 !    Solver: source file
 !       Main function           
 !
 !       The integration is based on a RK4 solver with adaptative
 !       step and local extrapolation. See header file. 
 !
 !                  - Lucas Cornetta 09/OCT/2021
 **********************************************************************/

#include "solver.hpp"
#include <iostream>
#include <fstream>
#include <string.h>

int main(int argc, char * argv[]) {

    double M = atof(argv[1])*M0;
    if (M <= .0) { 
        std::cout << "Black-hole with M less or equal to 0. Abort...";
        exit(EXIT_FAILURE);
    }   
    //double TF = (M_PI)*(1.46e82*pow((M/M0),3)); 
    //double magnitudeFactor = M;
    double Q = atof(argv[2])*M0;
    double QoverM2 = pow((Q/M),2);

    double t = .0; 
    double dt = 1e-20;

    // RK4 constants and function approximators
    double a14 = 0.1666;
    double a23 = 0.3333;
    double K11, K21, K31, K41;
    double L11, L21, L31, L41;
    double K12, K22, K32, K42;
    double L12, L22, L32, L42;
    double K13, K23, K33, K43;
    double L13, L23, L33, L43;

    // formatting output
    std::ofstream fout (strcat(argv[3],".dat"));
    if (!fout) {
        std::cout << "\n error: could not open file" << std::endl;
        exit(EXIT_FAILURE);
    }   
    fout << std::scientific;
    fout.precision(15);
    fout << "#Time (s)" << "\t\t\t\t" << "Mass (g. u.)" << "\t\t\t" <<
             "Charge (g. u.)" << "\t\t\t" << "Q^2/M^2" << "\n";
    fout << t << "\t" << M << "\t" << Q << "\t" << QoverM2 << "\n";

    while (M > .0) {
    
        if (extreme(M,Q)) {
            double dMtmp = dM(M,Q);
            double dQtmp = dQ(M,Q);
            double h = 1./fabs(dMtmp - dQtmp);
            while ((Q + h*dQtmp)/(M + h*dMtmp) >= 1.) {
                h *= 1 - 1e-8;
            }
            M += h*dMtmp; M = std::max(M,.0);
            Q += h*dQtmp; Q = std::max(Q,.0);
            dt = h; t += h;
        }

        else {
            double ratio = 1.0 + 1e-10;
            double dM1, dM2, dM3, dMtmp = .0; 
            double dQ1, dQ2, dQ3, dQtmp = .0;    
            while (ratio >= 1.0 + 1e-10) {

                dt /= ratio;

                // 1:
                // First dt step
                K11 = dt*dM(M,Q);
                L11 = dt*dQ(M,Q);
                K21 = dt*dM(M + 0.5*K11, Q + 0.5*L11);
                L21 = dt*dQ(M + 0.5*K11, Q + 0.5*L11);
                K31 = dt*dM(M + 0.5*K21, Q + 0.5*L21);
                L31 = dt*dQ(M + 0.5*K21, Q + 0.5*L21);
                K41 = dt*dM(M + K31, Q + L31);
                L41 = dt*dQ(M + K31, Q + L31);

                dM1 = a14*(K11 + K41) + a23*(K21 + K31);
                dQ1 = a14*(L11 + L41) + a23*(L21 + L31);

                // Second dt step
                K12 = dt*dM(M + dM1,Q + dQ1);
                L12 = dt*dQ(M + dM1,Q + dQ1);
                K22 = dt*dM(M + dM1 + 0.5*K12, Q + dQ1 + 0.5*L12);
                L22 = dt*dQ(M + dM1 + 0.5*K12, Q + dQ1 + 0.5*L12);
                K32 = dt*dM(M + dM1 + 0.5*K22, Q + dQ1 + 0.5*L22);
                L32 = dt*dQ(M + dM1 + 0.5*K22, Q + dQ1 + 0.5*L22);
                K42 = dt*dM(M + dM1 + K32, Q + dQ1 + L32);
                L42 = dt*dQ(M + dM1 + K32, Q + dQ1 + L32);

                dM2 = dM1 + a14*(K12 + K42) + a23*(K22 + K32);
                dQ2 = dQ1 + a14*(L12 + L42) + a23*(L22 + L32);

                // 2:
                // Unique 2dt step
                K13 = 2*dt*dM(M,Q);
                L13 = 2*dt*dQ(M,Q);
                K23 = 2*dt*dM(M + 0.5*K13, Q + 0.5*L13);
                L23 = 2*dt*dQ(M + 0.5*K13, Q + 0.5*L13);
                K33 = 2*dt*dM(M + 0.5*K23, Q + 0.5*L23);
                L33 = 2*dt*dQ(M + 0.5*K23, Q + 0.5*L23);
                K43 = 2*dt*dM(M + K33, Q + L33);
                L43 = 2*dt*dQ(M + K33, Q + L33);

                dM3 = a14*(K13 + K43) + a23*(K23 + K33);
                dQ3 = a14*(L13 + L43) + a23*(L23 + L33);

                // Error evaluation
                double epsM = (dM2 - dM3)/30;
                double epsQ = (dQ2 - dQ3)/30;
                double error = std::max(fabs(epsM),fabs(epsQ));
                ratio = pow(error/precision(Q,M),0.25);
            }

            t += 2*dt;
            dt = std::min(dt/(ratio + 1e-30),2*dt);
            dMtmp = dM2 + (dM2 - dM3)/15;
            dQtmp = dQ2 + (dQ2 - dQ3)/15;
            M += dMtmp; M = std::max(M,.0);
            Q += dQtmp; Q = std::max(Q,.0);
        }

        if (M <= .0) break;
        QoverM2 = pow((Q/M),2);
        fout << t/timeConversion << "\t" << M << "\t" << Q << "\t" << QoverM2 << "\n";
    }

    t += 2*dt;
    fout << t/timeConversion << "\t" << M << "\t" << Q << "\t" << QoverM2 << "\n";
    return 0;
}
