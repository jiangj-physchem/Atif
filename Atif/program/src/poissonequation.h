//declaration the functions for solving poisson euqation**********//
//****************************************************************//
#ifndef POISSONEQUATION_H_
#define POISSONEQUATION_H_
void PoissonEquation(float* Z,double** rho,double phiL,double phiR,double* Psi);
void PoissonEquation(float* Z,double** rho,double* Psi);
void PoissonEquation(float* Z,double** rho,double phi,double* Psi);
void PoissonEquation(float* Z,double** rho,double* Psi1,double* Psi);
void PoissonEquation(double sigma,float* Z,double** rho,double* Psi);
#endif //POISSONEQUATION_H_
