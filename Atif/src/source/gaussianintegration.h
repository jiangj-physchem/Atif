//declaration the functions for GAUSSIAN Inergration//
//*****************************************************************//
#ifndef GAUSSIANINTEGRAL_H_
#define GAUSSIANINTEGRAL_H_
double GaussianIntegrationZ0(double* integrandF,int aa,int bb,int ia,int ib);
double GaussianIntegrationZ1(double* integrandF,int aa,int bb,int ia,int ib);
double GaussianIntegrationZ2(double* integrandF,int aa,int bb,int ia,int ib);
double GaussianIntegration(double* integrandF1,double* integrandF2,int aa,int bb,int ia,int ib);
double GaussianIntegration(double* integrandF1,double* integrandF2,int orient0,int ia,int ib);
double GaussianIntegrationShell(double* integrandF,double* Psi_val,int aa,int bb,int ia,int ib);
#endif //GAUSSIANINTEGRAL_H_
