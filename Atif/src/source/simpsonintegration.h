//declaration the functions for Simpson Inergration//
//*****************************************************************//
#ifndef SIMPSONINTEGRAL_H_
#define SIMPSONINTEGRAL_H_
double SimpsonIntegration(double* integrandF,int aa,int bb,int ia,int ib);
double SimpsonIntegration(double* integrandF,int aa,int bb,int ia,int ib,int i_R);
double SimpsonIntegration(double* integrandF,int aa,int bb,int ia,int ib,double D2_j,int i_R);
double SimpsonIntegration(int aa,int bb,int ia,int ib,double D2_j);
double SimpsonIntegration(double* integrandF1,double* integrandF2,int aa,int bb,int ia,int ib);
double SimpsonIntegration(double* integrandF1,double* integrandF2,int orient0,int ia,int ib);
#endif //SIMPSONINTEGRAL_H_
