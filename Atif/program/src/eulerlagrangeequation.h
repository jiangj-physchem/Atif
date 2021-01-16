//*****declaration the functions for solving Euler Lagrange*******//
//****************************************************************//
#ifndef EULERLAGRANGE_H_
#define EULERLAGRANGE_H_
#include "clibrary.h"
using namespace std;
void EulerLagrange(double phi,double f,double eta_t,int* LLI,int* ULI,float* D,float* Z,double* etar,
                   double* Psi,double** pairEner,double* mu,double* rhoB,double** Ext,
                   double*** BesselZero,double** rho,double** rho1,string* MODEL);
void EulerLagrange(double sigma,double f,double eta_t,double& deltaPhi,double err,int* LLI,int* ULI,float* D,
                   double* etar,float* Z,double* Psi,double** pairEner,double* mu,double* rhoB,double** Ext,
                   double*** BesselZero,double** rho,double** rho1,string* MODEL);
void EulerLagrange(int* LLI,int* ULI,float* D,float* Z,double** pairEner,double* mu,
                   double* rhoB,double*** BesselZero,string* MODEL);
#endif//EULERLAGRANGE_H_
