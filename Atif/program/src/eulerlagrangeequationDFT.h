//*****declaration the functions for solving Euler Lagrange*******//
//****************************************************************//
#ifndef EULERLAGRANGEDFT_H_
#define EULERLAGRANGEDFT_H_
#include "clibrary.h"
using namespace std;
void EulerLagrangeDFT(double gama,double phi,double f,double eta_t,int* LLI,int* ULI,float* D,double* BB,
                      float* Z,double* etar,double* Psi,double** pairEner,double* mu,double* rhoB,double** Ext,
                      double** ATT,double*** BesselZero,double*** Psi_IJ,double** rho,double** rho1,string* MODEL);
void EulerLagrangeDFT(double gama,double sigma,double f,double eta_t,double& deltaPhi,double err,int* LLI,int* ULI,float* D,
                   double* BB,double* etar,float* Z,double* Psi,double** pairEner,double* mu,double* rhoB,double** Ext,
                   double** ATT,double*** BesselZero,double*** Psi_IJ,double** rho,double** rho1,string* MODEL);
void EulerLagrangeDFT(double gama,int* LLI,int* ULI,float* D,double* BB,float* Z,double** pairEner,double* mu,
                   double* rhoB,double** ATT,double*** BesselZero,double*** Psi_IJ,string* MODEL);
#endif//EULERLAGRANGEDFT_H_
