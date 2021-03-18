//********declaration the functions for PB initialization********//
//****************************************************************//
#ifndef IDEALINITIALIZATION_H_
#define IDEALINITIALIZATION_H_
#include "clibrary.h"
using namespace std;
void Initialization(double gama,double f,double eta_t,int* LLI, int* ULI,float* D,double* BB,double* etar,float* Z,
                    double* rhoB,double* mu,double** Ext,double** ATT,double*** BesselZero,double*** Psi_IJ,double** rho,
                    double** pairEner,string* MODEL,string* MethG);
void Initialization(double& sigma_t,double& dcharge,int& cstep,int& iter,double** rho,ifstream& inFile);
void Initialization(int* LLI, int* ULI,double* rhoB, double** rho);
void Initialization(double& sigma_t,double& dcharge,int& cstep,int& iter,int* LLI, int* ULI,double* rhoB,
                    double** rho,ifstream& inFile);
#endif//IDEALINITIALIZATION_H_
