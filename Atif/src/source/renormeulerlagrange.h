//*****declaration the functions for solving Euler Lagrange*******//
//****************************************************************//
#ifndef RENORMEULERLAGRANGE_H_
#define RENORMEULERLAGRANGE_H_
#include "clibrary.h"
using namespace std;
double RenormEulerLagrange(double ratio,double rhoBM1,double rhoBM2,float* D,short* MB1,short* MB2,float* Z,
                           double** ff1,double** ff2,double** rho1,double** rho2,double*** BesselZero,string* MODEL);
#endif//RENORMEULERLAGRANGE_H_
