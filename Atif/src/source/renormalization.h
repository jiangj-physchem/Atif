//*******declaration the functions for charge neutrality**********//
//****************************************************************//
#ifndef RENORMALIZATION_H_
#define RENORMALIZATION_H_
void RenormDensity(double sigma,double& deltaPhi,double err,double rhoBM1,double rhoBM2,float* D,
                   short* MB1,short* MB2,float* Z,double** ff1,double** ff2,double** rho1,
                   double*** BesselZero,string* MODEL);
void RenormDensity_Newton(double sigma,double& deltaPhi,double err,double rhoBM1,double rhoBM2,float* D,
                          short* MB1,short* MB2,float* Z,double** ff1,double** ff2,double** rho1,
                          double*** BesselZero,string* MODEL);
void RenormDensity_Newton(double sigma,double err,double rhoBM1,double rhoBM2,float* Z,
                          double** rho1);
void Renorm_Newton_Downhill(double sigma,double err,double rhoBM1,double rhoBM2,float* Z,
                            double** rho1);
#endif//RENORMALIZATION_H_
