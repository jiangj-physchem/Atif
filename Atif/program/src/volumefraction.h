//declaration the functions for calculating volume fraction ******//
//****************************************************************//
#ifndef VOLUMEFRACTION_H_
#define VOLUMEFRACTION_H_
void VolumeFraction(double eta_t,float* D,double* etar,double** rho);
void VolumeFraction(int i,int* LLI,int* ULI,float* D,double** rho,double* eta);
void VolumeFractionDFT(double eta_t,float* D,int* LLI,int* ULI,double* etar,double** rho);
#endif //VOLUMEFRACTION_H_
