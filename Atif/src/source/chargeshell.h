//declaration the functions for electrostatistic correlation******//
//****************************************************************//
#ifndef CHARGESHELL_H_
#define CHARGESHELL_H_
void ChargeShell(int* LLI,int* ULI,double* BB,double** rho,double*** PsiC_IJ,double** DSH);
void ChargeShell(double* rhoB,double* BB,float* Z,double*** PsiC_IJ,double* DSHB);
void ChargeShellDirk(int* LLI,int* ULI,float* D,double** rho,double*** Psi_IJ,double** DSH);
void ChargeShellDirk(double* rhoB,float* D,float* Z,double*** Psi_IJ,double* DSHB);
#endif //CHARGESHELL_H_
