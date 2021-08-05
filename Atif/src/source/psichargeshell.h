//declaration the functions for electrostatistic correlation******//
//****************************************************************//
#ifndef PSICHARGESHELL_N_H_
#define PSICHARGESHELL_N_H_
void PsiChargeShell(double* B,float* Z,double*** Psi_IJ);
void PsiChargeShellJJ(double* B,double* TB,float* Z,double*** Psi_IJ);
void PsiChargeShellDirk(double gama,float* D,float* Z,double*** Psi_IJ);
#endif //PSICHARGESHELL_N_H_
