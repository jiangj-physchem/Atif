//*********declaration the functions for energy calculation*******//
//****************************************************************//
#ifndef ENERGYCALCULATION_H_
#define ENERGYCALCULATION_H_
void EnergyCalculation(double sigma,double gammab,double f,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                       double* BB,double** pairEner,double** ATT,double** rho,double* Psi,double*** Psi_IJ,
                       double& Ener_tot);
void EnergyCalculation(double sigma,double f,double eta,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                       double** pairEner,double** rho,double* Psi,double& Ener_tot);
#endif //ENERGYCALCULATION_H_
