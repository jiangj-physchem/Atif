//declaration the functions for calculating weighted density******//
//************** from FMT (see Rosenfeld, PRL, 1989 **************//
//****************************************************************//
#ifndef WEIGHTEDDENSITY_H_
#define WEIGHTEDDENSITY_H_
void WeightedDensity(int i,int* LLI,int* ULI,float* D,double** rho,double** NI);
void WeightedDensity(float* D,double* rhoB,double** NIB);
void WeightedDensity(int i,int* LLI,int* ULI,double* B,double** H,double** rho,double** QI);
void WeightedDensity(double gammab,double* rhoB,double* B,double** H,double** QIB);
void WeightedDensityDirk(int i,int* LLI,int* ULI,double* B,double** H,double** rho,double** QI);
void WeightedDensityDirk(double gammab,double* rhoB,double* B,double** H,double** QIB);
void WeightedDensityV(int i,int* LLI,int* ULI,double* B,double** H,double** rho,double** QI);
void WeightedDensityV(double gammab,double* rhoB,double* B,double** H,double** QIB);
void WeightedDensityDirk(int i,int* LLI,int* ULI,double* B,double** rho,double* QI);
void WeightedDensityDirk(double gammab,double* rhoB,double* B,double* QIB);
void WeightedDensityV(int i,int* LLI,int* ULI,double* B,double** rho,double* QI);
void WeightedDensityV(double gammab,double* rhoB,double* B,double* QIB);

#endif //WEIGHTEDDENSITY_H_
