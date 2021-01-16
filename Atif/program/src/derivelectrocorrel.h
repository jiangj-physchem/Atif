//declaration the function for calculating inhomogeneous chemical potential//
//****************** for electrostatic correlation  ***********************//
#ifndef DERIVELECTROCORREL_H_
#define DERIVELECTROCORREL_H_
void DerivElectroCorrel(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                        double** rho,double** DES);
void DerivElectroCorrel(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN);
void DerivElectroCorrelDirk(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                        double** rho,double** DES);
void DerivElectroCorrelDirk(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN);
void DerivElectroCorrelV(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                        double** rho,double** DES);
void DerivElectroCorrelV(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN);
#endif //DERIVELECTROCORREL_H_
