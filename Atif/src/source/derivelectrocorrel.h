//declaration the function for calculating inhomogeneous chemical potential//
//****************** for electrostatic correlation  ***********************//
#ifndef DERIVELECTROCORREL_H_
#define DERIVELECTROCORREL_H_
void DerivElectroCorrel(double gammab,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                        double** ATT,double** rho,double** DES);
void DerivElectroCorrel(double* rhoB,double gammab,float* D,float* Z,double** ATT,
                        double* DESB,double& ES_EN);
#endif //DERIVELECTROCORREL_H_
