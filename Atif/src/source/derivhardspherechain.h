//declaration the function for calculating inhomogeneous chemical potential//
//********************* for hard sphere chain  ****************************//
//*****************************FMT + TPT1**********************************//
#ifndef DERIVHARDSPHERECHAIN_H_
#define DERIVHARDSPHERECHAIN_H_
void DerivHardSphereChain(int MAXR,int* LLI,int* ULI,double* rhoB,float* D,float* Z,double** ATT,
                          double** rho,double** DCH);
void DerivHardSphereChain(double* rhoB,float* D,float* Z,double** ATT,double* DCHB,double& HCH_EN);
#endif //DERIVHARDSPHERECHAIN_H_
