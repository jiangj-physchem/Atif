//declaration the functions for calculating van del Waal**********//
//****************************************************************//
#ifndef INHOMVANDELWAAL_H_
#define INHOMVANDELWAAL_H_
void InhomVanDelWaal(int* LLI,int* ULI,float* D,double** rho,
                     double** pairEner,double** VanDW);

void InhomVanDelWaal(double* rhoB,float* D,double** pairEner,double* VanDWB,double& Van_EN);
#endif //INHOMVANDELWAAL_H_
