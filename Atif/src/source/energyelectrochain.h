//declaration the function for calculating free energy for electrostatic//
//****************** for electrostatic correlation  ***********************//
#ifndef ENERGYELECTROCHAIN_H_
#define ENERGYELECTROCHAIN_H_
void EnergyElectroChain(double gammab,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                        double** ATT,double** rho,double& EN_EC);
#endif //ENERGYELECTROCHAIN_H_
