//declaration the function for calculating bulk chemical potential//
//***********************FMT + MSA + TPT1*************************//
#ifndef BULKCHEMPOTENTIALDFT_H_
#define BULKCHEMPOTENTIALDFT_H_
void BulkChemPotentialDFT(double gama,double* rhoB,float* D,double* B,double* BB,float* Z,
                          double** pairEner,double** ATT,double*** Psi_IJ,double* mu,double& P_bulk);
void BulkPureSolvent(double rho_s0,float DS,double pEner,double& mu_s0,double& p_s0);
#endif //BULKCHEMPOTENTIALDFT_H_
