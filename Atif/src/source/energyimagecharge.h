//*******declaration the functions for image charge energy*********//
//*****************************************************************//
#ifndef ENERGYIMAGECHARGE_H_
#define ENERGYIMAGECHARGE_H_
void ImageChargeEnergy(int i,double f,float* Z,double** rho,double& f_im);
void ImageChargeEnergy(int i,double f,double lamadaKJ,double& f_im);
void FunctionJ(int i,double f,double alpkapa,double& jka);
#endif //ENERGYIMAGECHARGE_H_
