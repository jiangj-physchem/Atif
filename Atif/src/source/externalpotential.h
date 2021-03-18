//declaration the function for calculating external potential//
//***********************************************************//
#ifndef EXTERNALPOTENTIAL_H_
#define EXTERNALPOTENTIAL_H_
void ExternalPotential(float* uWall,int** depthLR,double** Ext);
void ExternalPotential(double alpha,float* D,float* uWall,int** depthLR,double** Ext);
#endif //EXTERNALPOTENTIAL_H_
