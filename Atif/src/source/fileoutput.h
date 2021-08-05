//**********declaration the functions for file output*************//
//****************************************************************//
#ifndef FILEOUTPUT_H_
#define FILEOUTPUT_H_
#include "clibrary.h"
using namespace std;
void FileOpen(int num_min,int num_max,ofstream& vari_press,ofstream& converge_rec,ofstream& parameters_s,
              ofstream& contact_theorem,string& fileIterative,string& fileILoop,string& fileTempDen,
              string& fileConverge,string& filePath,string nSurf);
void FileClose(ofstream& vari_press,ofstream& converge_rec,ofstream& parameters_s,ofstream& contact_theorem,string nSurf);
void FileOutPut(int num_file,double lunit,double coe3,double* rhoB,double** rho,double* Psi,string filePath);
void FileOutPut(double& sigma_t,double& dcharge,int& cstep,short& iCode,int iter,
                float* Z,double** rho,string& fileTemp);
void FileOutPut(double& sigma_t,double& dcharge,int& cstep,short& iCode,int iter,int i_num,float* Z,
                double** rho,string& fileILoop,string& fileTemp);
void FileOutPut(double* Psi,double** rho,string& fileTemp);
#endif //FILEOUTPUT_H_
