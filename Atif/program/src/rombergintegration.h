//declaration the functions for Romberg Inergration//
//*****************************************************************//
#ifndef ROMBERGINTEGRAL_H_
#define ROMBERGINTEGRAL_H_
double RombergIntegration(double* integrandF,int aa,int bb,int ia,int ib);
double RombergIntegration(double* integrandF1,double* integrandF2,int aa,int bb,int ia,int ib);
double RombergIntegration(double* integrandF,int aa,int bb,int ia,int ib,int i_R);
double RombergIntegration(double* integrandF,int aa,int bb,int ia,int ib,double D2_j,int i_R);
double RombergIntegration(int aa,int bb,int ia,int ib,double D2_j);
double RombergIntegration(double* integrandF,double basic_dr,short num_deep,int basic_n,int iaa,int ibb);
double RombergIntegration(double* integrandF1,double* integrandF2,double basic_dr,short num_deep,int basic_n,int iaa,int ibb);
double RombergIntegration(double* integrandF,double basic_dr,short num_deep,int basic_n,int iaa,int ibb,int i_R);
double RombergIntegration(double* integrandF,double basic_dr,short num_deep,int basic_n,int iaa,int ibb,double D2_j,int i_R);
double RombergIntegration(double basic_dr,short num_deep,int basic_n,int iaa,int ibb,double D2_j);

#endif //ROMBERGINTEGRAL_H_
