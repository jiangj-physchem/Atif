//****************energy for image charge**********************//
#include "clibrary.h"
#include "energyimagecharge.h"
#include "constantnum.h"

extern double dr;
extern double BJ;
extern double g_size;
extern short  nspecies;

void ImageChargeEnergy(int i,double f,float* Z,double** rho,double& f_im)
{
    double  jka,alpkapa,dalpha,rhot,alpha0,kapax,kapax2;
    int     nalpha;
    short   hspecies;
    
    hspecies= nspecies - 1;
    
    
    rhot = 0;
    for(short j=0; j<hspecies; ++j)
    {
        rhot    = rhot + Z[j]*Z[j]*rho[j][i];
    }
    kapax2= BJ*rhot;
    kapax = sqrt(4*Pi*kapax2);
    
    dalpha = 0.01;
    nalpha = round(1.0/dalpha);
    
    FunctionJ(i,f,kapax,jka);
    f_im = jka;
    for(int i1=1; i1<nalpha; ++i1)
    {
        alpha0 = i1*dalpha;
        alpkapa= kapax*alpha0;
        FunctionJ(i,f,alpkapa,jka);
        
        if(i1%2 != 0)
        {
            f_im = f_im + alpha0*jka*4.0;
        }
        else
        {
            f_im = f_im + alpha0*jka*2.0;
        }
    }
    
    f_im = kapax2*f_im*dalpha/3.0;

}



void ImageChargeEnergy(int i,double f,double lamadaKJ,double& f_im)
{
    double  jka,alpkapa,dalpha,alpha0,kapax,kapax2;
    int     nalpha;
    
    kapax = 1.0/lamadaKJ;
    kapax2= kapax*kapax/(4*Pi);
    
    
    dalpha = 0.01;
    nalpha = round(1.0/dalpha);
    
    FunctionJ(i,f,kapax,jka);
    f_im = jka;
    for(int i1=1; i1<nalpha; ++i1)
    {
        alpha0 = i1*dalpha;
        alpkapa= kapax*alpha0;
        FunctionJ(i,f,alpkapa,jka);
        
        if(i1%2 != 0)
        {
            f_im = f_im + alpha0*jka*4.0;
        }
        else
        {
            f_im = f_im + alpha0*jka*2.0;
        }
    }
    
    f_im = kapax2*f_im*dalpha/3.0;
    
}


void FunctionJ(int i,double f,double alpkapa,double& jka)
{
    double  errIm,u_im0,kxD1,kxD0,kxX1,kxX0,fm;
    double  R,expm,iterk;
    int     k;
    
    errIm = 1E-10;
    iterk = 1E5;
    
    
    R = dr*i;
    
    kxD0  = exp(alpkapa*g_size);
    kxD1  = 1.0/kxD0;
    kxX1  = exp(2.0*alpkapa*R);
    kxX0  = 1.0/kxX1;
    
    
    fm  = f;
    expm= kxD1;
    
    u_im0 = fm*(0.5*kxX0/R + 0.5*kxX1*expm*kxD1/(g_size-R));
    
    
    k   = 2;
    do
    {
        fm  = f*fm;
        expm= expm*kxD1;
        
        if(k%2 ==0)
        {
            jka = u_im0 + fm*expm/(k*g_size);
        }
        else
        {
            jka = u_im0 + 0.5*fm*expm*(kxX0*kxD0/((k-1)*g_size+2*R) + kxX1*kxD1/((k+1)*g_size-2*R));
        }
        
        if(fabs(jka-u_im0) < errIm) break;
        u_im0 = jka;
        ++k;
    }
    while(k < iterk);
    
    if(k >= iterk)
    {
        std::cerr<<"something wrong in image charge energy: J(ak) error"<<std::endl;
        exit(0);
    }
}
