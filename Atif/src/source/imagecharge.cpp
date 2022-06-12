//************************function for image charge potential**************************//
#include "clibrary.h"
#include "imagecharge.h"
#include "constantnum.h"

extern double dr;
extern double BJ;
extern double g_size;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int LLIM; //the minimum lower limit of intergral
extern short nspecies;

void ImageChargePotential(double f,float* Z,double** rho,double* u_im)
{
    double  kapax,rhot,errIm,u_im0,kxD1,kxD0,kxX1,kxX0,fm;
    double  R,expm,iterk;
    int     k;
    short   hspecies;
    
    
    
    errIm = 1E-8;
    iterk = 1E4;
    hspecies = nspecies - 1;
    
    for(int i=LLIM; i<=ngrid_m; ++i)
    {
        rhot = 0;
        for(short j=0; j<hspecies; ++j)
        {
            rhot    = rhot + Z[j]*Z[j]*rho[j][i];
        }
        kapax= sqrt(4*Pi*BJ*rhot);
        
        R = dr*i;
        
        kxD0  = exp(kapax*g_size);
        kxD1  = 1.0/kxD0;
        kxX1  = exp(2.0*kapax*R);
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
                u_im[i] = u_im0 + fm*expm/(k*g_size);
            }
            else
            {
                u_im[i] = u_im0 + 0.5*fm*expm*(kxX0*kxD0/((k-1)*g_size+2*R) + kxX1*kxD1/((k+1)*g_size-2*R));
            }
            
            if(fabs(u_im[i]-u_im0) < errIm) break;
            u_im0 = u_im[i];
            ++k;
        }
        while(k < iterk);
        
        if(k >= iterk)
        {
            std::cerr<<"something wrong in image charge: WKB error"<<std::endl;
            exit(0);
        }
        u_im[i] = u_im[i]*BJ;
        
    }
}


void ImageChargePotential(int i,double f,double lamadaKJ,double& u_im)
{
    double  errIm,u_im0,kxD1,kxD0,kxX1,kxX0,fm;
    double  R,expm,iterk;
    int     k;
    
    errIm = 1E-8;
    iterk = 1E4;
    
    R = dr*i;
    
    kxD0  = exp(g_size/lamadaKJ);
    kxD1  = 1.0/kxD0;
    kxX1  = exp(2.0*R/lamadaKJ);
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
            u_im = u_im0 + fm*expm/(k*g_size);
        }
        else
        {
            u_im = u_im0 + 0.5*fm*expm*(kxX0*kxD0/((k-1)*g_size+2*R) + kxX1*kxD1/((k+1)*g_size-2*R));
        }
        
        if(fabs(u_im-u_im0) < errIm) break;
        u_im0 = u_im;
        ++k;
    }
    while(k < iterk);
    
    if(k >= iterk)
    {
        std::cerr<<"something wrong in image charge: WKB error"<<std::endl;
        exit(0);
    }
    
    u_im = u_im*BJ;
    
}
