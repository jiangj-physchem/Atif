//***********solve the possion equation****************//
#include "clibrary.h"
#include "volumefraction.h"
#include "simpsonintegration.h"
#include "constantnum.h"

extern double dr;
extern short nspecies;
extern int LLIM; //the minimum lower limit of intergral
extern int ngrid; //the number of grids
extern int ngrid_m; //the number of grids: the middle between two surfaces

void VolumeFraction(double eta_t,float* D,double* etar,double** rho)
{
    short  hspecies;
    double eta0;
    double* D3;
    double* eta;
    
    hspecies = nspecies - 1;
    
    eta  = new double[hspecies]();
    D3   = new double[nspecies]();
    
    for(short i=0; i<nspecies; ++i)
    {
        D3[i] = D[i]*D[i]*D[i]*Pi/6.0;
    }
    
    for(int k=LLIM; k<=ngrid_m; ++k)
    {
        for(short i=0; i<hspecies; ++i)
        {
            eta[i] = rho[i][k]*D3[i];
        }

        
        eta0 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            eta0 += eta[i];
        }
        
        if(eta0 > eta_t)
        {
            for(short i=0; i<hspecies; ++i)
            {
                if(eta[i] > etar[i])
                {
                    rho[i][k] = rho[i][k]*etar[i]/eta[i];
                    rho[i][ngrid-k] = rho[i][k];
                }
            }
        }
        
        eta0= 0;
        for(short i=0; i<hspecies; ++i)
        {
            eta0 += rho[i][k]*D3[i];
        }
        rho[hspecies][k] = (eta_t - eta0)/D3[hspecies];
        rho[hspecies][ngrid-k] = rho[hspecies][k];
        
        if(rho[hspecies][k] < 0)
        {
            std::cerr<<"Error in VolumeFraction.cpp: solvent concentration is negative"<<std::endl;
            exit(0);
        }
    }
    
    delete [] eta;
    delete [] D3;
}

void VolumeFraction(int i,int* LLI,int* ULI,float* D,double** rho,double* eta)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double* N3I;
    double  R,rin,rip,D2,temp;
    int kin,kip,hspecies;
    
    hspecies = nspecies - 1;
    
    N3I   = new double[hspecies]();
    R   = i*dr;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - 0.5*D[j];
        rip = R + 0.5*D[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        D2  = D[j]*D[j]*0.25;
        
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        temp  =SimpsonIntegration(rho[j],0,ngrid,kin,kip,D2,i);
        N3I[j]=Pi*temp;
        eta[j] = N3I[j];
    }
    
    
    delete [] N3I;
}


void VolumeFractionDFT(double eta_t,float* D,int* LLI,int* ULI,double* etar,double** rho)
{
    short  hspecies;
    double eta0;
    double* eta;
    
    hspecies = nspecies - 1;
    
    eta  = new double[hspecies]();
    
    
    for(int k=LLIM; k<=ngrid_m; ++k)
    {
        VolumeFraction(k,LLI,ULI,D,rho,eta);
        
        
        eta0 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            eta0 += eta[i];
        }
        
        if(eta0 > eta_t)
        {
            for(short i=0; i<hspecies; ++i)
            {
                if(eta[i] > etar[i])
                {
                    rho[i][k] = rho[i][k]*etar[i]/eta[i];
                    rho[i][ngrid-k] = rho[i][k];
                }
            }
        }
    }
    
    delete [] eta;
}
