//***********potential from electrostatistic correlation****************//
#include "clibrary.h"
#include "chargeshell.h"
#include "gaussianintegration.h"
//#include "simpsonintegration.h"
#include "constantnum.h"

extern double dr;
extern int    DMAX;
extern int    DMAXg;
extern int    ngrid_m; //the number of grids: the middle between two surfaces
extern short  nspecies;

using namespace std;

void ChargeShell(int* LLI,int* ULI,double* BB,double** rho,double*** Psi_IJ,double** DSH)
{
    double    temp1;
    double*   IntegF;
    int**     BIJ;
    short   hspecies;
    int     kin,kip,RanSH,iR_p;
    
    hspecies = nspecies - 1;
    RanSH    = 2*DMAXg;
    
    IntegF  = new double[RanSH+1]();
    BIJ     = new int*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        BIJ[i]  = new int[hspecies]();
    }

    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = BB[i] + BB[j];
            BIJ[i][j] = round(temp1/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(int k=0; k<=ngrid_m; ++k)
        {
            DSH[i][k] = 0;
            for(short j=0; j<hspecies; ++j)
            {
                kin = LLI[j] - k;
                if(kin < (-BIJ[i][j])) kin = -BIJ[i][j];
                kip = ULI[j] - k;
                if(kip > BIJ[i][j]) kip = BIJ[i][j];
                
                kin = kin + DMAXg;
                kip = kip + DMAXg;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    iR_p = k + k1 - DMAXg;
                    IntegF[k1] = rho[j][iR_p];
                    //IntegF[k1] = rho[j][iR_p]*Psi_IJ[i][j][k1];
                }
                temp1 = GaussianIntegrationShell(IntegF,Psi_IJ[i][j],0,RanSH,kin,kip);
                DSH[i][k] += temp1;
            }
        }
    }
    
    delete [] IntegF;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] BIJ[i];
    }
    delete [] BIJ;
}




void ChargeShell(double* rhoB,double* BB,float* Z,double*** Psi_IJ,double* DSHB)
{
    double   temp1;
    double*  IntegF;
    int**    BIJ;
    short   hspecies;
    int     kin,kip,RanSH;
    
    hspecies = nspecies - 1;
    RanSH    = 2*DMAXg;
    
    IntegF  = new double[RanSH+1]();
    BIJ     = new int*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        BIJ[i] = new int[hspecies]();
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = BB[i] + BB[j];
            BIJ[i][j] = round(temp1/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        DSHB[i] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            kin = -BIJ[i][j];
            kip = BIJ[i][j];
            
            kin = kin + DMAXg;
            kip = kip + DMAXg;
            for(int k1=kin; k1<=kip; ++k1)
            {
                IntegF[k1] = rhoB[j];
            }
            
            
            temp1 = GaussianIntegrationShell(IntegF,Psi_IJ[i][j],0,RanSH,kin,kip);
            DSHB[i] += temp1;
        }
        
    }
    delete [] IntegF;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] BIJ[i];
    }
    delete [] BIJ;
}





void ChargeShellDirk(int* LLI,int* ULI,float* D,double** rho,double*** Psi_IJ,double** DSH)
{
    double    temp1;
    double*   IntegF;
    int**     BIJ;
    short   hspecies;
    int     kin,kip,RanSH,iR_p;
    
    hspecies = nspecies - 1;
    RanSH    = 2*DMAX;
    
    IntegF  = new double[RanSH+1]();
    BIJ     = new int*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        BIJ[i]  = new int[hspecies]();
    }

    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = (D[i] + D[j])*0.5;
            BIJ[i][j] = round(temp1/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(int k=0; k<=ngrid_m; ++k)
        {
            DSH[i][k] = 0;
            for(short j=0; j<hspecies; ++j)
            {
                kin = LLI[j] - k;
                if(kin < (-BIJ[i][j])) kin = -BIJ[i][j];
                kip = ULI[j] - k;
                if(kip > BIJ[i][j]) kip = BIJ[i][j];
                
                kin = kin + DMAX;
                kip = kip + DMAX;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    iR_p = k + k1 - DMAX;
                    IntegF[k1] = rho[j][iR_p];
                    //IntegF[k1] = rho[j][iR_p]*Psi_IJ[i][j][k1];
                }
                temp1 = GaussianIntegrationShell(IntegF,Psi_IJ[i][j],0,RanSH,kin,kip);
                DSH[i][k] += temp1;
            }
        }
    }
    
    delete [] IntegF;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] BIJ[i];
    }
    delete [] BIJ;
}




void ChargeShellDirk(double* rhoB,float* D,float* Z,double*** Psi_IJ,double* DSHB)
{
    double   temp1;
    double*  IntegF;
    int**    BIJ;
    short   hspecies;
    int     kin,kip,RanSH;
    
    hspecies = nspecies - 1;
    RanSH    = 2*DMAX;
    
    IntegF  = new double[RanSH+1]();
    BIJ     = new int*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        BIJ[i] = new int[hspecies]();
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = (D[i] + D[j])*0.5;
            BIJ[i][j] = round(temp1/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        DSHB[i] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            kin = -BIJ[i][j];
            kip = BIJ[i][j];
            
            kin = kin + DMAX;
            kip = kip + DMAX;
            for(int k1=kin; k1<=kip; ++k1)
            {
                IntegF[k1] = rhoB[j];
            }
            
            
            temp1 = GaussianIntegrationShell(IntegF,Psi_IJ[i][j],0,RanSH,kin,kip);
            DSHB[i] += temp1;
        }
        
    }
    delete [] IntegF;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] BIJ[i];
    }
    delete [] BIJ;
}


/*
void ChargeShell(int* LLI,int* ULI,float* D,double** rho,double*** Psi_IJ,double** DSH)
{
    double  temp1;
    double*  IntegF;
    int**     BIJ;
    short   hspecies;
    int     kin,kip,RanSH,iR_p;
    
    hspecies = nspecies - 1;
    RanSH    = 2*DMAX;
    
    IntegF  = new double[RanSH+1]();
    BIJ     = new int*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        BIJ[i]  = new int[hspecies]();
    }

    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = (D[i] + D[j])*0.5;
            BIJ[i][j] = round(temp1/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(int k=0; k<=ngrid_m; ++k)
        {
            DSH[i][k] = 0;
            for(short j=0; j<hspecies; ++j)
            {
                kin = LLI[j] - k;
                if(kin < (-BIJ[i][j])) kin = -BIJ[i][j];
                kip = ULI[j] - k;
                if(kip > BIJ[i][j]) kip = BIJ[i][j];
                
                kin = kin + DMAX;
                kip = kip + DMAX;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    iR_p = k + k1 - DMAX;
                    IntegF[k1] = rho[j][iR_p]*Psi_IJ[i][j][k1];
                }
                temp1 = SimpsonIntegration(IntegF,0,RanSH,kin,kip);
                
                DSH[i][k] += temp1;
            }
        }
        
    }
    
    
    delete [] IntegF;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] BIJ[i];
    }
    delete [] BIJ;
}




void ChargeShell(double* rhoB,float* D,float* Z,double*** Psi_IJ,double* DSHB)
{
    double   temp1;
    double*  IntegF;
    int**     BIJ;
    short   hspecies;
    int     kin,kip,RanSH;
    
    hspecies = nspecies - 1;
    RanSH    = 2*DMAX;
    
    IntegF  = new double[RanSH+1]();
    BIJ     = new int*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        BIJ[i] = new int[hspecies]();
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = (D[i] + D[j])*0.5;
            BIJ[i][j] = round(temp1/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        DSHB[i] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            kin = -BIJ[i][j];
            kip = BIJ[i][j];
            
            kin = kin + DMAX;
            kip = kip + DMAX;
            for(int k1=kin; k1<=kip; ++k1)
            {
                IntegF[k1] = rhoB[j]*Psi_IJ[i][j][k1];
            }
            temp1 = SimpsonIntegration(IntegF,0,RanSH,kin,kip);
            
            DSHB[i] += temp1;
        }
        
    }
    
    
    delete [] IntegF;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] BIJ[i];
    }
    delete [] BIJ;
}
*/
