//*********** calculating van del Waal ****************//
#include "clibrary.h"
#include "inhomvandelwaal.h"
#include "rombergintegration.h"
#include "constantnum.h"

extern double dr;
extern double size;
extern double BJ;
extern double depthW;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int LLIM; //the minimum lower limit of intergral
extern short nspecies;

void InhomVanDelWaal(int* LLI,int* ULI,float* D,double** rho,double** pairEner,
                     double** VanDW)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double R,R1,lGau,rGau,coe,rhoD;
    double rin1,rin2,rin3,rin4,temp0,SWD2,SWD,SIJ2,pEner,tempv;
    double** sIJ;
    int kin1,kin2,kin3,kin4;
    short   hspecies;
    
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    coe = Pi*dr*0.5;
    hspecies = nspecies - 1;
    
    sIJ    = new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        sIJ[j]  = new double[hspecies]();
    }
    
    for(short j1=0; j1<hspecies; ++j1)
    {
        for(short j2=0; j2<hspecies; ++j2)
        {
            sIJ[j1][j2] = 0.5*(D[j1]+ D[j2]);
        }
    }
    
    for(int i=LLIM; i<=ngrid_m; ++i)
    {
        R   = i*dr;
        //depthW= 1.2;
        for(short j1=0; j1<hspecies; ++j1)
        {
            VanDW[j1][i] = 0;
            for(short j2=0; j2<hspecies; ++j2)
            {
                pEner = pairEner[j1][j2];
                tempv = 0;
                if(pEner != 0)
                {
                    SWD  = depthW*sIJ[j1][j2];
                    rin1 = R - SWD;
                    rin2 = R - sIJ[j1][j2];
                    rin3 = R + sIJ[j1][j2];
                    rin4 = R + SWD;
                    kin1 = round(rin1/dr);
                    kin2 = round(rin2/dr);
                    kin3 = round(rin3/dr);
                    kin4 = round(rin4/dr);
                    
                    if(kin1 < LLI[j2]) kin1 = LLI[j2];
                    if(kin2 < LLI[j2]) kin2 = LLI[j2];
                    if(kin2 > ULI[j2]) kin2 = ULI[j2];
                    if(kin3 < LLI[j2]) kin3 = LLI[j2];
                    if(kin3 > ULI[j2]) kin3 = ULI[j2];
                    if(kin4 > ULI[j2]) kin4 = ULI[j2];
                    
                    
                    //integration from z - sIJ to z + sIJ
                    temp0 = 0;
                    SWD2  = depthW*depthW;
                    SIJ2  = sIJ[j1][j2]*sIJ[j1][j2];
                    for(int k=kin2; k<kin3; ++k)
                    {
                        //R1  = (k + lGau)*dr;
                        rhoD= rho[j2][k] + lGau*(rho[j2][k+1] - rho[j2][k]);
                        temp0 = temp0 + rhoD;
                        
                        //R1  = (k + rGau)*dr;
                        rhoD= rho[j2][k] + rGau*(rho[j2][k+1] - rho[j2][k]);
                        temp0 = temp0 + rhoD;
                    }
                    temp0 = temp0*SIJ2*(SWD2 - 1);
                    tempv = tempv + temp0;
                    
                    
                    //integration from z + sIJ to z + SWD
                    temp0 = 0;
                    SWD2  = SWD*SWD;
                    for(int k=kin3; k<kin4; ++k)
                    {
                        R1  = (k + lGau)*dr;
                        rhoD= rho[j2][k] + lGau*(rho[j2][k+1] - rho[j2][k]);
                        temp0 = temp0 + rhoD*(SWD2 - (R1-R)*(R1-R));
                        
                        R1  = (k + rGau)*dr;
                        rhoD= rho[j2][k] + rGau*(rho[j2][k+1] - rho[j2][k]);
                        temp0 = temp0 + rhoD*(SWD2 - (R1-R)*(R1-R));
                    }
                    tempv = tempv + temp0;
                    
                    
                    //integration from z - SWD to z - sIJ
                    temp0 = 0;
                    for(int k=kin1; k<kin2; ++k)
                    {
                        R1  = (k + lGau)*dr;
                        rhoD= rho[j2][k] + lGau*(rho[j2][k+1] - rho[j2][k]);
                        temp0 = temp0 + rhoD*(SWD2 - (R1-R)*(R1-R));
                        
                        R1  = (k + rGau)*dr;
                        rhoD= rho[j2][k] + rGau*(rho[j2][k+1] - rho[j2][k]);
                        temp0 = temp0 + rhoD*(SWD2 - (R1-R)*(R1-R));
                    }
                    tempv = tempv + temp0;
                }
                VanDW[j1][i] = VanDW[j1][i] + tempv*pEner;
                
                
            }
            VanDW[j1][i] = VanDW[j1][i]*coe;
        }
        
    }
    
    
    for(short j=0; j<hspecies; ++j)
    {
        delete [] sIJ[j];
    }
    delete [] sIJ;
}


void InhomVanDelWaal(double* rhoB,float* D,double** pairEner,double* VanDWB,double& Van_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double R1,lGau,rGau,coe;
    double rin1,rin2,rin3,rin4,temp0,SWD2,SWD,SIJ2,pEner,tempv;
    double** sIJ;
    int kin1,kin2,kin3,kin4;
    short   hspecies;
    
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    coe = Pi*dr*0.5;
    hspecies = nspecies - 1;
    
    sIJ    = new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        sIJ[j]  = new double[hspecies]();
    }
    
    for(short j1=0; j1<hspecies; ++j1)
    {
        for(short j2=0; j2<hspecies; ++j2)
        {
            sIJ[j1][j2] = 0.5*(D[j1]+ D[j2]);
        }
    }
    
    
    Van_EN = 0;
    for(short j1=0; j1<hspecies; ++j1)
    {
        VanDWB[j1] = 0;
        for(short j2=0; j2<hspecies; ++j2)
        {
            pEner = pairEner[j1][j2];
            tempv = 0;
            if(pEner != 0)
            {
                SWD  = depthW*sIJ[j1][j2];
                rin1 = -SWD;
                rin2 = -sIJ[j1][j2];
                rin3 = sIJ[j1][j2];
                rin4 = SWD;
                kin1 = round(rin1/dr);
                kin2 = round(rin2/dr);
                kin3 = round(rin3/dr);
                kin4 = round(rin4/dr);
                
                
                //integration from z - sIJ to z + sIJ
                temp0 = 0;
                SWD2  = depthW*depthW;
                SIJ2  = sIJ[j1][j2]*sIJ[j1][j2];
                for(int k=kin2; k<kin3; ++k)
                {
                    //R1  = (k + lGau)*dr;
                    temp0 = temp0 + rhoB[j2];
                    
                    //R1  = (k + rGau)*dr;
                    temp0 = temp0 + rhoB[j2];
                }
                temp0 = temp0*SIJ2*(SWD2 - 1);
                tempv = tempv + temp0;
                
                
                //integration from z + sIJ to z + SWD
                temp0 = 0;
                SWD2  = SWD*SWD;
                for(int k=kin3; k<kin4; ++k)
                {
                    R1  = (k + lGau)*dr;
                    temp0 = temp0 + rhoB[j2]*(SWD2 - R1*R1);
                    
                    R1  = (k + rGau)*dr;
                    temp0 = temp0 + rhoB[j2]*(SWD2 - R1*R1);
                }
                tempv = tempv + temp0;
                
                
                //integration from z - SWD to z - sIJ
                temp0 = 0;
                for(int k=kin1; k<kin2; ++k)
                {
                    R1  = (k + lGau)*dr;
                    temp0 = temp0 + rhoB[j2]*(SWD2 - R1*R1);
                    
                    R1  = (k + rGau)*dr;
                    temp0 = temp0 + rhoB[j2]*(SWD2 - R1*R1);
                }
                tempv = tempv + temp0;
            }
            VanDWB[j1] = VanDWB[j1] + tempv*pEner;
            
            
        }
        VanDWB[j1] = VanDWB[j1]*coe;
        
        Van_EN    += (VanDWB[j1]*rhoB[j1]*0.5);
    }
    
    
    for(short j=0; j<hspecies; ++j)
    {
        delete [] sIJ[j];
    }
    delete [] sIJ;
}
