//***********calculating bulk chemical potential****************//
#include "clibrary.h"
#include "bulkchempotential.h"
#include "constantnum.h"
extern double dr;
extern double depthW;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern short  nspecies;

void BulkChemPotential(double* rhoB,float* D,float* Z,double** pairEner,
                       double* mu)
{
    double*  muAT;
    double*  muID;
    double** sIJ;
    double   lGau,rGau,rin1,R1,tempv,rhoBM1,rhoBM2,VS,VI;
    double   rin2,rin3,rin4,SWD,SWD2,SIJ2;
    double   pEner,coe,temp;
    short    nblocks,hspecies;
    int      kin1,kin2,kin3,kin4;
    
    
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr

    
    nblocks = nb[0] + nb[1];
    hspecies= nspecies - 1;
    muAT   = new double[hspecies]();
    muID   = new double[hspecies]();
    sIJ    = new double*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        sIJ[i] = new double[hspecies]();
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            sIJ[i][j] = 0.5*(D[i] + D[j]);
        }
    }
    
    
    rhoBM1= 0.0;
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1 = rhoBM1 + rhoB[i];
    }
    rhoBM2= 0.0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        rhoBM2 = rhoBM2 + rhoB[i];
    }
    

    
    //Start: bulk chemical potential for ideal gas
    for(short i=0; i<hspecies; ++i)
    {
        muID[i] = -1E20;
        if((i < nb[0]) && (rhoBM1 >  1E-15)) //if(i==0 || i==1)
        {
            muID[i] = log(rhoBM1/((double) mp[0]))/((double) mp[0]);
        }
        else if((i >= nb[0]) && (i < nblocks) && (rhoBM2 >  1E-15)) //else if(i==2 || i==3)
        {
            muID[i] = log(rhoBM2/((double) mp[1]))/((double) mp[1]);
        }
        else if((i >= nblocks) && (rhoB[i] >  1E-15))
        {
            muID[i] = log(rhoB[i]);
        }
    }
    
    //van del Waals
    //depthW= 1.2;
    coe  = Pi*dr*0.5;
    for(short j1=0; j1<hspecies; ++j1)
    {
        muAT[j1] = 0;
        for(short j2=0; j2<hspecies; ++j2)
        {
            tempv = 0;
            pEner = pairEner[j1][j2];
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
                
                
                //integration from - sIJ to sIJ
                temp = 0;
                SWD2  = depthW*depthW;
                SIJ2  = sIJ[j1][j2]*sIJ[j1][j2];
                for(int k=kin2; k<kin3; ++k)
                {
                    //R1  = (k + lGau)*dr;
                    temp = temp + rhoB[j2];
                    
                    //R1  = (k + rGau)*dr;
                    temp = temp + rhoB[j2];
                }
                temp = temp*SIJ2*(SWD2 - 1);
                tempv = tempv + temp;
                
                
                //integration from sIJ to SWD
                temp = 0;
                SWD2  = SWD*SWD;
                for(int k=kin3; k<kin4; ++k)
                {
                    R1  = (k + lGau)*dr;
                    temp = temp + rhoB[j2]*(SWD2 - R1*R1);
                    
                    R1  = (k + rGau)*dr;
                    temp = temp + rhoB[j2]*(SWD2 - R1*R1);
                }
                tempv = tempv + temp;
                
                //integration from - SWD to - sIJ
                temp = 0;
                for(int k=kin1; k<kin2; ++k)
                {
                    R1  = (k + lGau)*dr;
                    temp = temp + rhoB[j2]*(SWD2 - R1*R1);
                    
                    R1  = (k + rGau)*dr;
                    temp = temp + rhoB[j2]*(SWD2 - R1*R1);
                }
                tempv = tempv + temp;
            }
            muAT[j1] = muAT[j1] + tempv*pEner;

        }
        muAT[j1] = muAT[j1]*coe;
    }
    
    
    VS = D[hspecies]*D[hspecies]*D[hspecies];
    for(short i=0; i<hspecies; ++i)
    {
        VI = D[i]*D[i]*D[i];
        mu[i] = muID[i] + muAT[i]-(VI/VS)*log(rhoB[hspecies]);
    }
    
    

    delete [] muID;
    delete [] muAT;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] sIJ[i];
    }
    delete [] sIJ;
}
