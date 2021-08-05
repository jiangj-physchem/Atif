//*******************solve the coefficience of the propagator ***************//
#include "renormGLGR.h"
#include "simpsonintegration.h"
#include "constantnum.h"

extern double dr;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block

extern double* DL1;
extern double* DR1;
extern double* DL2;
extern double* DR2;
extern double* GRC1;
extern double* GRC2;

void RenormGLGR(float* D,double* rhoB,double*** BesselZero,string* MODEL)
{
    short*   MB1;
    short*   MB2;
    double*  f_1;
    
    
    double rhoBM1,rhoBM2; //,coe
    
    short  mp1,mp2,nblocks,i0,j0;
    int    nBessel,nBessel0;

    nblocks = nb[0] + nb[1];
    
    nBessel0= 2*round(1.0/dr);
    nBessel = nBessel0 + 1;

    
    mp1 = mp[0];
    mp2 = mp[1];
    
    MB1  = new short[nb[0]+1]();
    MB2  = new short[nb[1]+1]();
    f_1  = new double[nBessel]();
    
    rhoBM1 = 0.0;
    MB1[0] = 0;
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1   = rhoBM1 + rhoB[i];
        MB1[i+1] = MB1[i] + mb[0][i];
    }
    rhoBM2 = 0.0;
    MB2[0] = 0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        j0        = i - nb[0];
        rhoBM2    = rhoBM2 + rhoB[i];
        MB2[j0+1] = MB2[j0] + mb[1][j0];
    }
    

    //loop-2
    if(rhoBM1 > 1E-10)
    {
        for(short i=0; i<nb[0]; ++i)
        {
            for(short j=(MB1[i]+1); j<(MB1[i+1]-1); ++j)
            {
                DL1[j] = D[i];
                DR1[j] = D[i];
            }
            
            j0 = MB1[i];
            if(j0 < mp1)
            {
                DR1[j0] = D[i];
                if((j0+1) == MB1[i+1]) DR1[j0] = 0.5*(D[i] + D[i+1]);
                if(j0 > 0) DL1[j0] = 0.5*(D[i] + D[i-1]);
            }
            
            j0      = MB1[i+1] - 1;
            DL1[j0] = D[i];
            if((j0 == MB1[i]) && (j0 > 0)) DL1[j0] = 0.5*(D[i] + D[i-1]);
            if(j0 < (mp1-1)) DR1[j0] = 0.5*(D[i] + D[i+1]);
        }
        
        
        
        
        if(MODEL[0] == "semiflexible" || MODEL[0] == "semi-flexible")//stiff chain
        {
            for(int k1=0; k1<nBessel; ++k1)
            {
                for(int k2=0; k2<nBessel; ++k2)
                {
                    f_1[k2]  = BesselZero[0][k2][k1];
                }
                GRC1[k1] = SimpsonIntegration(f_1,0,nBessel0,0,nBessel0);
            }
        }//stiff chain
    }
    //loop-2
    
    //loop-3
    if(rhoBM2 > 1E-10)
    {
        for(short i=0; i<nb[1]; ++i)
        {
            i0 = i + nb[0];
            for(short j=(MB2[i]+1); j<(MB2[i+1]-1); ++j)
            {
                DL2[j] = D[i0];
                DR2[j] = D[i0];
            }
            
            j0 = MB2[i];
            if(j0 < mp2)
            {
                DR2[j0] = D[i0];
                if((j0+1) == MB2[i+1]) DR2[j0] = 0.5*(D[i0] + D[i0+1]);
                if(j0 > 0) DL2[j0] = 0.5*(D[i0] + D[i0-1]);
            }
            j0      = MB2[i+1] - 1;
            DL2[j0] = D[i0];
            if((j0 == MB2[i]) && (j0 > 0)) DL2[j0] = 0.5*(D[i0] + D[i0-1]);
            if(j0 < (mp2-1)) DR2[j0] = 0.5*(D[i0] + D[i0+1]);
        }
        
        if(MODEL[1] == "semiflexible" || MODEL[1] == "semi-flexible")//stiff chain
        {
            for(int k1=0; k1<nBessel; ++k1)
            {//loop for Z_2
                for(int k2=0; k2<nBessel; ++k2)
                {//loop for Z_0
                    f_1[k2]  = BesselZero[1][k2][k1];
                }//loop for Z_0
                GRC2[k1] = SimpsonIntegration(f_1,0,nBessel0,0,nBessel0);

            }//loop for Z_2
        }//stiff chain
    }
    //loop-3
    
    
    delete [] MB1;
    delete [] MB2;
    delete [] f_1;
}
