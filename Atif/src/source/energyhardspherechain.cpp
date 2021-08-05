//**calculating free energy for hard sphere chain: no electrostatics**//
#include "clibrary.h"
#include "weighteddensity.h"
#include "simpsonintegration.h"
#include "constantnum.h"
#include "energyhardspherechain.h"

extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern int ngrid;
extern short  nspecies;
//using namespace std;

void EnergyHardSphereChain(int MAXR,int* LLI,int* ULI,double* rhoB,float* D,float* Z,double** ATT,
                          double** rho,double& Energy_HC)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  NV1I;
    double*  NV2I;
    double** sIJ;
    double** Y;
    double** G;
    double** NI;
    double*  f_EN;
    double   N0,N1,N2,N3,NV1,NV2,N31,LN31,AN0M1,AN0M2,rhoBM1,rhoBM2,temp,En_HS,En_CH,DSIJ;
    short    i1,p_max,nblocks,hspecies;
    
    nblocks = nb[0] + nb[1];
    hspecies= nspecies - 1;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NV1I  = new double[hspecies]();
    NV2I  = new double[hspecies]();
    f_EN  = new double[ngrid+1]();



    

    sIJ    = new double*[nblocks]();
    Y      = new double*[nblocks]();
    G      = new double*[nblocks]();
    NI     = new double*[hspecies]();
    
    for(short i=0; i<hspecies; ++i)
    {
        NI[i]  = new double[6]();
    }

    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
        G[i]    = new double[nblocks]();
        sIJ[i]  = new double[nblocks]();
    }
        
    
    for(short i=0; i<nblocks; ++i)
    {
        for(short j=0; j<nblocks; ++j)
        {
            sIJ[i][j] = 0.5*(D[i] + D[j]);
        }
    }
    
    
    rhoBM1 = 0.0;
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1 = rhoBM1 + rhoB[i];
    }
    rhoBM2 = 0.0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        rhoBM2 = rhoBM2 + rhoB[i];
    }
    

    p_max = 0;
    if(rhoBM1 > 1.0e-12) p_max = p_max + 1;
    if(rhoBM2 > 1.0e-12) p_max = p_max + 1;
    
    //big loop
    for(int k=0; k<=ngrid; ++k)
    {
        WeightedDensity(k,LLI,ULI,D,rho,NI);
        for(short i=0; i<hspecies; ++i)
        {
            N0I[i] = NI[i][0];
            N1I[i] = NI[i][1];
            N2I[i] = NI[i][2];
            N3I[i] = NI[i][3];
            NV1I[i]= NI[i][4];
            NV2I[i]= NI[i][5];
        }
        N0 = 0.0;
        N1 = 0.0;
        N2 = 0.0;
        N3 = 0.0;
        NV1= 0.0;
        NV2= 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            N0 = N0 + N0I[i];
            N1 = N1 + N1I[i];
            N2 = N2 + N2I[i];
            N3 = N3 + N3I[i];
            NV1= NV1 + NV1I[i];
            NV2= NV2 + NV2I[i];
            
        }
        
        if(N3 < 1E-20)
        {
            f_EN[k] = 0;
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            N31  = 1 - N3;
            LN31 = log(N31);

            En_HS = (N1*N2-NV1*NV2)/N31 - N0*LN31 + (N2*N2*N2-3*N2*NV2*NV2)*(N3+N31*N31*LN31)/(36*Pi*N3*N3*N31*N31);
            
            //Start: bulk chemical potential for chain connectivity
            //Here Y is log(Y) actually
            for(short p=0; p<p_max; ++p)
            {
                i1 = p*nb[1];
                for(short i=(p*nb[0]); i<(nb[0]+i1); ++i)
                {
                    DSIJ = D[i]*D[i]/sIJ[i][i];
                    G[i][i] = 1/N31 + N2*DSIJ/(4*N31*N31) + N2*N2*DSIJ*DSIJ/(72*N31*N31*N31);
                    Y[i][i] = log(G[i][i]);
                    
                    if(i < (nb[0]+i1-1))
                    {
                        DSIJ = D[i]*D[i+1]/sIJ[i][i+1];
                        G[i][i+1] = 1/N31 + N2*DSIJ/(4*N31*N31) + N2*N2*DSIJ*DSIJ/(72*N31*N31*N31);
                        Y[i][i+1] = log(G[i][i+1]);
                    }
                }
            }
            
            //The energy of excess free energy density for chain connectivity
            En_CH = 0;
            if(rhoBM1 > 1.0e-12)
            {
                AN0M1 = 0.0;
                for(short i=0; i<nb[0]; ++i)
                {
                    AN0M1 = AN0M1 + N0I[i];
                }
                
                temp = 0.0;
                for(short i0=0; i0<nb[0]; ++i0)
                {
                    temp = temp - ATT[i0][i0]*Y[i0][i0];
                    if(i0 < (nb[0]-1)) temp = temp - ATT[i0][i0+1]*Y[i0][i0+1];
                }
                
                En_CH += (temp*AN0M1);
                
            }
            if(rhoBM2 > 1.0e-12)
            {
                AN0M2 = 0.0;
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    AN0M2 = AN0M2 + N0I[i];
                }
                
                temp = 0.0;
                for(short i0=nb[0]; i0<nblocks; ++i0)
                {
                    temp = temp - ATT[i0][i0]*Y[i0][i0];
                    if(i0 < (nblocks-1)) temp = temp - ATT[i0][i0+1]*Y[i0][i0+1];
                }
                
                En_CH += (temp*AN0M2);
            }
            
            f_EN[k] = En_CH + En_HS;
            //End: bulk chemical potential for chain connectivity
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    Energy_HC = SimpsonIntegration(f_EN,0,ngrid,0,ngrid);
    
        
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] NV1I;
    delete [] NV2I;
    delete [] f_EN;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] NI[i];
    }
    delete [] NI;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] Y[i];
        delete [] G[i];
        delete [] sIJ[i];
    }
    delete [] Y;
    delete [] G;
    delete [] sIJ;
}
