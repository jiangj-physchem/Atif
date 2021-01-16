//**calculating inhomogeneous chemical potential for hard sphere chain: no electrostatics**//
#include "clibrary.h"
#include "weighteddensity.h"
#include "simpsonintegration.h"
#include "constantnum.h"
#include "derivhardspherechain.h"

extern double dr;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int DMAX;
extern int LLIM; //the minimum lower limit of intergral
extern short  nspecies;
//using namespace std;

void DerivHardSphereChain(int MAXR,int* LLI,int* ULI,double* rhoB,float* D,float* Z,double** ATT,
                          double** rho,double** DCH)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  RI;
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  NV1I;
    double*  NV2I;
    double*  DDH;
    double*  D2;
    double*  DCC1;
    double*  tempm;
    double*  PiD;
    double** DDC;
    double** sIJ;
    double** Y;
    double** G;
    double** NI;
    double*** DD;
    double*** DY;
    double   N0,N1,N2,N3,NV1,NV2,N31,LN31,AN0M1,AN0M2,rhoBM1,rhoBM2,temp,Pi2;
    short    i1,p_max,nblocks,hspecies;
    
    double   R,rin,rip,DSIJ;
    int      kin,kip;
    
    nblocks = nb[0] + nb[1];
    hspecies= nspecies - 1;
    
    RI    = new double[hspecies]();
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NV1I  = new double[hspecies]();
    NV2I  = new double[hspecies]();
    D2    = new double[hspecies]();
    PiD   = new double[hspecies]();
    DDH   = new double[6]();
    DCC1  = new double[3]();
    tempm = new double[2]();


    

    sIJ    = new double*[nblocks]();
    Y      = new double*[nblocks]();
    G      = new double*[nblocks]();
    NI     = new double*[hspecies]();
    DDC    = new double*[hspecies]();
    
    DD     = new double**[hspecies]();
    
    Pi2    = 2.0*Pi;
    for(short i=0; i<hspecies; ++i)
    {
        DD[i] = new double*[3]();
        DDC[i]= new double[4]();
        for(short j=0; j<3; ++j)
        {
            DD[i][j] = new double[MAXR+1]();
        }
        
        NI[i]  = new double[6]();
        D2[i]  = D[i]*D[i]*0.25;
        RI[i]  = 1.0/(Pi2*D[i]);
        PiD[i] = Pi*D[i];
    }

    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
        G[i]    = new double[nblocks]();
        sIJ[i]  = new double[nblocks]();
    }
    
    DY  = new double**[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        DY[i]  = new double*[nblocks]();
        for(short j=0; j<nblocks; ++j)
        {
            DY[i][j]  = new double[2]();
        }
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
    for(int k=0; k<=MAXR; ++k)
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
            for(short j1=0; j1<hspecies; ++j1)
            {
                for(short j2=0; j2<3; ++j2)
                {
                    DD[j1][j2][k] = 0;
                }
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            N31  = 1 - N3;
            LN31 = log(N31);
            DDH[0] = -LN31;
            DDH[1] = N2/N31;
            DDH[2] = N1/N31 + (N2*N2-NV2*NV2)*(LN31/N3 + 1/(N31*N31))/(12*Pi*N3);
            DDH[3] = N0/N31 + (N1*N2-NV1*NV2)/(N31*N31) - (N2*N2*N2-3*N2*NV2*NV2)*
                     (LN31/(18*Pi*N3*N3*N3) + 1/(36*Pi*N3*N3*N31) + (1-3*N3)/(36*
                     Pi*N3*N3*N31*N31*N31));
            DDH[4] = -(NV2/N31);
            DDH[5] = -(NV1/N31) - (LN31/N3 + 1/(N31*N31))*N2*NV2/(6*Pi*N3);
            
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
                    
                    G[i][i] = G[i][i]*N31*N31;
                    DY[i][i][0] = (DSIJ*0.25 + N2*DSIJ*DSIJ/(36*N31))/G[i][i];
                    DY[i][i][1] = (1 + N2*DSIJ/(2*N31) + N2*N2*DSIJ*DSIJ/(24*N31*N31))/G[i][i];
                    
                    if(i < (nb[0]+i1-1))
                    {
                        DSIJ = D[i]*D[i+1]/sIJ[i][i+1];
                        G[i][i+1] = 1/N31 + N2*DSIJ/(4*N31*N31) + N2*N2*DSIJ*DSIJ/(72*N31*N31*N31);
                        Y[i][i+1] = log(G[i][i+1]);
                        
                        G[i][i+1] = G[i][i+1]*N31*N31;
                        DY[i][i+1][0] = (DSIJ*0.25 + N2*DSIJ*DSIJ/(36*N31))/G[i][i+1];
                        DY[i][i+1][1] = (1 + N2*DSIJ/(2*N31) + N2*N2*DSIJ*DSIJ/(24*N31*N31))/G[i][i+1];
                    }
                }
            }
            
            //The derivation of excess free energy density for chain connectivity
            for(short i=0; i<hspecies; ++i)
            {
                for(short j=0; j<4; ++j)
                {
                    DDC[i][j] = 0;
                }
            }
            
            if(rhoBM1 > 1.0e-12)
            {
                AN0M1 = 0.0;
                for(short i=0; i<nb[0]; ++i)
                {
                    AN0M1 = AN0M1 + N0I[i];
                }
                
                for(short j=0; j<2; ++j)
                {
                    tempm[j] = 0.0;
                    for(short i0=0; i0<nb[0]; ++i0)
                    {
                        tempm[j] = tempm[j] - ATT[i0][i0]*DY[i0][i0][j];
                        if(i0 < (nb[0]-1)) tempm[j] = tempm[j] - ATT[i0][i0+1]*DY[i0][i0+1][j];
                    }
                    //DDC[i][j] = DDC[i][j] - (N0I[0]+N0I[1])*(ATT[0][0]*DY[0][0][j-2] + ATT[0][1]*DY[0][1][j-2] + ATT[1][1]*DY[1][1][j-2]);
                }
                for(short i=0; i<nb[0]; ++i)
                {
                    DDC[i][2]  += AN0M1*tempm[0];
                    DDC[i][3]  += AN0M1*tempm[1];
                }
                
                for(short i=nblocks; i<hspecies; ++i)
                {
                    DDC[i][2]  += AN0M1*tempm[0];
                    DDC[i][3]  += AN0M1*tempm[1];
                }
                
                temp = 0.0;
                for(short i0=0; i0<nb[0]; ++i0)
                {
                    temp = temp - ATT[i0][i0]*Y[i0][i0];
                    if(i0 < (nb[0]-1)) temp = temp - ATT[i0][i0+1]*Y[i0][i0+1];
                }
                
                for(short i=0; i<nb[0]; ++i)
                {
                    DDC[i][0] += temp;
                    //DDC[i][0] = DDC[i][0] - (ATT[0][0]*Y[0][0] + ATT[0][1]*Y[0][1] + ATT[1][1]*Y[1][1]);
                }
            }
            if(rhoBM2 > 1.0e-12)
            {
                AN0M2 = 0.0;
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    AN0M2 = AN0M2 + N0I[i];
                }
                
                for(short j=0; j<2; ++j)
                {
                    tempm[j] = 0.0;
                    for(short i0=nb[0]; i0<nblocks; ++i0)
                    {
                        tempm[j] = tempm[j] - ATT[i0][i0]*DY[i0][i0][j];
                        if(i0 < (nblocks-1)) tempm[j] = tempm[j] - ATT[i0][i0+1]*DY[i0][i0+1][j];
                    }
                    //DDC[i][j] = DDC[i][j] - (N0I[2]+N0I[3])*(ATT[2][2]*DY[2][2][j-2] + ATT[2][3]*DY[2][3][j-2] + ATT[3][3]*DY[3][3][j-2]);
                }
                
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    DDC[i][2] += AN0M2*tempm[0];
                    DDC[i][3] += AN0M2*tempm[1];
                }
                
                for(short i=nblocks; i<hspecies; ++i)
                {
                    DDC[i][2] += AN0M2*tempm[0];
                    DDC[i][3] += AN0M2*tempm[1];
                }
                
                temp = 0.0;
                for(short i0=nb[0]; i0<nblocks; ++i0)
                {
                    temp = temp - ATT[i0][i0]*Y[i0][i0];
                    if(i0 < (nblocks-1)) temp = temp - ATT[i0][i0+1]*Y[i0][i0+1];
                }
                
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    DDC[i][0] += temp;
                    //DDC[i][0] = DDC[i][0] - (ATT[2][2]*Y[2][2] + ATT[2][3]*Y[2][3] + ATT[3][3]*Y[3][3]);
                }
            }

            for(short i=0; i<hspecies; ++i)
            {
                DD[i][0][k] = DDH[2] + DDC[i][2] + DDH[1]*RI[i] + 2.0*RI[i]*(DDH[0]+DDC[i][0])/D[i];
                DD[i][1][k] = DDH[3] + DDC[i][3];
                DD[i][2][k] = DDH[5] + DDH[4]*RI[i];
            }
            
            
            //End: bulk chemical potential for chain connectivity
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=LLIM; i<=ngrid_m; ++i)
    {
        R = i*dr;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - D[j]*0.5;
            rip = R + D[j]*0.5;
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > MAXR) kip = MAXR;
            
            for(short j1=0; j1<3; ++j1)
            {
                DCC1[j1] = 0;
            }
                        
            DCC1[0] = SimpsonIntegration(DD[j][0],0,MAXR,kin,kip);
            DCC1[1] = SimpsonIntegration(DD[j][1],0,MAXR,kin,kip,D2[j],i);
            DCC1[2] = SimpsonIntegration(DD[j][2],0,MAXR,kin,kip,i);
            
            
            DCC1[0] = DCC1[0]*PiD[j];
            DCC1[1] = DCC1[1]*Pi;
            DCC1[2] = DCC1[2]*Pi2;
            
            DCH[j][i] = 0;
            for(short j1=0; j1<3; ++j1)
            {
                DCH[j][i] = DCH[j][i] + DCC1[j1];
            }
        }
    }
        
    delete [] RI;
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] NV1I;
    delete [] NV2I;
    delete [] DCC1;
    delete [] PiD;
    delete [] D2;
    delete [] DDH;
    delete [] tempm;
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<3; ++j)
        {
            delete [] DD[i][j];
        }
        delete [] DD[i];
        delete [] NI[i];
        delete [] DDC[i];
    }
    delete [] DD;
    delete [] NI;
    delete [] DDC;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] Y[i];
        delete [] G[i];
        delete [] sIJ[i];
    }
    delete [] Y;
    delete [] G;
    delete [] sIJ;

    for(short i=0; i<nblocks; ++i)
    {
        for(short j=0; j<nblocks; ++j)
        {
            delete [] DY[i][j];
        }
        delete [] DY[i];
    }
    delete [] DY;

}




void DerivHardSphereChain(double* rhoB,float* D,float* Z,double** ATT,double* DCHB,double& HCH_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  RI;
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  D2;
    double*  tempm;
    double*  PiD;
    double*  DDB1;
    double*  DDB2;
    double** DDC;
    double** sIJ;
    double** Y;
    double** G;
    double** NI;
    double** DD1;
    double** DD2;
    double*** DY;
    double   N0,N1,N2,N3,N31,LN31,AN0M1,AN0M2,rhoBM1,rhoBM2,temp,Pi2,DCC1,DCC2;
    short    i1,p_max,nblocks,hspecies;
    
    double   rin,rip,DSIJ,f_ch;
    int      kin,kip,RMAX;
    
    nblocks = nb[0] + nb[1];
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    RI    = new double[hspecies]();
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    D2    = new double[hspecies]();
    PiD   = new double[hspecies]();
    DDH   = new double[4]();
    tempm = new double[2]();


    

    sIJ    = new double*[nblocks]();
    Y      = new double*[nblocks]();
    G      = new double*[nblocks]();
    NI     = new double*[hspecies]();
    DDC    = new double*[hspecies]();
    
    DD1    = new double*[hspecies]();
    DD2    = new double*[hspecies]();
    DDB1   = new double[hspecies]();
    DDB2   = new double[hspecies]();
    
    Pi2    = 2.0*Pi;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i]= new double[DMAX+1]();
        DD2[i]= new double[DMAX+1]();
        DDC[i]= new double[4]();
        
        NI[i]  = new double[4]();
        D2[i]  = D[i]*D[i]*0.25;
        RI[i]  = 1.0/(Pi2*D[i]);
        PiD[i] = Pi*D[i];
    }

    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
        G[i]    = new double[nblocks]();
        sIJ[i]  = new double[nblocks]();
    }
    
    DY  = new double**[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        DY[i]  = new double*[nblocks]();
        for(short j=0; j<nblocks; ++j)
        {
            DY[i][j]  = new double[2]();
        }
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
    
    WeightedDensity(D,rhoB,NI);
    for(short i=0; i<hspecies; ++i)
    {
        N0I[i] = NI[i][0];
        N1I[i] = NI[i][1];
        N2I[i] = NI[i][2];
        N3I[i] = NI[i][3];
    }
    N0 = 0.0;
    N1 = 0.0;
    N2 = 0.0;
    N3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        N0 = N0 + N0I[i];
        N1 = N1 + N1I[i];
        N2 = N2 + N2I[i];
        N3 = N3 + N3I[i];
    }
    
    
    HCH_EN = 0;
    
    if(N3 < 1E-20)
    {
        for(short j1=0; j1<hspecies; ++j1)
        {
            DDB1[j1] = 0;
            DDB2[j1] = 0;
        }
    }
    else //small loop
    {
        //Start: derivation of the hard sphere contribution
        if(N3 > 0.99) N3 = 0.99;
        
        N31  = 1 - N3;
        LN31 = log(N31);
        DDH[0] = -LN31;
        DDH[1] = N2/N31;
        DDH[2] = N1/N31 + (N2*N2)*(LN31/N3 + 1/(N31*N31))/(12*Pi*N3);
        DDH[3] = N0/N31 + (N1*N2)/(N31*N31) - (N2*N2*N2)*(LN31/(18*Pi*N3*N3*N3) + 1/(36*Pi*N3*N3*N31) + (1-3*N3)/(36*
                 Pi*N3*N3*N31*N31*N31));
        
        
        HCH_EN = N1*N2/N31 - N0*LN31 + N2*N2*N2*(N3+N31*N31*LN31)/(36*Pi*N3*N3*N31*N31);
        
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
                
                G[i][i] = G[i][i]*N31*N31;
                DY[i][i][0] = (DSIJ*0.25 + N2*DSIJ*DSIJ/(36*N31))/G[i][i];
                DY[i][i][1] = (1 + N2*DSIJ/(2*N31) + N2*N2*DSIJ*DSIJ/(24*N31*N31))/G[i][i];
                
                if(i < (nb[0]+i1-1))
                {
                    DSIJ = D[i]*D[i+1]/sIJ[i][i+1];
                    G[i][i+1] = 1/N31 + N2*DSIJ/(4*N31*N31) + N2*N2*DSIJ*DSIJ/(72*N31*N31*N31);
                    Y[i][i+1] = log(G[i][i+1]);
                    
                    G[i][i+1] = G[i][i+1]*N31*N31;
                    DY[i][i+1][0] = (DSIJ*0.25 + N2*DSIJ*DSIJ/(36*N31))/G[i][i+1];
                    DY[i][i+1][1] = (1 + N2*DSIJ/(2*N31) + N2*N2*DSIJ*DSIJ/(24*N31*N31))/G[i][i+1];
                }
            }
        }
        
        //The derivation of excess free energy density for chain connectivity
        for(short i=0; i<hspecies; ++i)
        {
            for(short j=0; j<4; ++j)
            {
                DDC[i][j] = 0;
            }
        }
        
        f_ch = 0;
        if(rhoBM1 > 1.0e-12)
        {
            AN0M1 = 0.0;
            for(short i=0; i<nb[0]; ++i)
            {
                AN0M1 = AN0M1 + N0I[i];
            }
            
            for(short j=0; j<2; ++j)
            {
                tempm[j] = 0.0;
                for(short i0=0; i0<nb[0]; ++i0)
                {
                    tempm[j] = tempm[j] - ATT[i0][i0]*DY[i0][i0][j];
                    if(i0 < (nb[0]-1)) tempm[j] = tempm[j] - ATT[i0][i0+1]*DY[i0][i0+1][j];
                }
                //DDC[i][j] = DDC[i][j] - (N0I[0]+N0I[1])*(ATT[0][0]*DY[0][0][j-2] + ATT[0][1]*DY[0][1][j-2] + ATT[1][1]*DY[1][1][j-2]);
            }
            for(short i=0; i<nb[0]; ++i)
            {
                DDC[i][2]  += AN0M1*tempm[0];
                DDC[i][3]  += AN0M1*tempm[1];
            }
            
            for(short i=nblocks; i<hspecies; ++i)
            {
                DDC[i][2]  += AN0M1*tempm[0];
                DDC[i][3]  += AN0M1*tempm[1];
            }
            
            temp = 0.0;
            for(short i0=0; i0<nb[0]; ++i0)
            {
                temp = temp - ATT[i0][i0]*Y[i0][i0];
                if(i0 < (nb[0]-1)) temp = temp - ATT[i0][i0+1]*Y[i0][i0+1];
            }
            
            for(short i=0; i<nb[0]; ++i)
            {
                DDC[i][0] += temp;
                f_ch      += N0I[i]*temp;
                //DDC[i][0] = DDC[i][0] - (ATT[0][0]*Y[0][0] + ATT[0][1]*Y[0][1] + ATT[1][1]*Y[1][1]);
            }
        }
        if(rhoBM2 > 1.0e-12)
        {
            AN0M2 = 0.0;
            for(short i=nb[0]; i<nblocks; ++i)
            {
                AN0M2 = AN0M2 + N0I[i];
            }
            
            for(short j=0; j<2; ++j)
            {
                tempm[j] = 0.0;
                for(short i0=nb[0]; i0<nblocks; ++i0)
                {
                    tempm[j] = tempm[j] - ATT[i0][i0]*DY[i0][i0][j];
                    if(i0 < (nblocks-1)) tempm[j] = tempm[j] - ATT[i0][i0+1]*DY[i0][i0+1][j];
                }
                //DDC[i][j] = DDC[i][j] - (N0I[2]+N0I[3])*(ATT[2][2]*DY[2][2][j-2] + ATT[2][3]*DY[2][3][j-2] + ATT[3][3]*DY[3][3][j-2]);
            }
            
            for(short i=nb[0]; i<nblocks; ++i)
            {
                DDC[i][2] += AN0M2*tempm[0];
                DDC[i][3] += AN0M2*tempm[1];
            }
            
            for(short i=nblocks; i<hspecies; ++i)
            {
                DDC[i][2] += AN0M2*tempm[0];
                DDC[i][3] += AN0M2*tempm[1];
            }
            
            temp = 0.0;
            for(short i0=nb[0]; i0<nblocks; ++i0)
            {
                temp = temp - ATT[i0][i0]*Y[i0][i0];
                if(i0 < (nblocks-1)) temp = temp - ATT[i0][i0+1]*Y[i0][i0+1];
            }
            
            for(short i=nb[0]; i<nblocks; ++i)
            {
                DDC[i][0] += temp;
                f_ch      += N0I[i]*temp;
                //DDC[i][0] = DDC[i][0] - (ATT[2][2]*Y[2][2] + ATT[2][3]*Y[2][3] + ATT[3][3]*Y[3][3]);
            }
        }

        for(short i=0; i<hspecies; ++i)
        {
            DDB1[i] = DDH[2] + DDC[i][2] + DDH[1]*RI[i] + 2.0*RI[i]*(DDH[0]+DDC[i][0])/D[i];
            DDB2[i] = DDH[3] + DDC[i][3];
        }
        
        HCH_EN = HCH_EN + f_ch;
        //End: bulk chemical potential for chain connectivity
        
    }//small loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=DMAX; ++k)
        {
            DD1[j][k] = DDB1[j];
            DD2[j][k] = DDB2[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = RMAX*dr - D[j]*0.5;
        rip = RMAX*dr + D[j]*0.5;
        kin = round(rin/dr);
        kip = round(rip/dr);
        

                                
        DCC1 = SimpsonIntegration(DD1[j],0,DMAX,kin,kip);
        DCC2 = SimpsonIntegration(DD2[j],0,DMAX,kin,kip,D2[j],RMAX);
        
        
        DCC1 = DCC1*PiD[j];
        DCC2 = DCC2*Pi;
        
        DCHB[j] = DCC1 + DCC2;
    }
            
    delete [] RI;
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] PiD;
    delete [] D2;
    delete [] DDH;
    delete [] tempm;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] DD2[i];
        delete [] NI[i];
        delete [] DDC[i];
    }
    delete [] DD1;
    delete [] DD2;
    delete [] DDB1;
    delete [] DDB2;
    delete [] NI;
    delete [] DDC;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] Y[i];
        delete [] G[i];
        delete [] sIJ[i];
    }
    delete [] Y;
    delete [] G;
    delete [] sIJ;

    for(short i=0; i<nblocks; ++i)
    {
        for(short j=0; j<nblocks; ++j)
        {
            delete [] DY[i][j];
        }
        delete [] DY[i];
    }
    delete [] DY;

}
