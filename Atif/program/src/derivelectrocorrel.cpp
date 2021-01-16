//**calculating inhomogeneous chemical potential for electrostatic correlation**//
#include "clibrary.h"
#include "weighteddensity.h"
#include "simpsonintegration.h"
#include "constantnum.h"
#include "calculategama.h"
#include "derivelectrocorrel.h"

extern double dr;
extern double BJ;
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int DMAX;
extern int LLIM; //the minimum lower limit of intergral
extern short  nspecies;
//using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// 1, R, S, V ///////////////////////////////////////////////////
void DerivElectroCorrel(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                          double** rho,double** DES)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  RI;
    double*  RI2;
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  B2;
    double*  PiB;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double** H;
    double** NI;
    double** DD1;
    double** DD2;
    double   temp1,temp2,temp3,N3,Pi2,Pi4,gammar,gamma_k,BJGA,Chi,gamma_s;
    short    hspecies;
    
    double   R,rin,rip;
    int      kin,kip,igammar,k1,gMAXR;
    
    hspecies= nspecies - 1;
    
    RI    = new double[hspecies]();
    RI2   = new double[hspecies]();
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B2    = new double[hspecies]();
    PiB   = new double[hspecies]();
    B     = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    DDH   = new double[4]();


    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    DD2    = new double*[hspecies]();
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    gMAXR  = igammar + MAXR;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[gMAXR+1]();
        DD2[i] = new double[gMAXR+1]();
        
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i]/B[i];
        H[0][i] = 1.0;
        H[1][i] = temp1;
        H[2][i] = temp1*temp1;
        H[3][i] = H[2][i]*temp1;
        /////////////////////////////////////
        
        B2[i]  = B[i]*B[i];
        RI[i]  = 1.0/(Pi4*B[i]);
        RI2[i] = RI[i]/B[i];
        PiB[i] = Pi2*B[i];
    }
    
    //big loop
    for(int k=0; k<=gMAXR; ++k)
    {
        k1 = k - igammar;
        WeightedDensity(k1,LLI,ULI,B,H,rho,NI);
        gamma_k=CalculateGama(Z,D,NI);
        
        N3 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            N0I[i] = NI[i][0];
            N1I[i] = NI[i][1];
            N2I[i] = NI[i][2];
            N3I[i] = NI[i][3];
            
            N3 += N3I[i];
        }
        
        BJGA  = BJ*gamma_k;
        
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gamma_k*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            
            temp1 += (Z[i]*N1I[i]*GammaS1[i]);
            temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
            temp3 += (N3I[i]*GammaS1[i]);
        }
        Chi     = 1.0/(1-N3+3*temp3);
        gamma_s = temp2*Chi;
        
        if(N3 < 1E-20)
        {
            for(short j1=0; j1<hspecies; ++j1)
            {
                DD1[j1][k] = 0;
                DD2[j1][k] = 0;
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;

            for(short i=0; i<hspecies; ++i)
            {
                DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
                DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
                DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
                DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
                
                
                DD1[i][k] = BJGA*(DDH[2]*H[2][i] + DDH[1]*RI[i]*H[1][i] + RI2[i]*DDH[0]);
                DD2[i][k] = BJGA*DDH[3];
            }
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=(LLIM+igammar); i<=(ngrid_m+igammar); ++i)
    {
        R = i*dr;
        k1= i - igammar;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - B[j];
            rip = R + B[j];
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > gMAXR) kip = gMAXR;
                        
            temp1 = SimpsonIntegration(DD1[j],0,gMAXR,kin,kip);
            temp2 = SimpsonIntegration(DD2[j],0,gMAXR,kin,kip,B2[j],i);
            
            
            temp1 = temp1*PiB[j];
            temp2 = temp2*Pi*H[3][j];
            
            DES[j][k1] = temp1 + temp2;
        }
    }
        
    delete [] RI;
    delete [] RI2;
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] PiB;
    delete [] B2;
    delete [] DDH;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] DD2[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] DD2;
    delete [] NI;

}


void DerivElectroCorrel(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  RI;
    double*  RI2;
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  B2;
    double*  PiB;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double** NI;
    double*  DDB1;
    double*  DDB2;
    double** H;
    double** DD1;
    double** DD2;
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,gamma_s,Chi;
    short    hspecies;
    
    double   rin,rip;
    int      kin,kip,RMAX,BMAX,BMAX1,igammar;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    RI    = new double[hspecies]();
    RI2   = new double[hspecies]();
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B2    = new double[hspecies]();
    PiB   = new double[hspecies]();
    B     = new double[hspecies]();
    DDB1  = new double[hspecies]();
    DDB2  = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    DDH   = new double[4]();
    
    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    DD2    = new double*[hspecies]();
    
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    BMAX   = RMAX + igammar;
    BMAX1  = 2*BMAX;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[BMAX1+1]();
        DD2[i] = new double[BMAX1+1]();
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp4   = 0.5*D[i]/B[i];
        H[0][i] = 1.0;
        H[1][i] = temp4;
        H[2][i] = temp4*temp4;
        H[3][i] = H[2][i]*temp4;
        /////////////////////////////////////
        
        B2[i]  = B[i]*B[i];
        RI[i]  = 1.0/(Pi4*B[i]);
        RI2[i] = RI[i]/B[i];
        PiB[i] = Pi2*B[i];
    }
    
    WeightedDensity(gammab,rhoB,B,H,NI);
    
    gamma_k=CalculateGama(Z,D,NI);
    
    N3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        N0I[i] = NI[i][0];
        N1I[i] = NI[i][1];
        N2I[i] = NI[i][2];
        N3I[i] = NI[i][3];
        
        N3 += N3I[i];
    }
    BJGA  = BJ*gamma_k;
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gamma_k*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        
        temp1 += (Z[i]*N1I[i]*GammaS1[i]);
        temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
        temp3 += (N3I[i]*GammaS1[i]);
        temp4 += (Z[i]*Z[i]*N0I[i]*GammaS1[i]);
    }
    Chi     = 1.0/(1-N3+3*temp3);
    gamma_s = temp2*Chi;
    
    
    
    ES_EN = 0;
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

        for(short i=0; i<hspecies; ++i)
        {
            DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
            DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
            DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
            DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
            
            
            DDB1[i] = DDH[2]*H[2][i] + DDH[1]*RI[i]*H[1][i] + RI2[i]*DDH[0];
            DDB2[i] = DDH[3];
        }
        
        ES_EN = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp4+gamma_s*temp1)*BJGA;
        
    }//small loop
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            DD1[j][k] = DDB1[j];
            DD2[j][k] = DDB2[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
                    
        temp1 = SimpsonIntegration(DD1[j],0,BMAX1,kin,kip);
        temp2 = SimpsonIntegration(DD2[j],0,BMAX1,kin,kip,B2[j],BMAX);
        
        
        temp1 = temp1*PiB[j];
        temp2 = temp2*Pi*H[3][j];
        
        DESB[j] = BJGA*(temp1 + temp2);
    }
    delete [] RI;
    delete [] RI2;
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] PiB;
    delete [] B2;
    delete [] DDH;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] DDB1;
    delete [] DDB2;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] DD2[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] DD2;
    delete [] NI;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////     1      ///////////////////////////////////////////////////
void DerivElectroCorrelDirk(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                        double** rho,double** DES)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double** H;
    double** NI;
    double** DD1;
    double   temp1,temp2,temp3,N3,Pi2,Pi4,gammar,gamma_k,BJGA,Chi,gamma_s;
    short    hspecies;
    
    double   R,rin,rip;
    int      kin,kip,igammar,k1,gMAXR;
    
    hspecies= nspecies - 1;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B     = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    DDH   = new double[4]();
    
    
    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    gMAXR  = igammar + MAXR;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[gMAXR+1]();
        
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;

        H[0][i] = 1.0;
        H[1][i] = temp1;
        H[2][i] = temp2;
        H[3][i] = temp3;
        /////////////////////////////////////
    }
    
    //big loop
    for(int k=0; k<=gMAXR; ++k)
    {
        k1 = k - igammar;
        
        WeightedDensityDirk(k1,LLI,ULI,B,H,rho,NI);
        gamma_k=CalculateGama(Z,D,NI);
        
        N3 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            N0I[i] = NI[i][0];
            N1I[i] = NI[i][1];
            N2I[i] = NI[i][2];
            N3I[i] = NI[i][3];
            
            N3 += N3I[i];
        }
        
        BJGA  = BJ*gamma_k;
        
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gamma_k*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            
            temp1 += (Z[i]*N1I[i]*GammaS1[i]);
            temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
            temp3 += (N3I[i]*GammaS1[i]);
        }
        Chi     = 1.0/(1-N3+3*temp3);
        gamma_s = temp2*Chi;
        
        if(N3 < 1E-20)
        {
            for(short j1=0; j1<hspecies; ++j1)
            {
                DD1[j1][k] = 0;
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            for(short i=0; i<hspecies; ++i)
            {
                DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
                DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
                DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
                DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
                
                
                DD1[i][k] = BJGA*(DDH[3]*H[3][i] + DDH[2]*H[2][i] + DDH[1]*H[1][i] + DDH[0]);
            }
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=(LLIM+igammar); i<=(ngrid_m+igammar); ++i)
    {
        R = i*dr;
        k1= i - igammar;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - B[j];
            rip = R + B[j];
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > gMAXR) kip = gMAXR;
            
            temp1 = SimpsonIntegration(DD1[j],0,gMAXR,kin,kip);
            temp1 = temp1/(2*B[j]);
            
            DES[j][k1] = temp1;
        }
    }
    
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] DDH;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] NI;
    
}


void DerivElectroCorrelDirk(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double** NI;
    double*  DDB1;
    double** H;
    double** DD1;
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,gamma_s,Chi;
    short    hspecies;
    
    double   rin,rip;
    int      kin,kip,RMAX,BMAX,BMAX1,igammar;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B     = new double[hspecies]();
    DDB1  = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    DDH   = new double[4]();
    
    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    BMAX   = RMAX + igammar;
    BMAX1  = 2*BMAX;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[BMAX1+1]();
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        
        H[0][i] = 1.0;
        H[1][i] = temp1;
        H[2][i] = temp2;
        H[3][i] = temp3;
        /////////////////////////////////////
    }
    
    WeightedDensityDirk(gammab,rhoB,B,H,NI);
    
    gamma_k=CalculateGama(Z,D,NI);
    
    N3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        N0I[i] = NI[i][0];
        N1I[i] = NI[i][1];
        N2I[i] = NI[i][2];
        N3I[i] = NI[i][3];
        
        N3 += N3I[i];
    }
    BJGA  = BJ*gamma_k;
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gamma_k*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        
        temp1 += (Z[i]*N1I[i]*GammaS1[i]);
        temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
        temp3 += (N3I[i]*GammaS1[i]);
        temp4 += (Z[i]*Z[i]*N0I[i]*GammaS1[i]);
    }
    Chi     = 1.0/(1-N3+3*temp3);
    gamma_s = temp2*Chi;
    
    
    
    ES_EN = 0;
    if(N3 < 1E-20)
    {
        for(short j1=0; j1<hspecies; ++j1)
        {
            DDB1[j1] = 0;
        }
    }
    else //small loop
    {
        //Start: derivation of the hard sphere contribution
        if(N3 > 0.99) N3 = 0.99;
        
        for(short i=0; i<hspecies; ++i)
        {
            DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
            DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
            DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
            DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
            
            
            DDB1[i] = DDH[3]*H[3][i] + DDH[2]*H[2][i] + DDH[1]*H[1][i] + DDH[0];
        }
        
        ES_EN = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp4+gamma_s*temp1)*BJGA;
        
    }//small loop
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            DD1[j][k] = DDB1[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        temp1 = SimpsonIntegration(DD1[j],0,BMAX1,kin,kip);
        temp1 = temp1/(2*B[j]);
        
        DESB[j] = BJGA*temp1;
    }
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] DDH;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] DDB1;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] NI;
    
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////     V      ///////////////////////////////////////////////////
void DerivElectroCorrelV(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                        double** rho,double** DES)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  B2;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double** H;
    double** NI;
    double** DD1;
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,Chi,gamma_s;
    short    hspecies;
    
    double   R,rin,rip;
    int      kin,kip,igammar,k1,gMAXR;
    
    hspecies= nspecies - 1;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B2    = new double[hspecies]();
    B     = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    DDH   = new double[4]();
    
    
    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    gMAXR  = igammar + MAXR;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[gMAXR+1]();
        
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        temp4   = Pi4*B[i]*B[i]*B[i]/3.0;
        
        H[0][i] = 1.0/temp4;
        H[1][i] = temp1/temp4;
        H[2][i] = temp2/temp4;
        H[3][i] = temp3/temp4;
        /////////////////////////////////////
        
        B2[i]  = B[i]*B[i];
    }
    
    //big loop
    for(int k=0; k<=MAXR; ++k)
    {
        k1 = k - igammar;
        
        WeightedDensityV(k1,LLI,ULI,B,H,rho,NI);
        gamma_k=CalculateGama(Z,D,NI);
        
        N3 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            N0I[i] = NI[i][0];
            N1I[i] = NI[i][1];
            N2I[i] = NI[i][2];
            N3I[i] = NI[i][3];
            
            N3 += N3I[i];
        }
        
        BJGA  = BJ*gamma_k;
        
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gamma_k*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            
            temp1 += (Z[i]*N1I[i]*GammaS1[i]);
            temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
            temp3 += (N3I[i]*GammaS1[i]);
        }
        Chi     = 1.0/(1-N3+3*temp3);
        gamma_s = temp2*Chi;
        
        if(N3 < 1E-20)
        {
            for(short j1=0; j1<hspecies; ++j1)
            {
                DD1[j1][k] = 0;
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            for(short i=0; i<hspecies; ++i)
            {
                DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
                DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
                DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
                DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
                
                
                DD1[i][k] = BJGA*(DDH[3]*H[3][i] + DDH[2]*H[2][i] + DDH[1]*H[1][i] + DDH[0]*H[0][i]);
            }
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=(LLIM+igammar); i<=(ngrid_m+igammar); ++i)
    {
        R = i*dr;
        k1= i - igammar;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - B[j];
            rip = R + B[j];
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > gMAXR) kip = gMAXR;
            
            temp1 = SimpsonIntegration(DD1[j],0,gMAXR,kin,kip,B2[j],i);
            
            DES[j][k1] = temp1*Pi;
        }
    }
    
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] B2;
    delete [] DDH;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] NI;
    
}


void DerivElectroCorrelV(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  DDH;
    double*  B2;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double** NI;
    double*  DDB1;
    double** H;
    double** DD1;
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,gamma_s,Chi;
    short    hspecies;
    
    double   rin,rip;
    int      kin,kip,RMAX,BMAX,BMAX1,igammar;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B2    = new double[hspecies]();
    B     = new double[hspecies]();
    DDB1  = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    DDH   = new double[4]();
    
    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    BMAX   = RMAX + igammar;
    BMAX1  = 2*BMAX;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[BMAX1+1]();
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        temp4   = Pi4*B[i]*B[i]*B[i]/3.0;
        
        H[0][i] = 1.0/temp4;
        H[1][i] = temp1/temp4;
        H[2][i] = temp2/temp4;
        H[3][i] = temp3/temp4;
        /////////////////////////////////////
        
        B2[i]  = B[i]*B[i];
    }
    
    WeightedDensityV(gammab,rhoB,B,H,NI);
    
    gamma_k=CalculateGama(Z,D,NI);
    
    N3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        N0I[i] = NI[i][0];
        N1I[i] = NI[i][1];
        N2I[i] = NI[i][2];
        N3I[i] = NI[i][3];
        
        N3 += N3I[i];
    }
    BJGA  = BJ*gamma_k;
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gamma_k*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        
        temp1 += (Z[i]*N1I[i]*GammaS1[i]);
        temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
        temp3 += (N3I[i]*GammaS1[i]);
        temp4 += (Z[i]*Z[i]*N0I[i]*GammaS1[i]);
    }
    Chi     = 1.0/(1-N3+3*temp3);
    gamma_s = temp2*Chi;
    
    
    
    ES_EN = 0;
    if(N3 < 1E-20)
    {
        for(short j1=0; j1<hspecies; ++j1)
        {
            DDB1[j1] = 0;
        }
    }
    else //small loop
    {
        //Start: derivation of the hard sphere contribution
        if(N3 > 0.99) N3 = 0.99;
        
        for(short i=0; i<hspecies; ++i)
        {
            DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
            DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
            DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
            DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
            
            
            DDB1[i] = DDH[3]*H[3][i] + DDH[2]*H[2][i] + DDH[1]*H[1][i] + DDH[0]*H[0][i];
        }
        
        ES_EN = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp4+gamma_s*temp1)*BJGA;
        
    }//small loop
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            DD1[j][k] = DDB1[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        temp1 = SimpsonIntegration(DD1[j],0,BMAX1,kin,kip,B2[j],BMAX);
        
        DESB[j] = BJGA*temp1*Pi;
    }
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] B2;
    delete [] DDH;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] DDB1;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] NI;
    
}




/*
///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////     1      ///////////////////////////////////////////////////
void DerivElectroCorrelDirk(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                            double** rho,double** DES)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double*  NI;
    double** H;
    double** DD1;
    double   temp1,temp2,temp3,N3,Pi2,Pi4,gammar,gamma_k,BJGA,Chi,gamma_s,DDH;
    short    hspecies;
    
    double   R,rin,rip;
    int      kin,kip,igammar;
    
    hspecies= nspecies - 1;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NI    = new double[hspecies]();
    B     = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    
    H      = new double*[3]();
    
    DD1    = new double*[hspecies]();
    for(short i=0; i<3; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[MAXR+1]();
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        
        H[0][i] = temp1;
        H[1][i] = temp2;
        H[2][i] = temp3;
        /////////////////////////////////////
    }
    
    //big loop
    for(int k=0; k<=MAXR; ++k)
    {
        WeightedDensityDirk(k,LLI,ULI,B,rho,NI);
        gamma_k=CalculateGama(Z,D,H,NI);
        
        N3 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            N0I[i] = NI[i];
            N1I[i] = NI[i]*H[0][i];
            N2I[i] = NI[i]*H[1][i];
            N3I[i] = NI[i]*H[2][i];
            
            N3 += N3I[i];
        }
        
        BJGA  = BJ*gamma_k;
        
        temp2 = 0;
        temp3 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gamma_k*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            
            temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
            temp3 += (N3I[i]*GammaS1[i]);
        }
        Chi     = 1.0/(1-N3+3*temp3);
        gamma_s = temp2*Chi;
        
        if(N3 < 1E-20)
        {
            for(short j1=0; j1<hspecies; ++j1)
            {
                DD1[j1][k] = 0;
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            for(short i=0; i<hspecies; ++i)
            {
                temp1 = -(gamma_s*Z[i]*D[i]*GammaS1[i]);
                temp2 = -(Z[i]*Z[i]*GammaS1[i]);
                temp3 = gamma_s*gamma_s*gamma_k*(2-GammaS[i])*GammaS1[i]*H[2][i]/Pi2;
                
                DDH   = temp3 + temp2 + temp1;
                
                DD1[i][k] = BJGA*DDH;
            }
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=LLIM; i<=ngrid_m; ++i)
    {
        R = i*dr;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - B[j];
            rip = R + B[j];
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > MAXR) kip = MAXR;
            
            temp1 = SimpsonIntegration(DD1[j],0,MAXR,kin,kip);
            temp1 = temp1/(2*B[j]);
            
            DES[j][i] = temp1;
        }
    }
    
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    
    for(short i=0; i<3; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
    }
    delete [] DD1;
    delete [] NI;
    
}


void DerivElectroCorrelDirk(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double*  NI;
    double*  DDB1;
    double** H;
    double** DD1;
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,gamma_s,Chi,DDH;
    short    hspecies;
    
    double   rin,rip;
    int      kin,kip,RMAX,BMAX,BMAX1,igammar;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NI    = new double[hspecies]();
    B     = new double[hspecies]();
    DDB1  = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    
    H      = new double*[3]();
    DD1    = new double*[hspecies]();
    
    for(short i=0; i<3; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    BMAX   = RMAX + igammar;
    BMAX1  = 2*BMAX;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[BMAX1+1]();
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        
        H[0][i] = temp1;
        H[1][i] = temp2;
        H[2][i] = temp3;
        /////////////////////////////////////
    }
    
    WeightedDensityDirk(gammab,rhoB,B,NI);
    
    gamma_k=CalculateGama(Z,D,H,NI);
    
    N3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        N0I[i] = NI[i];
        N1I[i] = NI[i]*H[0][i];
        N2I[i] = NI[i]*H[1][i];
        N3I[i] = NI[i]*H[2][i];
        
        N3 += N3I[i];
    }
    BJGA  = BJ*gamma_k;
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gamma_k*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        
        temp1 += (Z[i]*N1I[i]*GammaS1[i]);
        temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
        temp3 += (N3I[i]*GammaS1[i]);
        temp4 += (Z[i]*Z[i]*N0I[i]*GammaS1[i]);
    }
    Chi     = 1.0/(1-N3+3*temp3);
    gamma_s = temp2*Chi;
    
    
    
    ES_EN = 0;
    if(N3 < 1E-20)
    {
        for(short j1=0; j1<hspecies; ++j1)
        {
            DDB1[j1] = 0;
        }
    }
    else //small loop
    {
        //Start: derivation of the hard sphere contribution
        if(N3 > 0.99) N3 = 0.99;
        
        for(short i=0; i<hspecies; ++i)
        {
            temp2 = -(gamma_s*Z[i]*D[i]*GammaS1[i]);
            temp3 = gamma_s*gamma_s*gamma_k*(2-GammaS[i])*GammaS1[i]*H[2][i]/Pi2;
            
            DDH   = temp3 + temp2 - (Z[i]*Z[i]*GammaS1[i]);
            
            DDB1[i] = DDH;
        }
        
        ES_EN = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp4+gamma_s*temp1)*BJGA;
        
    }//small loop
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            DD1[j][k] = DDB1[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        temp1 = SimpsonIntegration(DD1[j],0,BMAX1,kin,kip);
        temp1 = temp1/(2*B[j]);
        
        DESB[j] = BJGA*temp1;
    }
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] DDB1;
    
    for(short i=0; i<3; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
    }
    delete [] DD1;
    delete [] NI;
    
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////     V      ///////////////////////////////////////////////////
void DerivElectroCorrelV(int MAXR,double gammab,int* LLI,int* ULI,float* D,float* Z,
                         double** rho,double** DES)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  B2;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double*  NI;
    double** H;
    double** DD1;
    double   temp1,temp2,temp3,N3,Pi2,Pi4,gammar,gamma_k,BJGA,Chi,gamma_s,DDH;
    short    hspecies;
    
    double   R,rin,rip;
    int      kin,kip,igammar;
    
    hspecies= nspecies - 1;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NI    = new double[hspecies]();
    B2    = new double[hspecies]();
    B     = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    
    
    H      = new double*[3]();
    DD1    = new double*[hspecies]();
    for(short i=0; i<3; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[MAXR+1]();
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        
        H[0][i] = temp1;
        H[1][i] = temp2;
        H[2][i] = temp3;
        /////////////////////////////////////
        
        B2[i]  = B[i]*B[i];
    }
    
    //big loop
    for(int k=0; k<=MAXR; ++k)
    {
        WeightedDensityV(k,LLI,ULI,B,rho,NI);
        gamma_k=CalculateGama(Z,D,H,NI);
        
        N3 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            N0I[i] = NI[i];
            N1I[i] = NI[i]*H[0][i];
            N2I[i] = NI[i]*H[1][i];
            N3I[i] = NI[i]*H[2][i];
            
            N3 += N3I[i];
        }
        
        BJGA  = BJ*gamma_k;
        
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gamma_k*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            
            //temp1 += (Z[i]*N1I[i]*GammaS1[i]);
            temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
            temp3 += (N3I[i]*GammaS1[i]);
        }
        Chi     = 1.0/(1-N3+3*temp3);
        gamma_s = temp2*Chi;
        
        if(N3 < 1E-20)
        {
            for(short j1=0; j1<hspecies; ++j1)
            {
                DD1[j1][k] = 0;
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            for(short i=0; i<hspecies; ++i)
            {
                temp1 = -(gamma_s*Z[i]*D[i]*GammaS1[i]);
                temp2 = -(Z[i]*Z[i]*GammaS1[i]);
                temp3 = gamma_s*gamma_s*gamma_k*(2-GammaS[i])*GammaS1[i]*H[2][i]/Pi2;
                
                DDH   = temp3 + temp2 + temp1;
                
                
                DD1[i][k] = BJGA*DDH;
            }
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=LLIM; i<=ngrid_m; ++i)
    {
        R = i*dr;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - B[j];
            rip = R + B[j];
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > MAXR) kip = MAXR;
            
            temp1 = SimpsonIntegration(DD1[j],0,MAXR,kin,kip,B2[j],i);
            
            //DES[j][i] = temp1*Pi;
            DES[j][i] = temp1*3/(4*B2[j]*B[j]);
        }
    }
    
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] B2;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    
    for(short i=0; i<3; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
    }
    delete [] DD1;
    delete [] NI;
    
}


void DerivElectroCorrelV(double* rhoB,double gammab,float* D,float* Z,double* DESB,double& ES_EN)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  B2;
    double*  B;
    double*  GammaS;
    double*  GammaS1;
    double*  NI;
    double*  DDB1;
    double** H;
    double** DD1;
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,gamma_s,Chi,DDH;
    short    hspecies;
    
    double   rin,rip;
    int      kin,kip,RMAX,BMAX,BMAX1,igammar;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NI    = new double[hspecies]();
    B2    = new double[hspecies]();
    B     = new double[hspecies]();
    DDB1  = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    
    H      = new double*[3]();
    DD1    = new double*[hspecies]();
    
    for(short i=0; i<3; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    BMAX   = RMAX + igammar;
    BMAX1  = 2*BMAX;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[BMAX1+1]();
        B[i] = 0.5*D[i] + igammar*dr;
        /////////////////////////////////////
        temp1   = 0.5*D[i];
        temp2   = Pi4*temp1*temp1;
        temp3   = temp2*temp1/3.0;
        
        H[0][i] = temp1;
        H[1][i] = temp2;
        H[2][i] = temp3;
        /////////////////////////////////////
        
        B2[i]  = B[i]*B[i];
    }
    
    WeightedDensityV(gammab,rhoB,B,NI);
    
    gamma_k=CalculateGama(Z,D,H,NI);
    
    N3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        N0I[i] = NI[i];
        N1I[i] = NI[i]*H[0][i];
        N2I[i] = NI[i]*H[1][i];
        N3I[i] = NI[i]*H[2][i];
        
        N3 += N3I[i];
    }
    BJGA  = BJ*gamma_k;
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gamma_k*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        
        temp1 += (Z[i]*N1I[i]*GammaS1[i]);
        temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
        temp3 += (N3I[i]*GammaS1[i]);
        temp4 += (Z[i]*Z[i]*N0I[i]*GammaS1[i]);
    }
    Chi     = 1.0/(1-N3+3*temp3);
    gamma_s = temp2*Chi;
    
    
    
    ES_EN = 0;
    if(N3 < 1E-20)
    {
        for(short j1=0; j1<hspecies; ++j1)
        {
            DDB1[j1] = 0;
        }
    }
    else //small loop
    {
        //Start: derivation of the hard sphere contribution
        if(N3 > 0.99) N3 = 0.99;
        
        for(short i=0; i<hspecies; ++i)
        {
            temp2 = -(gamma_s*Z[i]*D[i]*GammaS1[i]);
            temp3 = gamma_s*gamma_s*gamma_k*(2-GammaS[i])*GammaS1[i]*H[2][i]/Pi2;
            
            DDH   = temp3 + temp2 - (Z[i]*Z[i]*GammaS1[i]);
            
            
            DDB1[i] = DDH;
        }
        
        ES_EN = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp4+gamma_s*temp1)*BJGA;
        
    }//small loop
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            DD1[j][k] = DDB1[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        temp1 = SimpsonIntegration(DD1[j],0,BMAX1,kin,kip,B2[j],BMAX);
        
        //DESB[j] = BJGA*temp1*Pi;
        DESB[j] = BJGA*temp1*3/(4*B2[j]*B[j]);
    }
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] B2;
    delete [] B;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] DDB1;
    
    for(short i=0; i<3; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
    }
    delete [] DD1;
    delete [] NI;
    
}
 
 */

