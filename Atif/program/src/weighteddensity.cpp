//***********solve the possion equation****************//
#include "clibrary.h"
#include "weighteddensity.h"
#include "simpsonintegration.h"
#include "constantnum.h"

extern double dr;
extern int ngrid;
extern int DMAX;
extern short  nspecies;


void WeightedDensity(int i,int* LLI,int* ULI,float* D,double** rho,double** NI)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double* N0I;
    double* N1I;
    double* N2I;
    double* N3I;
    double* NV1I;
    double* NV2I;
    double  R,rin,rip,D2,Pi2,Pi4,coe;
    int kin,kip;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    NV1I  = new double[hspecies]();
    NV2I  = new double[hspecies]();
    
    
    R   = i*dr;
    
    Pi2 = 2*Pi;
    Pi4 = 2*Pi2;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - 0.5*D[j];
        rip = R + 0.5*D[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        D2  = D[j]*D[j]*0.25;
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        if(kin < kip)
        {
            
            N2I[j] = SimpsonIntegration(rho[j],0,ngrid,kin,kip);
            N3I[j] = SimpsonIntegration(rho[j],0,ngrid,kin,kip,D2,i);
            NV2I[j]= SimpsonIntegration(rho[j],0,ngrid,kin,kip,i);
            NV2I[j]= -NV2I[j];
            
            coe     = 1.0/(Pi2*D[j]);
            
            N2I[j]  = N2I[j]*Pi*D[j];
            N3I[j]  = N3I[j]*Pi;
            NV2I[j] = NV2I[j]*Pi2;
            N0I[j]  = N2I[j]/(Pi4*D2);
            N1I[j]  = coe*N2I[j];
            NV1I[j] = coe*NV2I[j];
            
            
            NI[j][0] = N0I[j];
            NI[j][1] = N1I[j];
            NI[j][2] = N2I[j];
            NI[j][3] = N3I[j];
            NI[j][4] = NV1I[j];
            NI[j][5] = NV2I[j];
            
        }
        else
        {
            NI[j][0] = 0;
            NI[j][1] = 0;
            NI[j][2] = 0;
            NI[j][3] = 0;
            NI[j][4] = 0;
            NI[j][5] = 0;
        }
    }
    

    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] NV1I;
    delete [] NV2I;
}


void WeightedDensity(float* D,double* rhoB,double** NIB)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double* N0I;
    double* N1I;
    double* N2I;
    double* N3I;
    double** rho_bk;
    double  rin,rip,D2,Pi2,Pi4;
    int   kin,kip,RMAX;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    
    rho_bk= new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        rho_bk[j] = new double[DMAX+1]();
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=DMAX; ++k)
        {
            rho_bk[j][k] = rhoB[j];
        }
    }
    
    
    Pi2 = 2*Pi;
    Pi4 = 2*Pi2;
    for(short j=0; j<hspecies; ++j)
    {
        rin = RMAX*dr - 0.5*D[j];
        rip = RMAX*dr + 0.5*D[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        D2  = D[j]*D[j]*0.25;
        
        if(kin < kip)
        {
            
            N2I[j] = SimpsonIntegration(rho_bk[j],0,DMAX,kin,kip);
            N3I[j] = SimpsonIntegration(rho_bk[j],0,DMAX,kin,kip,D2,RMAX);
            
            N2I[j]  = N2I[j]*Pi*D[j];
            N3I[j]  = N3I[j]*Pi;
            N0I[j]  = N2I[j]/(Pi4*D2);
            N1I[j]  = N2I[j]/(Pi2*D[j]);
            
            
            NIB[j][0] = N0I[j];
            NIB[j][1] = N1I[j];
            NIB[j][2] = N2I[j];
            NIB[j][3] = N3I[j];

            
        }
        else
        {
            NIB[j][0] = 0;
            NIB[j][1] = 0;
            NIB[j][2] = 0;
            NIB[j][3] = 0;
        }
    }
    

    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    for(short j=0; j<hspecies; ++j)
    {
        delete [] rho_bk[j];
    }
    delete [] rho_bk;
}


////////////////////weighted for electrostatic correlation: 1,R,S,V///////////////////////////////
void WeightedDensity(int i,int* LLI,int* ULI,double* B,double** H,double** rho,double** QI)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double* Q0I;
    double* Q1I;
    double* Q2I;
    double* Q3I;
    double  R,rin,rip,B2,Pi2,Pi4;
    int kin,kip;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    Q0I   = new double[hspecies]();
    Q1I   = new double[hspecies]();
    Q2I   = new double[hspecies]();
    Q3I   = new double[hspecies]();
    
    
    R   = i*dr;
    Pi2 = 2*Pi;
    Pi4 = 2*Pi2;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - B[j];
        rip = R + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        B2  = B[j]*B[j];
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        if(kin < kip)
        {
            
            Q2I[j] = SimpsonIntegration(rho[j],0,ngrid,kin,kip);
            Q3I[j] = SimpsonIntegration(rho[j],0,ngrid,kin,kip,B2,i);
            
            Q2I[j]  = Q2I[j]*Pi2*B[j];
            Q3I[j]  = Q3I[j]*Pi;
            Q0I[j]  = Q2I[j]/(Pi4*B2);
            Q1I[j]  = Q2I[j]/(Pi4*B[j]);
            
            
            QI[j][0] = Q0I[j]*H[0][j];
            QI[j][1] = Q1I[j]*H[1][j];
            QI[j][2] = Q2I[j]*H[2][j];
            QI[j][3] = Q3I[j]*H[3][j];
            
        }
        else
        {
            QI[j][0] = 0;
            QI[j][1] = 0;
            QI[j][2] = 0;
            QI[j][3] = 0;
        }
    }
    
    
    delete [] Q0I;
    delete [] Q1I;
    delete [] Q2I;
    delete [] Q3I;
}



void WeightedDensity(double gammab,double* rhoB,double* B,double** H,double** QIB)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double* Q0I;
    double* Q1I;
    double* Q2I;
    double* Q3I;
    double** rho_bk;
    double  rin,rip,B2,Pi2,Pi4;
    int kin,kip,RMAX,BMAX,igammar,BMAX1;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    igammar=round(0.5/(dr*gammab));
    BMAX   = RMAX + igammar;
    BMAX1  = BMAX*2;
    
    Q0I   = new double[hspecies]();
    Q1I   = new double[hspecies]();
    Q2I   = new double[hspecies]();
    Q3I   = new double[hspecies]();
    
    rho_bk= new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        rho_bk[j] = new double[BMAX1+1]();
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            rho_bk[j][k] = rhoB[j];
        }
    }
    
    
    Pi2 = 2*Pi;
    Pi4 = 2*Pi2;
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        B2  = B[j]*B[j];
        
        if(kin < kip)
        {
            
            Q2I[j] = SimpsonIntegration(rho_bk[j],0,BMAX1,kin,kip);
            Q3I[j] = SimpsonIntegration(rho_bk[j],0,BMAX1,kin,kip,B2,BMAX);
            
            Q2I[j]  = Q2I[j]*Pi2*B[j];
            Q3I[j]  = Q3I[j]*Pi;
            Q0I[j]  = Q2I[j]/(Pi4*B2);
            Q1I[j]  = Q2I[j]/(Pi4*B[j]);
            
            
            QIB[j][0] = Q0I[j]*H[0][j];
            QIB[j][1] = Q1I[j]*H[1][j];
            QIB[j][2] = Q2I[j]*H[2][j];
            QIB[j][3] = Q3I[j]*H[3][j];
            
        }
        else
        {
            QIB[j][0] = 0;
            QIB[j][1] = 0;
            QIB[j][2] = 0;
            QIB[j][3] = 0;
        }
    }
    
    
    delete [] Q0I;
    delete [] Q1I;
    delete [] Q2I;
    delete [] Q3I;
    
    for(short j=0; j<hspecies; ++j)
    {
        delete [] rho_bk[j];
    }
    delete [] rho_bk;

}



////////////////////weighted for electrostatic correlation: 1 ///////////////////////////////
void WeightedDensityDirk(int i,int* LLI,int* ULI,double* B,double** H,double** rho,double** QI)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double  R,rin,rip,temp;
    int kin,kip;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    
    R   = i*dr;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - B[j];
        rip = R + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        if(kin < kip)
        {
            temp   = SimpsonIntegration(rho[j],0,ngrid,kin,kip);
            temp   = temp/(2*B[j]);
            
            
            QI[j][0] = temp*H[0][j];
            QI[j][1] = temp*H[1][j];
            QI[j][2] = temp*H[2][j];
            QI[j][3] = temp*H[3][j];
            
        }
        else
        {
            QI[j][0] = 0;
            QI[j][1] = 0;
            QI[j][2] = 0;
            QI[j][3] = 0;
        }
    }
    
}



void WeightedDensityDirk(double gammab,double* rhoB,double* B,double** H,double** QIB)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double** rho_bk;
    double  rin,rip,temp;
    int kin,kip,RMAX,BMAX,igammar,BMAX1;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    igammar=round(0.5/(dr*gammab));
    BMAX   = RMAX + igammar;
    BMAX1  = BMAX*2;
        
    rho_bk= new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        rho_bk[j] = new double[BMAX1+1]();
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            rho_bk[j][k] = rhoB[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        if(kin < kip)
        {
            
            temp = SimpsonIntegration(rho_bk[j],0,BMAX1,kin,kip);
            temp = temp/(2*B[j]);
            
            QIB[j][0] = temp*H[0][j];
            QIB[j][1] = temp*H[1][j];
            QIB[j][2] = temp*H[2][j];
            QIB[j][3] = temp*H[3][j];
            
        }
        else
        {
            QIB[j][0] = 0;
            QIB[j][1] = 0;
            QIB[j][2] = 0;
            QIB[j][3] = 0;
        }
    }
    
        
    for(short j=0; j<hspecies; ++j)
    {
        delete [] rho_bk[j];
    }
    delete [] rho_bk;
    
}


////////////////////weighted for electrostatic correlation: V ///////////////////////////////
void WeightedDensityV(int i,int* LLI,int* ULI,double* B,double** H,double** rho,double** QI)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double  R,rin,rip,B2,temp;
    int kin,kip;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    R   = i*dr;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - B[j];
        rip = R + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        B2  = B[j]*B[j];
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        if(kin < kip)
        {
            temp = SimpsonIntegration(rho[j],0,ngrid,kin,kip,B2,i);
            temp = temp*Pi;
            
            QI[j][0] = temp*H[0][j];
            QI[j][1] = temp*H[1][j];
            QI[j][2] = temp*H[2][j];
            QI[j][3] = temp*H[3][j];
            
        }
        else
        {
            QI[j][0] = 0;
            QI[j][1] = 0;
            QI[j][2] = 0;
            QI[j][3] = 0;
        }
    }
    
}



void WeightedDensityV(double gammab,double* rhoB,double* B,double** H,double** QIB)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr
    double** rho_bk;
    double  rin,rip,B2,temp;
    int kin,kip,RMAX,BMAX,igammar,BMAX1;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    igammar=round(0.5/(dr*gammab));
    BMAX   = RMAX + igammar;
    BMAX1  = BMAX*2;
        
    rho_bk= new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        rho_bk[j] = new double[BMAX1+1]();
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            rho_bk[j][k] = rhoB[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        B2  = B[j]*B[j];
        
        if(kin < kip)
        {
            temp = SimpsonIntegration(rho_bk[j],0,BMAX1,kin,kip,B2,BMAX);
            temp = temp*Pi;

            QIB[j][0] = temp*H[0][j];
            QIB[j][1] = temp*H[1][j];
            QIB[j][2] = temp*H[2][j];
            QIB[j][3] = temp*H[3][j];
            
        }
        else
        {
            QIB[j][0] = 0;
            QIB[j][1] = 0;
            QIB[j][2] = 0;
            QIB[j][3] = 0;
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        delete [] rho_bk[j];
    }
    delete [] rho_bk;
    
}




void WeightedDensityDirk(int i,int* LLI,int* ULI,double* B,double** rho,double* QI)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double  R,rin,rip,temp;
    int kin,kip;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    R   = i*dr;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - B[j];
        rip = R + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        if(kin < kip)
        {
            temp   = SimpsonIntegration(rho[j],0,ngrid,kin,kip);
            temp   = temp/(2*B[j]);
            
            
            QI[j] = temp;
            
        }
        else
        {
            QI[j] = 0;
        }
    }
}



void WeightedDensityDirk(double gammab,double* rhoB,double* B,double* QIB)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double** rho_bk;
    double  rin,rip,temp;
    int kin,kip,RMAX,BMAX,igammar,BMAX1;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    igammar=round(0.5/(dr*gammab));
    BMAX   = RMAX + igammar;
    BMAX1  = BMAX*2;
        
    rho_bk= new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        rho_bk[j] = new double[BMAX1+1]();
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            rho_bk[j][k] = rhoB[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        if(kin < kip)
        {
            
            temp = SimpsonIntegration(rho_bk[j],0,BMAX1,kin,kip);
            temp = temp/(2*B[j]);
            
            QIB[j] = temp;
            
        }
        else
        {
            QIB[j] = 0;
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        delete [] rho_bk[j];
    }
    delete [] rho_bk;
    
}


////////////////////weighted for electrostatic correlation: V ///////////////////////////////
void WeightedDensityV(int i,int* LLI,int* ULI,double* B,double** rho,double* QI)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double  R,rin,rip,B2,temp,temp1;
    int kin,kip;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    R   = i*dr;
    for(short j=0; j<hspecies; ++j)
    {
        rin = R - B[j];
        rip = R + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        B2   = B[j]*B[j];
        temp1= 4.0*B2*B[j]/3.0;
        if(kin < LLI[j]) kin = LLI[j];
        if(kip > ULI[j]) kip = ULI[j];
        
        if(kin < kip)
        {
            temp = SimpsonIntegration(rho[j],0,ngrid,kin,kip,B2,i);
            //temp = temp*Pi;
            
            QI[j] = temp/temp1;
        }
        else
        {
            QI[j] = 0;
        }
    }
    
}



void WeightedDensityV(double gammab,double* rhoB,double* B,double* QIB)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double** rho_bk;
    double  rin,rip,B2,temp,temp1;
    int kin,kip,RMAX,BMAX,igammar,BMAX1;
    short hspecies;
    
    hspecies= nspecies - 1;
    
    RMAX  = DMAX/2;
    
    igammar=round(0.5/(dr*gammab));
    BMAX   = RMAX + igammar;
    BMAX1  = BMAX*2;
    
    
    rho_bk= new double*[hspecies]();
    for(short j=0; j<hspecies; ++j)
    {
        rho_bk[j] = new double[BMAX1+1]();
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        for(int k=0; k<=BMAX1; ++k)
        {
            rho_bk[j][k] = rhoB[j];
        }
    }
    
    for(short j=0; j<hspecies; ++j)
    {
        rin = BMAX*dr - B[j];
        rip = BMAX*dr + B[j];
        kin = round(rin/dr);
        kip = round(rip/dr);
        
        B2   = B[j]*B[j];
        temp1= 4.0*B2*B[j]/3.0;
        
        if(kin < kip)
        {
            temp = SimpsonIntegration(rho_bk[j],0,BMAX1,kin,kip,B2,BMAX);
            //temp = temp*Pi;

            QIB[j] = temp/temp1;
        }
        else
        {
            QIB[j] = 0;
        }
    }
        
    for(short j=0; j<hspecies; ++j)
    {
        delete [] rho_bk[j];
    }
    delete [] rho_bk;
    
}

