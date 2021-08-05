//***********solve the possion equation****************//
#include "clibrary.h"
#include "poissonequation.h"
//#include "rombergintegration.h"
#include "simpsonintegration.h"
#include "constantnum.h"
using namespace std;

extern double dr;
extern double kapa;
extern double BJ;
extern double boundary;
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int ngrid_b; //the number of grids: the middle between two surfaces
extern short nspecies;

//In this code, phiR means the middle potential of the system
void PoissonEquationSingle(float* Z,double** rho,double phiL,double phiR,double* Psi)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //int     ngridm1; //,i0
    double  R1,lGau,rGau,coe,rhoD,phi0,sigma0,sigmai,size_b,R_L,dr_size_b;
    double* rhoZ;
    double* phii;
    short   hspecies;
    
    hspecies = nspecies - 1;
    
    //coe = 2*Pi*BJ*dr;
    coe = 2*Pi*BJ*dr;
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    size_b   = ngrid_b*dr;
    //ngridm1 = ngrid_m + 1;
    dr_size_b= dr/size_b;
    //if(ngrid%2 != 0) size_mid = ngrid_m*dr + 0.5*dr;
    //ngrid0 = ngrid;
    //if(ngrid%2 != 0) ngrid0 = ngrid - 1;
    
    rhoZ = new double[ngrid_m+1]();
    phii = new double[ngrid_m+1]();
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        rhoZ[k] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            if(Z[j] != 0) rhoZ[k] += Z[j]*rho[j][k];
        }
    }
    
    //sigma0 = RombergIntegration(rhoZ,0,ngridm1,0,ngrid_m);
    sigma0 = SimpsonIntegration(rhoZ,0,ngrid_m,0,ngrid_b);
    
    
    phi0    = 0;
    phii[0] = 0;
    for(int k=0; k<ngrid_b; ++k)
    {
        R1    = (k + lGau)*dr;
        rhoD  = rhoZ[k] + lGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        R1    = (k + rGau)*dr;
        rhoD  = rhoZ[k] + rGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        phii[k+1] = phi0;
    }
    
    Psi[0]    = phiL;
    //Psi[ngrid]= phiL;
    for(int i=1; i<=ngrid_b; ++i)
    {
        R_L    = i*dr_size_b;
        //sigmai = RombergIntegration(rhoZ,0,ngridm1,0,i);
        sigmai = SimpsonIntegration(rhoZ,0,ngrid_m,0,i);
        
        Psi[i] = phiL + (phiR-phiL)*R_L + 2*coe*i*(sigma0-sigmai) + coe*(phii[i]-phii[ngrid_b]*R_L);
        Psi[ngrid-i]= Psi[i];
    }
    
    for(int i=(ngrid_b+1); i<=ngrid_m; ++i)
    {
        Psi[i] = 0;
        
        Psi[ngrid-i]= 0;
    }
    
    delete [] rhoZ;
    delete [] phii;
}


//In this code, phiR means the middle potential of the system
void PoissonEquationSingle(float* Z,double** rho,double* Psi)
{
    double  R1,lGau,rGau,coe,rhoD,phi0,sigma0,sigmai;
    double* rhoZ;
    double* phii;
    
    short   hspecies;
    
    hspecies = nspecies - 1;
    
    //coe = 2*Pi*BJ*dr;
    coe = 2*Pi*BJ*dr;
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    
    rhoZ = new double[ngrid_m+1]();
    phii = new double[ngrid_m+1]();
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        rhoZ[k] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            if(Z[j] != 0) rhoZ[k] += Z[j]*rho[j][k];
        }
    }
    
    sigma0 = SimpsonIntegration(rhoZ,0,ngrid_m,0,ngrid_b);
    
    
    phi0    = 0;
    phii[0] = 0;
    for(int k=0; k<ngrid_b; ++k)
    {
        R1    = (k + lGau)*dr;
        rhoD  = rhoZ[k] + lGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        R1    = (k + rGau)*dr;
        rhoD  = rhoZ[k] + rGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        phii[k+1] = phi0;
    }
    
    for(int i=0; i<=ngrid_b; ++i)
    {
        sigmai = SimpsonIntegration(rhoZ,0,ngrid_m,0,i);
        
        Psi[i] = 2*coe*i*(sigma0-sigmai) + coe*(phii[i]-phii[ngrid_b]);
        
        Psi[ngrid-i]= Psi[i];
    }
    
    for(int i=(ngrid_b+1); i<=ngrid_m; ++i)
    {
        Psi[i] = 0;
        
        Psi[ngrid-i]= 0;
    }
    
    delete [] rhoZ;
    delete [] phii;
}

void PoissonEquationTwo(float* Z,double** rho,double* Psi)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    double  R1,lGau,rGau,coe,rhoD,phi0,sigma0,sigmai;
    double* rhoZ;
    double* phii;
    int     ngridm1;
    short   hspecies;
    
    hspecies = nspecies - 1;
    ngridm1  = ngrid_m + 1;

    
    coe = 2*Pi*BJ*dr;
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    
    rhoZ = new double[ngridm1+1]();
    phii = new double[ngrid_m+1]();
    
    for(int k=0; k<=ngridm1; ++k)
    {
        rhoZ[k] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            if(Z[j] != 0) rhoZ[k] += Z[j]*rho[j][k];
        }
    }
    
    sigma0 = SimpsonIntegration(rhoZ,0,ngridm1,0,ngrid_m);

    sigma0 = 2*sigma0;
    
    phi0    = 0;
    phii[0] = 0;
    for(int k=0; k<ngrid_m; ++k)
    {
        R1    = (k + lGau)*dr;
        rhoD  = rhoZ[k] + lGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        R1    = (k + rGau)*dr;
        rhoD  = rhoZ[k] + rGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        phii[k+1] = phi0;
    }
    
    for(int i=0; i<=ngrid_m; ++i)
    {
        sigmai = SimpsonIntegration(rhoZ,0,ngridm1,0,i);
        
        Psi[i] = 2*coe*i*(sigma0-sigmai)+coe*(phii[i]-sigma0*i);
        Psi[ngrid-i] = Psi[i];
    }
    
    delete [] rhoZ;
    delete [] phii;
}



void PoissonEquationTwo(float* Z,double** rho,double phi,double* Psi)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    int     ngridm1; //,i0
    double  R1,lGau,rGau,coe,rhoD,phi0,sigma0,sigmai;
    double* rhoZ;
    double* phii;
    
    short   hspecies;
    
    hspecies = nspecies - 1;
    
    coe = 2*Pi*BJ*dr;
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    ngridm1= ngrid_m + 1;
    
    //ngrid0 = ngrid;
    //if(ngrid%2 != 0) ngrid0 = ngrid - 1;
    
    rhoZ = new double[ngridm1+1]();
    phii = new double[ngrid_m+1]();
    
    for(int k=0; k<=ngridm1; ++k)
    {
        for(short j=0; j<hspecies; ++j)
        {
            if(Z[j] != 0) rhoZ[k] += Z[j]*rho[j][k];
        }
    }
    
    //sigma0 = RombergIntegration(rhoZ,0,ngridm1,0,ngrid_m);
    sigma0 = SimpsonIntegration(rhoZ,0,ngridm1,0,ngrid_m);
    sigma0 = 2*sigma0;
    
    phi0    = 0;
    phii[0] = 0;
    for(int k=0; k<ngrid_m; ++k)
    {
        R1    = (k + lGau)*dr;
        rhoD  = rhoZ[k] + lGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        R1    = (k + rGau)*dr;
        rhoD  = rhoZ[k] + rGau*(rhoZ[k+1] - rhoZ[k]);
        phi0 += rhoD*R1;
        
        phii[k+1] = phi0;
    }
    
    
    for(int i=0; i<=ngrid_m; ++i)
    {
        //sigmai = RombergIntegration(rhoZ,0,ngridm1,0,i);
        sigmai = SimpsonIntegration(rhoZ,0,ngridm1,0,i);
        
        Psi[i] = 2*coe*i*(sigma0-sigmai)+coe*(phii[i]-sigma0*i);
    }
    
    for(int i=0; i<=ngrid_m; ++i)
    {
        Psi[i]      = Psi[i] + phi;
        Psi[ngrid-i]= Psi[i];
    }
    
    delete [] rhoZ;
    delete [] phii;
}



//using finite difference method and chasing method
void PoissonEquation(float* Z,double** rho,double* Psi1,double* Psi)
{
    double *A, *B, *C, *E, *L, *R, *Y, *rhoZ;
    double dr2,bi,kappa2,coe;
    int    i1;
    
    short   hspecies;
    
    hspecies = nspecies - 1;
    
    
    A     = new double[ngrid+1] ();
    B     = new double[ngrid+1] ();
    C     = new double[ngrid+1] ();
    E     = new double[ngrid+1] ();
    L     = new double[ngrid+1] ();
    R     = new double[ngrid] ();
    Y     = new double[ngrid+1] ();
    rhoZ  = new double[ngrid+1]();
    
    
    A[0]   = 0.0;
    B[0]   = -1.0;
    C[0]   = 1.0;
    E[0]   = boundary;
    A[ngrid] = 1.0;
    B[ngrid] = -1.0;
    C[ngrid] = 0.0;
    E[ngrid] = boundary;
    kappa2   = kapa*kapa;
    
    for(i1=0; i1<=ngrid; ++i1)
    {
        Psi1[i1] = Psi[i1];
    }
    
    for(i1=0; i1<=ngrid; ++i1)
    {
        for(short j=0; j<hspecies; ++j)
        {
            rhoZ[i1] = rhoZ[i1] - Z[j]*rho[j][i1];
        }
    }
    
    dr2 = dr*dr;
    bi  = dr2*kappa2;
    coe = 4*Pi*BJ*dr2;
    for(i1=1; i1<ngrid; i1++)
    {
        A[i1] = 1.0;
        B[i1] = -2.0-bi;
        C[i1] = 1.0;
        E[i1] = rhoZ[i1]*coe - bi*Psi1[i1];
    }
    
    
    //chasing method
    L[0] = B[0];
    Y[0] = E[0]/L[0];
    for(i1=1; i1<=ngrid; i1++)
    {
        R[i1-1] = C[i1-1]/L[i1-1];
        L[i1]   = B[i1] - A[i1]*R[i1-1];
        Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
    }
    
    Psi[ngrid] = Y[ngrid];
    for(i1=(ngrid-1); i1>=0; i1--)
    {
        Psi[i1] = Y[i1] - R[i1]*Psi[i1+1];
    }
    
    
    
    
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] E;
    delete [] L;
    delete [] R;
    delete [] Y;
    delete [] rhoZ;
}


