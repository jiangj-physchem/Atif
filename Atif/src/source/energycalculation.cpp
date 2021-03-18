//***********potential from electrostatistic correlation****************//
#include "clibrary.h"
#include "inhomvandelwaal.h"
#include "imagecharge.h"
#include "derivelectrocorrel.h"
#include "energyelectrochain.h"
#include "simpsonintegration.h"
#include "energyimagecharge.h"
#include "derivhardspherechain.h"
#include "energyhardspherechain.h"
#include "chargeshell.h"
#include "constantnum.h"
#include "energycalculation.h"

extern double dr;
extern double BJ;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int ngrid_b;
extern int DMAX;
extern short  nspecies;
extern short neutsys;

void EnergyCalculation(double sigma,double gammab,double f,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                       double* BB,double** pairEner,double** ATT,double** rho,double* Psi,double*** Psi_IJ,
                       double& Ener_tot)
{
    double*  f_im;
    double*  u_im;
    double*  f_TOT;
    double** DCH;
    double** DEC;
    double** DSH;
    double** VanDW;
    double** UUIM;
    
    int    MAXR;
    short  hspecies,nblocks;
    double EN_EX,f_im_k,EN_HC,EN_EC,EN_IM,temp,rhoBM1,rhoBM2;
    
    nblocks = nb[0] + nb[1];
    hspecies= nspecies - 1;
    MAXR    = DMAX + ngrid_m;
    
    u_im   = new double[ngrid+1]();
    f_im   = new double[ngrid+1]();
    f_TOT  = new double[ngrid+1]();
    
    DCH  = new double*[hspecies]();
    DEC  = new double*[hspecies]();
    DSH  = new double*[hspecies]();
    VanDW= new double*[hspecies]();
    UUIM = new double*[hspecies]();
    
    for(short i=0; i<hspecies; ++i)
    {
        DCH[i]  = new double[ngrid+1]();
        DEC[i]  = new double[ngrid+1]();
        DSH[i]  = new double[ngrid+1]();
        VanDW[i]= new double[ngrid+1]();
        UUIM[i] = new double[ngrid+1]();
        
    }
    
    if(neutsys == 1)
    {
        if(f != 0.0)
        {
            ImageChargePotential(f,Z,rho,u_im);
            for(int k=0; k<=ngrid_m; ++k)
            {
                ImageChargeEnergy(k,f,Z,rho,f_im_k);
                f_im[k] = f_im_k;
            }
        }
        //electrostatic correlation from MSA
        DerivElectroCorrel(gammab,LLI,ULI,D,Z,rhoB,ATT,rho,DEC);
        //electrostatic contributions from charge shell
        ChargeShell(LLI,ULI,BB,rho,Psi_IJ,DSH);
    }
    
    
    //Square-well potential
    InhomVanDelWaal(LLI,ULI,D,rho,pairEner,VanDW);
    
    //hard sphere + hard sphere chain
    DerivHardSphereChain(MAXR,LLI,ULI,rhoB,D,Z,ATT,rho,DCH);
    
    
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        f_im[ngrid-k] = f_im[k];
        for(short j=0; j<hspecies; ++j)
        {
            DCH[j][ngrid-k] = DCH[j][k];
            DEC[j][ngrid-k] = DEC[j][k];
            DSH[j][ngrid-k] = DSH[j][k];
            
            VanDW[j][ngrid-k] = VanDW[j][k];
            
            UUIM[j][k] = Z[j]*Z[j]*u_im[k];
            UUIM[j][ngrid-k] = UUIM[j][k];
        }
    }
    
    
    EN_EC = 0;
    EN_HC = 0;
    EN_IM = 0;
    if(neutsys == 1)
    {
        if(f != 0.0) EN_IM = SimpsonIntegration(f_im,0,ngrid,0,ngrid);
        EnergyElectroChain(gammab,LLI,ULI,D,Z,rhoB,ATT,rho,EN_EC);
    }
    EnergyHardSphereChain(MAXR,LLI,ULI,rhoB,D,Z,ATT,rho,EN_HC);
    
    
    
    
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
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        temp = 0;
        for(short j=nblocks; j<hspecies; ++j)
        {
            temp += (DCH[j][k]+DEC[j][k]+UUIM[j][k]+(VanDW[j][k]+Psi[k]*Z[j]+DSH[j][k])*0.5
                     + 1.0)*rho[j][k];
        }
        
        if(rhoBM1 > 1E-10)
        {
            for(short j=0; j<nb[0]; ++j)
            {
                temp += (DCH[j][k]+DEC[j][k]+UUIM[j][k]+(VanDW[j][k]+Psi[k]*Z[j]+DSH[j][k])*0.5
                         + 1.0/((double) mp[0]))*rho[j][k];
            }
        }
        
        if(rhoBM2 > 1E-10)
        {
            for(short j=nb[0]; j<nblocks; ++j)
            {
                temp += (DCH[j][k]+DEC[j][k]+UUIM[j][k]+(VanDW[j][k]+Psi[k]*Z[j]+DSH[j][k])*0.5
                         + 1.0/((double) mp[1]))*rho[j][k];
            }
        }
        
        
        f_TOT[k]       = temp;
        f_TOT[ngrid-k] = temp;
    }
    
    EN_EX = SimpsonIntegration(f_TOT,0,ngrid,0,ngrid);
    EN_EX = 0.5*sigma*(Psi[0]+Psi[ngrid]) - EN_EX;
    
    //energy calculation end: Part 2
    Ener_tot = EN_EC + EN_HC + EN_IM + EN_EX;
    
    
    delete [] u_im;
    delete [] f_im;
    delete [] f_TOT;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] UUIM[i];
        delete [] DCH[i];
        delete [] DSH[i];
        delete [] DEC[i];
        delete [] VanDW[i];
    }
    delete [] UUIM;
    delete [] DCH;
    delete [] DEC;
    delete [] DSH;
    delete [] VanDW;
}



void EnergyCalculation(double sigma,double f,double eta,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                       double** pairEner,double** rho,double* Psi,double& Ener_tot)
{
    double*  f_im;
    double*  lambda;
    double*  u_im;
    double*  f_TOT;
    double** VanDW;
    double** UUIM;
    
    short  hspecies,nblocks;
    double EN_EX,f_im_k,EN_VV,EN_IM,temp,rhoBM1,rhoBM2,VV;
    
    nblocks = nb[0] + nb[1];
    hspecies= nspecies - 1;
    
    u_im   = new double[ngrid+1]();
    f_im   = new double[ngrid+1]();
    f_TOT  = new double[ngrid+1]();
    lambda = new double[ngrid+1]();
    VanDW= new double*[hspecies]();
    UUIM = new double*[hspecies]();
    
    for(short i=0; i<hspecies; ++i)
    {
        VanDW[i]= new double[ngrid+1]();
        UUIM[i] = new double[ngrid+1]();
        
    }
    
    if(neutsys == 1)
    {
        if(f != 0.0)
        {
            ImageChargePotential(f,Z,rho,u_im);
            for(int k=0; k<=ngrid_m; ++k)
            {
                ImageChargeEnergy(k,f,Z,rho,f_im_k);
                f_im[k] = f_im_k;
            }
        }
    }
    
    //Square-well potential
    InhomVanDelWaal(LLI,ULI,D,rho,pairEner,VanDW);
    
    VV = 6.0/(D[hspecies]*D[hspecies]*D[hspecies]*Pi);
    for(int k=0; k<=ngrid_m; ++k)
    {
        lambda[k] = log(rho[hspecies][k])*VV;
        lambda[ngrid-k] = lambda[k];
        f_im[ngrid-k]   = f_im[k];
        
        for(short j=0; j<hspecies; ++j)
        {
            VanDW[j][ngrid-k] = VanDW[j][k];
            
            UUIM[j][k] = Z[j]*Z[j]*u_im[k];
            UUIM[j][ngrid-k] = UUIM[j][k];
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
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        temp = 0;
        for(short j=nblocks; j<hspecies; ++j)
        {
            temp += (UUIM[j][k]+(VanDW[j][k]+Psi[k]*Z[j])*0.5 + 1.0)*rho[j][k];
        }
        
        if(rhoBM1 > 1E-10)
        {
            for(short j=0; j<nb[0]; ++j)
            {
                temp += (UUIM[j][k]+(VanDW[j][k]+Psi[k]*Z[j])*0.5 + 1.0/((double) mp[0]))*rho[j][k];
            }
        }
        
        if(rhoBM2 > 1E-10)
        {
            for(short j=nb[0]; j<nblocks; ++j)
            {
                temp += (UUIM[j][k]+(VanDW[j][k]+Psi[k]*Z[j])*0.5 + 1.0/((double) mp[1]))*rho[j][k];
            }
        }
        
        
        f_TOT[k]       = temp;
        f_TOT[ngrid-k] = temp;
    }
    
    EN_EX = SimpsonIntegration(f_TOT,0,ngrid,0,ngrid);
    EN_EX = 0.5*sigma*(Psi[0]+Psi[ngrid]) - EN_EX;
    
    
    EN_IM = 0;
    EN_VV = 0;
    if(neutsys == 1)
    {
        if(f != 0.0) EN_IM = SimpsonIntegration(f_im,0,ngrid,0,ngrid);
    }

    EN_VV = SimpsonIntegration(lambda,0,ngrid,0,ngrid);
    EN_VV = eta*EN_VV;
    
    //energy calculation end: Part 2
    Ener_tot = EN_VV + EN_IM + EN_EX;
    
    
    delete [] u_im;
    delete [] f_im;
    delete [] f_TOT;
    delete [] lambda;
    for(short i=0; i<hspecies; ++i)
    {
        delete [] UUIM[i];
        delete [] VanDW[i];
    }
    delete [] UUIM;
    delete [] VanDW;
}
