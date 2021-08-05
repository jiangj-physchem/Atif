//***********calculating bulk chemical potential****************//
#include "clibrary.h"
#include "bulkchempotentialDFT.h"
#include "derivhardspherechain.h"
#include "derivelectrocorrel.h"
#include "chargeshell.h"
#include "inhomvandelwaal.h"
#include "constantnum.h"
extern double errTol;
extern double BJ; //The Bjerrum length
extern double dr;
extern double depthW;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern short  neutsys;
extern short  nspecies;

using namespace std;

void BulkChemPotentialDFT(double gama,double* rhoB,float* D,double* B,double* BB,float* Z,
                          double** pairEner,double** ATT,double*** Psi_IJ,double* mu,double& P_bulk)
{
    double*  muID;
    double*  muHC;
    double*  muAT;
    double*  muRe;
    double*  muSH;

    double   rhoBM1,rhoBM2;
    double   f_id,f_at,f_el,f_hc;
    double   p_id,p_el,p_at,p_hc,p_sh;
    short    nblocks,i2,hspecies;
    
    nblocks  = nb[0] + nb[1];
    hspecies = nspecies - 1;
    
 
    muID   = new double[hspecies]();
    muSH   = new double[hspecies]();
    muHC   = new double[hspecies]();
    muAT   = new double[hspecies]();
    muRe   = new double[hspecies]();
    
    
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
    
    //////////////////////////////////////////////////////////////////////////////
    if(rhoBM1 > 1.0e-12)
    {
        for(short i=0; i<nb[0]; ++i)
        {
            if(mb[0][i] > 0)
            {
                ATT[i][i]  = ((double) (mb[0][i]-1))/((double) mp[0]);
            }
            else
            {
                ATT[i][i]  = 0;
            }
            
            if(i < (nb[0]-1))
            {
                if(mb[0][i+1] > 0)
                {
                    ATT[i][i+1] = 1.0/((double) mp[0]);
                }
                else
                {
                    ATT[i][i+1] = 0;
                }
            }
        }
    }
    
    
    if(rhoBM2 > 1.0e-12)
    {
        for(short i=nb[0]; i<nblocks; ++i)
        {
            i2 = i - nb[0];
            if(mb[1][i2] > 0)
            {
                ATT[i][i]  = ((double) (mb[1][i2]-1))/((double) mp[1]);
            }
            else
            {
                ATT[i][i]  = 0;
            }
            
            if(i < (nblocks-1))
            {
                if(mb[1][i2+1] > 0)
                {
                    ATT[i][i+1] = 1.0/((double) mp[1]);
                }
                else
                {
                    ATT[i][i+1] = 0;
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////
    
    //Start: bulk chemical potential for ideal gas
    f_id   = 0;
    for(short i=0; i<hspecies; ++i)
    {
        muID[i] = -1E20;
        if((i < nb[0]) && (rhoBM1 >  1E-15)) //if(i==0 || i==1)
        {
            muID[i] = log(rhoBM1/((double) mp[0]))/((double) mp[0]);
            f_id = f_id + (rhoB[i]/((double) mp[0]))*(log(rhoBM1/((double) mp[0])) - 1.0);
        }
        else if((i >= nb[0]) && (i < nblocks) && (rhoBM2 >  1E-15)) //else if(i==2 || i==3)
        {
            muID[i] = log(rhoBM2/((double) mp[1]))/((double) mp[1]);
            f_id = f_id + (rhoB[i]/((double) mp[1]))*(log(rhoBM2/((double) mp[1])) - 1.0);
        }
        else if((i >= nblocks) && (rhoB[i] >  1E-15))
        {
            muID[i] = log(rhoB[i]);
            f_id = f_id + rhoB[i]*(log(rhoB[i]) - 1.0);
        }
    }
    
    f_el = 0;
    if(neutsys == 1)
    {
        ChargeShell(rhoB,BB,Z,Psi_IJ,muSH);
        
        //ChargeShellDirk(rhoB,D,Z,Psi_IJ,muSH);
        DerivElectroCorrel(rhoB,gama,D,Z,ATT,muRe,f_el);
        
        //DerivElectroCorrelDirk(rhoB,gama,D,Z,muRe,f_el);
    }
    DerivHardSphereChain(rhoB,D,Z,ATT,muHC,f_hc);
    InhomVanDelWaal(rhoB,D,pairEner,muAT,f_at);

    
    p_id = -f_id;
    p_hc = -f_hc;
    p_at = -f_at;
    p_el = -f_el;
    
    p_sh = 0;
    for(short i=0; i<hspecies; ++i)
    {
        //mu[i] = muID[i] + muHS[i] + muCH[i] + muAT[i];
        mu[i] = muID[i] + muHC[i] + muAT[i] + muRe[i] + muSH[i];
        
        p_id = p_id + muID[i]*rhoB[i];
        p_hc = p_hc + muHC[i]*rhoB[i];
        p_at = p_at + muAT[i]*rhoB[i];
        p_el = p_el + muRe[i]*rhoB[i];
        p_sh = p_sh + 0.5*muSH[i]*rhoB[i];
        //cout<<"i= "<<i<<" muHC="<<muHC[i]<<" muAT="<<muAT[i]<<" muRe="<<muRe[i]<<" muSH="<<muSH[i]<<endl;
    }
    P_bulk = p_id + p_hc + p_at + p_el +p_sh;

    /*
    double* muSH_t;
    double  temp1,temp2,temp3,temp4,temp5,temp22,temp32;
    
    muSH_t= new double[hspecies]();
    
    for(short i=0; i<hspecies; ++i)
    {
        temp1     = Z[i]*BJ*Pi;
        muSH_t[i] = 0;
        for(short j=0; j<hspecies; ++j)
        {
            temp2 = BB[i] + BB[j];
            temp3 = B[i] + B[j];
            temp22= (temp3-temp2)*(temp3-temp2)*(temp3-temp2);
            temp32= (temp3+3*temp2);
            
            temp4 = temp22*temp32/(12*B[i]*B[j]) -2*(B[i]*B[i]+B[j]*B[j])/3;
            temp5 = Z[j]*rhoB[j];
            
            muSH_t[i] += (temp5*temp4);
        }
        
        muSH_t[i] = muSH_t[i]*temp1;
    }
    delete [] muSH_t;
     */
    
    

    delete [] muID;
    delete [] muHC;
    delete [] muAT;
    delete [] muRe;
    delete [] muSH;

}



void BulkPureSolvent(double rho_s0,float DS,double pEner,double& mu_s0,double& p_s0)
{
    double   f_id,f_hs,f_at;
    double   p_id,p_hs,p_at;
    double   AN3,AN3V,AN3D,AN3S,AN31,LAN31,AN3I,muAT,muHS,muID,S,V;
    double   part1,part2,part3,part4,part5,part6,DHN3I;
    
    //double   lGau,rGau,rin,rip,rin1,R1;
    //double   rin2,rin3,rin4,SWD,SWD2,SIJ2,DD2,sIJ;
    //double   coe,coe3,temp,tempv;
    //int      kin,kip,kin1,kin2,kin3,kin4;
    //lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    //rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    
    
    //Start: bulk chemical potential for hard sphere contribution
    
    S    = DS*DS*Pi;
    V    = DS*DS*DS*Pi/6;
    AN3I = rho_s0*V;
    
    AN3 = AN3I;
    AN3D= AN3I/DS;
    AN3S= AN3I/S;
    AN3V= AN3I/V;
    
    
    AN31  = 1 - AN3;
    LAN31 = log(AN31);
    
    part1 = AN3V/AN31;
    part2 = 18*AN3D/AN31;
    part3 = 18*AN3S/AN31;
    part4 = 18*AN3D*AN3S/(AN31*AN31);
    part5 = 18*AN3D*AN3D*(LAN31/(AN3*AN3) + 1/(AN3*AN31*AN31))/Pi;
    part6 = 6*AN3D*AN3D*AN3D*(2*LAN31/(AN3*AN3*AN3) + (2-5*AN3+AN3*AN3)/(AN3*AN3*AN31*AN31*AN31))/Pi;
    
    
    //f_hs = AN1*AN2/AN31 + AN0*AN3/(AN31*AN31);
    f_hs = 18*AN3S*AN3D/AN31 - AN3V*LAN31 + 6*AN3D*AN3D*AN3D*(AN3+AN31*AN31*LAN31)/(Pi*AN3*AN3*AN31*AN31);
    
    //End: bulk chemical potential for hard sphere contribution

    
    DHN3I = part1 - LAN31/V + part2/S + part3/DS + part4 + part5/DS - part6;
    muHS  = DHN3I*V;
    
    
    //Start: bulk chemical potential for ideal gas
    f_id   = 0;
    if(rho_s0 >  1E-15)
    {
        muID = log(rho_s0);
        f_id = f_id + rho_s0*(log(rho_s0) - 1.0);
    }
    else
    {
        muID = -1E20;
    }
    
    muAT = Pi*4*rho_s0*pEner*DS*DS*DS*(depthW*depthW*depthW - 1)/3;
    f_at = muAT*rho_s0*0.5;
    
    p_id = -f_id;
    p_hs = -f_hs;
    p_at = -f_at;
    
    mu_s0 = muID + muHS + muAT;
    p_id = p_id + muID*rho_s0;
    p_hs = p_hs + muHS*rho_s0;
    p_at = p_at + muAT*rho_s0;
    
    p_s0 = p_id + p_hs + p_at;
    
    
}

/*
void BulkChemPotentialDFT(double gama,double* rhoB,float* D,float* Z,double** pairEner,
                          double** ATT,double*** Psi_IJ,double* mu,double& P_bulk)
{
    double*  AN0I;
    double*  AN1I;
    double*  AN2I;
    double*  AN3I;
    double*  muID;
    double*  muHC;
    double*  muAT;
    double*  muRe;
    double*  muSH;
    double*  tempm;
    double*  X;
    double*  gamaD;
    double*  D2;
    double** sIJ;
    double** Y;
    double** G;
    double** DDC;
    double*** DY;
    double   lGau,rGau,rin,rip,rin1,R1,C,tempv,AN0M1,AN0M2,rhoBM1,rhoBM2;
    double   rin2,rin3,rin4,SWD,SWD2,SIJ2,DD2;
    double   f_id,f_hs,f_ch,f_at,f_el,f_hc;
    double   p_id,p_el,p_at,p_hc;
    double   AN0,AN1,AN2,AN3,AN31,LAN31;
    double   DHN0I,DHN1I,DHN2I,DHN3I,gamma_s;
    double   temp1,pEner,coe,coe0,coe1,coe2,coe3,temp,temp2,temp3,temp4,DSIJ;
    short    i1,p_max,nblocks,i2,hspecies;
    int      kin,kip,kin1,kin2,kin3,kin4,iter;
    
    
    lGau= 0.2113248654052;   //(0.5-sqrt(3)/6)*dr
    rGau= 0.7886751345948;   //(0.5+sqrt(3)/6)*dr
    
    nblocks  = nb[0] + nb[1];
    hspecies = nspecies - 1;
    
    AN0I   = new double[hspecies]();
    AN1I   = new double[hspecies]();
    AN2I   = new double[hspecies]();
    AN3I   = new double[hspecies]();
    
    muID   = new double[hspecies]();
    muSH   = new double[hspecies]();
    muHC   = new double[hspecies]();
    muAT   = new double[hspecies]();
    muRe   = new double[hspecies]();
    gamaD  = new double[hspecies]();
    D2     = new double[hspecies]();
    X      = new double[hspecies]();
    tempm  = new double[2]();
    
    
    
    
    DDC    = new double*[hspecies]();
    sIJ    = new double*[hspecies]();
    Y      = new double*[nblocks]();
    G      = new double*[nblocks]();
    
    
    for(short i=0; i<hspecies; ++i)
    {
        DDC[i] = new double[4]();
        sIJ[i] = new double[hspecies]();
    }
    
    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
        G[i]    = new double[nblocks]();
        
    }
    
    DY  = new double**[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        DY[i]  = new double*[nblocks]();
        for(short j=0; j<nblocks; ++j)
        {
            DY[i][j]  = new double[4]();
        }
    }
    
    //Start: bulk chemical potential for hard sphere contribution
    for(short i=0; i<hspecies; ++i)
    {
        AN0I[i] = rhoB[i];
        AN1I[i] = rhoB[i]*D[i]/2;
        AN2I[i] = rhoB[i]*D[i]*D[i]*Pi;
        AN3I[i] = rhoB[i]*D[i]*D[i]*D[i]*Pi/6;
        gamaD[i]= 1.0/(1+gama*D[i]);
        D2[i]   = D[i]*D[i];
    }
    
    
    AN0M1 = 0.0;
    rhoBM1= 0.0;
    for(short i=0; i<nb[0]; ++i)
    {
        AN0M1  = AN0M1 + AN0I[i];
        rhoBM1 = rhoBM1 + rhoB[i];
    }
    AN0M2 = 0.0;
    rhoBM2= 0.0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        AN0M2  = AN0M2 + AN0I[i];
        rhoBM2 = rhoBM2 + rhoB[i];
    }
    
    //////////////////////////////////////////////////////////////////////////////
    if(rhoBM1 > 1.0e-12)
    {
        for(short i=0; i<nb[0]; ++i)
        {
            if(mb[0][i] > 0)
            {
                ATT[i][i]  = ((double) (mb[0][i]-1))/((double) mp[0]);
            }
            else
            {
                ATT[i][i]  = 0;
            }
            
            if(i < (nb[0]-1))
            {
                if(mb[0][i+1] > 0)
                {
                    ATT[i][i+1] = 1.0/((double) mp[0]);
                }
                else
                {
                    ATT[i][i+1] = 0;
                }
            }
        }
    }
    
    
    if(rhoBM2 > 1.0e-12)
    {
        for(short i=nb[0]; i<nblocks; ++i)
        {
            i2 = i - nb[0];
            if(mb[1][i2] > 0)
            {
                ATT[i][i]  = ((double) (mb[1][i2]-1))/((double) mp[1]);
            }
            else
            {
                ATT[i][i]  = 0;
            }
            
            if(i < (nblocks-1))
            {
                if(mb[1][i2+1] > 0)
                {
                    ATT[i][i+1] = 1.0/((double) mp[1]);
                }
                else
                {
                    ATT[i][i+1] = 0;
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////
    
    AN0 = 0.0;
    AN1 = 0.0;
    AN2 = 0.0;
    AN3 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        AN0 = AN0 + AN0I[i];
        AN1 = AN1 + AN1I[i];
        AN2 = AN2 + AN2I[i];
        AN3 = AN3 + AN3I[i];
    }
    
    AN31  = 1 - AN3;
    LAN31 = log(AN31);
    DHN0I = -LAN31;
    DHN1I = AN2/AN31;
    DHN2I = AN1/AN31 + AN2*AN2/(12*Pi)*(LAN31/(AN3*AN3) + 1/(AN31*AN31*AN3));
    DHN3I = AN0/AN31 + AN1*AN2/(AN31*AN31) - AN2*AN2*AN2/(36*Pi)*(2*LAN31/(AN3*AN3*AN3)
                                                                  + (2-5*AN3+AN3*AN3)/(AN3*AN3*AN31*AN31*AN31));
    
    
    //f_hs = AN1*AN2/AN31 + AN0*AN3/(AN31*AN31);
    f_hs = AN1*AN2/AN31 - AN0*LAN31 + AN2*AN2*AN2*(AN3+AN31*AN31*LAN31)/(36*Pi*AN3*AN3*AN31*AN31);
    //End: bulk chemical potential for hard sphere contribution
    
    //calculate the bulk gama
    if(neutsys == 1)
    {
        C = 0.5*Pi/AN31;
        temp1 = 0;
        temp2 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            temp1 += (rhoB[i]*D[i]*Z[i]*gamaD[i]);
            temp2 += (rhoB[i]*D[i]*D2[i]*gamaD[i]);
        }
        gamma_s = C*temp1/(1+C*temp2);
        
        //cout<<"gamma_s_bulk= "<<(2*gamma_s/gama)<<" temp1="<<(Pi*temp1/gama)<<" temp2="<<(Pi*temp2/6)<<endl;
        
        temp1 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            X[i] = (Z[i]-gamma_s*D2[i])*gamaD[i];
            
            temp1 += (rhoB[i]*Z[i]*(X[i]-Z[i])/D[i]);
            
            temp2   = (gamma_s*D[i]*((gamma_s*D2[i]-2*Z[i])*gamaD[i] - gamma_s*D2[i]/3.0));
            muRe[i] = BJ*(temp2 - gama*Z[i]*Z[i]*gamaD[i]);
        }
        temp1 = temp1*BJ;
        f_el  = temp1 + (gama*gama*gama)/(3.0*Pi);
        
        for(short i=0; i<hspecies; ++i)
        {
            temp2 = 0;
            for(short j=0; j<hspecies; ++j)
            {
                temp2 += (Z[j]*rhoB[j]*D2[j]);
                
            }
            temp2  = Pi*BJ*Z[i]*temp2/6.0;
            muSH[i]= -temp2;
        }
    }
    
    //Start: bulk chemical potential for chain connectivity
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            sIJ[i][j] = 0.5*(D[i] + D[j]);
        }
    }
    
    p_max = 0;
    if(rhoBM1 > 1.0e-12) p_max = p_max + 1;
    if(rhoBM2 > 1.0e-12) p_max = p_max + 1;
    
    //Here Y is log(Y) actually
    for(short p=0; p<p_max; ++p)
    {
        i1 = p*nb[1];
        for(short i=(p*nb[0]); i<(nb[0]+i1); ++i)
        {
            DSIJ = D[i]*D[i]/sIJ[i][i];
            G[i][i] = 1/AN31 + AN2*DSIJ/(4*AN31*AN31) + AN2*AN2*DSIJ*DSIJ/(72*AN31*AN31*AN31);
            Y[i][i] = log(G[i][i]);//log(g_hs)
            
            if(i < (nb[0]+i1-1))
            {
                DSIJ = D[i]*D[i+1]/sIJ[i][i+1];
                G[i][i+1] = 1/AN31 + AN2*DSIJ/(4*AN31*AN31) + AN2*AN2*DSIJ*DSIJ/(72*AN31*AN31*AN31);
                Y[i][i+1] = log(G[i][i+1]);//log(g_hs)
            }
            
        }
    }
    
    //The derivation of log(Y)
    for(short p=0; p<p_max; ++p)
    {
        i1 = p*nb[1];
        for(short i=(p*nb[0]); i<(nb[0]+i1); ++i)
        {
            G[i][i] = G[i][i]*AN31*AN31;
            DSIJ    = D[i]*D[i]/sIJ[i][i];
            
            DY[i][i][2] = (DSIJ*0.25 + AN2*DSIJ*DSIJ/(36*AN31))/G[i][i];
            DY[i][i][3] = (1 + AN2*DSIJ/(2*AN31) + AN2*AN2*DSIJ*DSIJ/(24*AN31*AN31))/G[i][i];
            
            if(i < (nb[0]+i1-1))
            {
                G[i][i+1] = G[i][i+1]*AN31*AN31;
                DSIJ      = D[i]*D[i+1]/sIJ[i][i+1];
                
                DY[i][i+1][2] = (DSIJ*0.25 + AN2*DSIJ*DSIJ/(36*AN31))/G[i][i+1];
                DY[i][i+1][3] = (1 + AN2*DSIJ/(2*AN31) + AN2*AN2*DSIJ*DSIJ/(24*AN31*AN31))/G[i][i+1];
            }
        }
    }
    
    //The derivation of excess free energy density for chain connectivity
    f_ch = 0;
    if(rhoBM1 > 1.0e-12)
    {
        for(short j=2; j<4; ++j)
        {
            tempm[j-2] = 0.0;
            for(short i0=0; i0<nb[0]; ++i0)
            {
                tempm[j-2] = tempm[j-2] - ATT[i0][i0]*DY[i0][i0][j];
                if(i0 < (nb[0]-1)) tempm[j-2] = tempm[j-2] - ATT[i0][i0+1]*DY[i0][i0+1][j];
            }
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
        
        
        tempm[0] = 0.0;
        for(short i0=0; i0<nb[0]; ++i0)
        {
            tempm[0] = tempm[0] - ATT[i0][i0]*Y[i0][i0];
            if(i0 < (nb[0]-1))
            {
                tempm[0] = tempm[0] - ATT[i0][i0+1]*Y[i0][i0+1];
            }
        }
        for(short i=0; i<nb[0]; ++i)
        {
            DDC[i][0] += tempm[0];
            f_ch      += AN0I[i]*tempm[0];
        }
    }
    
    
    if(rhoBM2  > 1.0e-12)
    {
        for(short j=2; j<4; ++j)
        {
            tempm[j-2] = 0.0;
            for(short i0=nb[0]; i0<nblocks; ++i0)
            {
                tempm[j-2] = tempm[j-2] - ATT[i0][i0]*DY[i0][i0][j];
                if(i0 < (nblocks-1)) tempm[j-2] = tempm[j-2] - ATT[i0][i0+1]*DY[i0][i0+1][j];
            }
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
        
        tempm[0] = 0.0;
        for(short i0=nb[0]; i0<nblocks; ++i0)
        {
            tempm[0] = tempm[0] - ATT[i0][i0]*Y[i0][i0];
            if(i0 < (nblocks-1))
            {
                tempm[0] = tempm[0] - ATT[i0][i0+1]*Y[i0][i0+1];
            }
        }
        for(short i=nb[0]; i<nblocks; ++i)
        {
            DDC[i][0] += tempm[0];
            f_ch      += AN0I[i]*tempm[0];
        }
    }
    
    
    for(short i=0; i<hspecies; ++i)
    {
        temp1 = DHN0I + 0.5*D[i]*DHN1I + Pi*D[i]*D[i]*DHN2I + DHN3I*Pi*D[i]*D[i]*D[i]/6;
        temp2 = DDC[i][0] + DDC[i][1]*0.5*D[i] + Pi*D2[i]*DDC[i][2] + Pi*D2[i]*D[i]*DDC[i][3]/6;
        
        muHC[i] = temp1 + temp2;
        
    }
    
    //cout<<"DHN0I: "<<DHN0I<<" "<<DHN1I<<" "<<DHN2I<<" "<<DHN3I<<" "<<endl;
    
    //Start: bulk chemical potential for ideal gas
    f_id   = 0;
    for(short i=0; i<hspecies; ++i)
    {
        muID[i] = -1E20;
        if((i < nb[0]) && (rhoBM1 >  1E-15)) //if(i==0 || i==1)
        {
            muID[i] = log(rhoBM1/((double) mp[0]))/((double) mp[0]);
            f_id = f_id + (rhoB[i]/((double) mp[0]))*(log(rhoBM1/((double) mp[0])) - 1.0);
        }
        else if((i >= nb[0]) && (i < nblocks) && (rhoBM2 >  1E-15)) //else if(i==2 || i==3)
        {
            muID[i] = log(rhoBM2/((double) mp[1]))/((double) mp[1]);
            f_id = f_id + (rhoB[i]/((double) mp[1]))*(log(rhoBM2/((double) mp[1])) - 1.0);
        }
        else if((i >= nblocks) && (rhoB[i] >  1E-15))
        {
            muID[i] = log(rhoB[i]);
            f_id = f_id + rhoB[i]*(log(rhoB[i]) - 1.0);
        }
    }
    
    
    //van del Waals
    //depthW= 1.2;
    coe  = Pi*dr*0.5;
    f_at = 0;
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
        
        f_at = f_at + muAT[j1]*rhoB[j1]*0.5;
    }
    
    //electrostatistic correlation
    
    
    f_hc = f_hs + f_ch;
    p_id = -f_id;
    p_hc = -f_hc;
    p_at = -f_at;
    p_el = -f_el;
    
    for(short i=0; i<hspecies; ++i)
    {
        //mu[i] = muID[i] + muHS[i] + muCH[i] + muAT[i];
        mu[i] = muID[i] + muHC[i] + muAT[i] + muRe[i] + muSH[i];
        
        p_id = p_id + muID[i]*rhoB[i];
        p_hc = p_hc + muHC[i]*rhoB[i];
        p_at = p_at + muAT[i]*rhoB[i];
        p_el = p_el + muRe[i]*rhoB[i];
        cout<<"i= "<<i<<" muHC="<<muHC[i]<<" muAT="<<muAT[i]<<" muRe="<<muRe[i]<<" muSH="<<muSH[i]<<endl;
    }
    P_bulk = p_id + p_hc + p_at + p_el;
    
    
    exit(0);
    delete [] AN0I;
    delete [] AN1I;
    delete [] AN2I;
    delete [] AN3I;
    
    delete [] muID;
    delete [] muHC;
    delete [] muAT;
    delete [] gamaD;
    delete [] D2;
    delete [] muRe;
    delete [] muSH;
    delete [] X;
    delete [] tempm;
    
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DDC[i];
        delete [] sIJ[i];
    }
    delete [] DDC;
    delete [] sIJ;
    
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] ATT[i];
        delete [] Y[i];
        delete [] G[i];
    }
    delete [] ATT;
    delete [] Y;
    delete [] G;
    
    
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
 */

 
