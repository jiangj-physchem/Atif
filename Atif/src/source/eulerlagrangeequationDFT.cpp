//*******************solve Euler Lagrange***************//
#include "eulerlagrangeequationDFT.h"
#include "derivhardspherechain.h"
#include "derivelectrocorrel.h"
#include "chargeshell.h"
#include "inhomvandelwaal.h"
#include "simpsonintegration.h"
#include "poissonequation.h"
#include "constantnum.h"
#include "volumefraction.h"
#include "renormalization.h"
#include "imagecharge.h"

#include "renormeulerlagrange.h"


extern double dr;
extern double errTol;
extern double C_psi;
extern int ngrid; //the number of grids
extern int ngrid_b; //the number of grids: the left bounday of bulk region
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern int ULIM; //the maxium uper limit of intergral
extern int LLIM; //the minimum lower limit of intergral
extern int DMAX;
extern short neutsys;
extern short nspecies;
extern double* rhoB_S;
extern double* temp_S;
extern double* DL1;
extern double* DR1;
extern double* DL2;
extern double* DR2;
extern double* GRC1;
extern double* GRC2;
extern float   ZMax;

void EulerLagrangeDFT(double gama,double phi,double f,double eta_t,int* LLI,int* ULI,float* D,double* BB,
                      float* Z,double* etar,double* Psi,double** pairEner,double* mu,double* rhoB,double** Ext,
                      double** ATT,double*** BesselZero,double*** Psi_IJ,double** rho,double** rho1,string* MODEL)
{
    double  R,rin,rip;
    double  ftemp1,ftemp2;//,coe;
    double  rhoBM1,rhoBM2,D_2,Dav_1,Dav_2; //thresh
    double* VI;
    double* temp;
    double* eta;
    double* u_im;
    double** DSH;
    double** DCH;
    double** DEC;
    double** VanDW;
    double** ff1;
    double** ff2;
    double** rho_s1;
    double** rho_s2;
    short*   MB1;
    short*   MB2;
    short mp1,mp2,nblocks,j0,i0,mp11,mp22,hspecies; //ithresh
    int kin,kip,MAXR1,igammar; //kip0
    
    struct InteLimit
    {
        int kin;    //lower limit
        int kip;    //upper limit
        
        InteLimit()
        {
            kin = 0;
            kip = 0;
        }
    };
    
    //coe = dr*0.5;
    mp1 = mp[0];
    mp2 = mp[1];
    hspecies = nspecies - 1;
    //thresh = 1.0E20;
    
    
    MAXR1 = DMAX + ngrid_m;
    if(neutsys == 1)
    {
        igammar = round(1.0/(dr*gama));
        MAXR1   = DMAX + igammar + ngrid_m;
    }
    if(MAXR1 > ngrid) MAXR1 = ngrid;
    
    nblocks= nb[0] + nb[1];

    VI   = new double[hspecies]();
    eta  = new double[hspecies]();
    u_im = new double[ngrid_m+1]();
    temp = new double[nblocks]();
    
    MB1  = new short[nb[0]+1]();
    MB2  = new short[nb[1]+1]();
    
    DCH  = new double*[nspecies]();
    DEC  = new double*[nspecies]();
    DSH  = new double*[nspecies]();
    VanDW= new double*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        VanDW[i]= new double[ngrid_m+1]();
        DCH[i]  = new double[ngrid_m+1]();
        DEC[i]  = new double[ngrid_m+1]();
        DSH[i]  = new double[ngrid_m+1]();
    }
    
    rhoBM1 = 0.0;
    MB1[0] = 0;
    Dav_1  = 0;
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1   = rhoBM1 + rhoB[i];
        MB1[i+1] = MB1[i] + mb[0][i];
        Dav_1   += D[i];
    }
    rhoBM2 = 0.0;
    MB2[0] = 0;
    Dav_2  = 0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        j0        = i - nb[0];
        rhoBM2    = rhoBM2 + rhoB[i];
        MB2[j0+1] = MB2[j0] + mb[1][j0];
        Dav_2    += D[i];
    }
    
    if(rhoBM1 > 1E-10)
    {
        ff1  = new double*[mp1]();
        rho_s1 = new double*[mp1]();
        for(short i=0; i<mp1; ++i)
        {
            ff1[i] = new double[ngrid+1]();
            rho_s1[i] = new double[ngrid+1]();
        }
    }
    
    if(rhoBM2 > 1E-10)
    {
        ff2  = new double*[mp2]();
        rho_s2 = new double*[mp2]();
        for(short i=0; i<mp2; ++i)
        {
            ff2[i] = new double[ngrid+1]();
            rho_s2[i] = new double[ngrid+1]();
        }
    }
    
    if(neutsys == 1)
    {
        //mean electric potential from Poisson equation
        
        if(ngrid_b < ngrid_m) PoissonEquationSingle(Z,rho,phi,0,Psi);
        if(ngrid_b== ngrid_m) PoissonEquationTwo(Z,rho,phi,Psi);
        if(f != 0.0) ImageChargePotential(f,Z,rho,u_im);
        
        //electrostatic correlation from MSA
        DerivElectroCorrel(gama,LLI,ULI,D,Z,rhoB,ATT,rho,DEC);
        //DerivElectroCorrelDirk(MAXR1,gama,LLI,ULI,D,Z,rho,DEC);
        
        //electrostatic contributions from charge shell
        ChargeShell(LLI,ULI,BB,rho,Psi_IJ,DSH);
        //ChargeShellDirk(LLI,ULI,D,rho,Psi_IJ,DSH);
    }
    
    
    //Square-well potential
    InhomVanDelWaal(LLI,ULI,D,rho,pairEner,VanDW);
    
    //hard sphere + hard sphere chain
    DerivHardSphereChain(MAXR1,LLI,ULI,rhoB,D,Z,ATT,rho,DCH);
    
    //loop-1
    for(int k=LLIM; k<=ngrid_m; ++k)
    {
        for(short j=nblocks; j<hspecies; ++j)
        {
            
            if((k>=LLI[j])&&(k<=ULI[j])&&(rhoB[j]>1E-10)&&(Ext[j][k]<1.0E10))
            {
                rho1[j][k] = exp(mu[j]-Ext[j][k]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                 -DCH[j][k]-Psi[k]*Z[j]-Z[j]*Z[j]*u_im[k]);
                if(rho1[j][k] > (1000*rhoB[j])) rho1[j][k] = rhoB[j]*1000;
            }
            else
            {
                rho1[j][k] = 0;
            }
            
            rho1[j][ngrid-k] = rho1[j][k];
        }
        
        
        if(rhoBM1 > 1E-10)
        {
            for(short j=0; j<nb[0]; ++j)
            {
                temp[j] = 0;
                if(Ext[j][k]<1.0E10)
                {
                    temp[j] = exp(mu[j]-Ext[j][k]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                  -DCH[j][k]-Psi[k]*Z[j]-Z[j]*Z[j]*u_im[k]);
                }
            }
            
            for(short i2=0; i2<nb[0]; ++i2)
            {
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    ff1[j2][k] = temp[i2];
                    //if(ff1[j2][k] > thresh) ff1[j2][k] = thresh;
                    ff1[j2][ngrid-k] = ff1[j2][k];
                }
            }
        }
        
        if(rhoBM2 > 1E-10)
        {
            for(short j=nb[0]; j<nblocks; ++j)
            {
                temp[j] = 0;
                if(Ext[j][k]<1.0E10)
                {
                    temp[j] = exp(mu[j]-Ext[j][k]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                  -DCH[j][k]-Psi[k]*Z[j]-Z[j]*Z[j]*u_im[k]);
                }
            }
            
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    ff2[j2][k] = temp[i0];
                    //if(ff2[j2][k] > thresh) ff2[j2][k] = thresh;
                    ff2[j2][ngrid-k] = ff2[j2][k];
                }
            }
        }
        
    }
    //loop-1
   
    
    //ithresh = 0;
    //loop-2
    if(rhoBM1 > 1E-10)
    {
        mp11 = mp1 - 1;
        if(MODEL[0] == "flexible")//flexible chain
        {
            double** fL1;
            double** fR1;
            fL1  = new double*[mp1]();
            fR1  = new double*[mp1]();
            for(short i=0; i<mp1; ++i)
            {
                fL1[i] = new double[ngrid+1]();
                fR1[i] = new double[ngrid+1]();
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                fL1[0][k]   = 1.0;
                fR1[mp11][k]= 1.0;
                
                fL1[0][ngrid-k]   = fL1[0][k];
                fR1[mp11][ngrid-k]= fR1[mp11][k];
            }
            
            for(short i1=1; i1<mp1; ++i1)
            {
                D_2 = 0.5/DL1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DL1[i1];
                    rip = R + DL1[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp1=RombergIntegration(ff1[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    ftemp1=SimpsonIntegration(ff1[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL1[i1][k] = ftemp1*D_2;
                    fL1[i1][ngrid-k] = fL1[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(short i1=(mp1-2); i1>=0; --i1)
            {
                D_2 = 0.5/DR1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DR1[i1];
                    rip = R + DR1[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp1=RombergIntegration(ff1[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    ftemp1=SimpsonIntegration(ff1[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR1[i1][k] = ftemp1*D_2;
                    fR1[i1][ngrid-k] = fR1[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    rho1[i2][k] = 0;
                    for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                    {
                        rho_s1[j2][k] = ff1[j2][k]*fL1[j2][k]*fR1[j2][k];
                        rho1[i2][k]   =  rho1[i2][k] + rho_s1[j2][k];
                    }
                    rho1[i2][ngrid-k] = rho1[i2][k];
                }
            }
            
            for(short i=0; i<mp1; ++i)
            {
                delete [] fL1[i];
                delete [] fR1[i];
            }
            delete [] fL1;
            delete [] fR1;
        }//flexible chain
        else if(MODEL[0] == "rod")//rod chain
        {
            int* N_z;
            double* DLL1;
            double* DRL1;
            double*** fL1;
            double*** fR1;
            double  slope_k,temp0;
            int     Nz2,Nz1,Nz0,Nzn,k2,k3,orient1,orient,ngNz0;
            
            InteLimit** inteBoud;
            Dav_1   = Dav_1/nb[0];
            orient1 = round(Dav_1/dr);
            orient  = 2*orient1;
            N_z  = new int[mp1]();
            DLL1 = new double[mp1]();
            DRL1 = new double[mp1]();
            fL1  = new double**[mp1]();
            fR1  = new double**[mp1]();

            inteBoud = new InteLimit*[mp1];
            for(short i=0; i<mp1; ++i)
            {
                fL1[i] = new double*[ngrid+1]();
                fR1[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL1[i][k] = new double[orient+1]();
                    fR1[i][k] = new double[orient+1]();
                }
                
                inteBoud[i] = new InteLimit[ngrid+1];
            }

            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(int k1=0; k1<=orient; ++k1)
                {
                    fL1[0][k][k1]    = 1.0;
                    fR1[mp11][k][k1] = 1.0;
                    
                    fL1[0][ngrid-k][k1]    = 1.0;
                    fR1[mp11][ngrid-k][k1] = 1.0;
                }
            }
            
            DLL1[0]   = 0;
            DRL1[mp11]= 0;
            for(short i1=1; i1<mp1; ++i1)
            {
                DLL1[i1]      =  DLL1[i1-1] + DL1[i1];
                DRL1[mp11-i1] =  DRL1[mp11-i1+1] + DR1[mp11-i1];
            }
            //orientation points to (the last bead - the first bead)
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {//loop for grids
                for(int k1=(-orient1); k1<=orient1; ++k1)
                {//loop for orientation
                    k2 = k1 + orient1;
                    k3 = orient1 - k1;

                    slope_k = k1/Dav_1;
                    for(short i1=0; i1<mp1; ++i1)
                    {
                        N_z[i1] = k + round(slope_k*DLL1[i1]);
                        if(N_z[i1] < 0) N_z[i1] = 0;
                        if(N_z[i1] > ngrid) N_z[i1] = ngrid;
                    }
                    //slope_kR = -slope_kL;
                    //Nzn = k - round(slope_kR*DRL1[0]);
                    Nz0  = N_z[0];
                    Nzn  = N_z[mp11];
                    ngNz0= ngrid - Nz0;
                    if((Nzn >= LLIM) && (Nzn <= ULIM))
                    {
                        temp0 = 1;
                        for(short i2=1; i2<mp1; ++i2)
                        {
                            Nz1    = N_z[i2-1];
                            Nz2    = N_z[i2];
                            
                            temp0 *= ff1[i2-1][Nz1];
                            fL1[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL1[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff1[i2+1][Nz2];
                            fR1[i2][Nz1][k2]       = temp0;
                            if(Nz0 != ngNz0) fR1[i2][ngrid-Nz1][k3] = temp0;
                        }
                        
                        if(k1 < 0)
                        {
                            for(short i2=0; i2<mp1; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kin += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kip += 1;
                            }
                        }
                        else if(k1 > 0)
                        {
                            for(short i2=0; i2<mp1; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kip += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kin += 1;
                            }
                        }
                    }
                    else
                    {
                        for(short i2=1; i2<mp1; ++i2)
                        {
                            Nz1 = N_z[i2];
                            fL1[i2][Nz1][k2]       = 0;
                            fL1[i2][ngrid-Nz1][k3] = 0;
                        }
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1 = N_z[i2];
                            fR1[i2][Nz1][k2]       = 0;
                            fR1[i2][ngrid-Nz1][k3] = 0;
                        }
                    }
                }//loop for orientation
                
            }//loop for grids
            
            D_2 = 0.5/Dav_1;
            for(short i1=0; i1<mp1; ++i1)
            {//loop for beads
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for grids
                    kin = orient1 - inteBoud[i1][k].kin;
                    kip = orient1 + inteBoud[i1][k].kip;
                    //ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],kin,kip,kin,kip);
                    ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],orient1,kin,kip);
                    rho_s1[i1][k] = ff1[i1][k]*ftemp1*D_2;
                }//loop for grids
            }//loop for beads
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    rho1[i2][k] = 0;
                    for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                    {
                        rho1[i2][k]   =  rho1[i2][k] + rho_s1[j2][k];
                    }
                    rho1[i2][ngrid-k] = rho1[i2][k];
                }
            }
            
            
            for(short i=0; i<mp1; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL1[i][k];
                    delete [] fR1[i][k];
                }
                delete [] fL1[i];
                delete [] fR1[i];
                delete [] inteBoud[i];
            }
            delete [] fL1;
            delete [] fR1;
            delete [] DLL1;
            delete [] DRL1;
            delete [] N_z;
            delete [] inteBoud;
            
        }//rod chain
        else if(MODEL[0] == "semiflexible" || MODEL[0] == "semi-flexible")//stiff chain
        {

            int*   nIntB;
            int*   delta_t;
            double*   f_0;
            double*   f_1;
            double*   f_2;
            double*** fL1;
            double*** fR1;
            double  nDR0,nDR1;
            int     nBessel,nBessel0,mp12,delta_Z1,delta_Z0,delta_Zt;
            int     kin1,kip1,midn,mIntB0,mIntB1,mIntBm,uppB,lowB;
            
            mp12    = mp11 - 1;
            nBessel = 2*round(1.0/dr) + 1;
            midn    = nBessel/2;
            nBessel0= midn*2;
            nIntB= new int[mp11]();
            f_0  = new double[nBessel]();
            f_1  = new double[nBessel]();
            f_2  = new double[nBessel]();
            fL1  = new double**[mp11]();
            fR1  = new double**[mp11]();
            delta_t = new int[MAXR1+1];
            for(short i=0; i<mp11; ++i)
            {
                nIntB[i] = round(DR1[i]/dr);
                fL1[i] = new double*[ngrid+1]();
                fR1[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL1[i][k] = new double[nBessel]();
                    fR1[i][k] = new double[nBessel]();
                }
            }
            
            
            mIntB0 = nIntB[0];
            mIntBm = nIntB[mp12];
            nDR0   = DR1[0];
            nDR1   = DR1[mp12];
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntB0;
                kip = k + mIntB0;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR0) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fL1[0][k][delta_Z0]    = 1.0;
                    fL1[0][ngrid-k][delta_Z1]    = 1.0;
                }
                
                
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR1) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fR1[mp12][k][delta_Z0] = 1.0;
                    fR1[mp12][ngrid-k][delta_Z1] = 1.0;
                }
            }
            
            for(short i1=1; i1<mp11; ++i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1-1];
                nDR0   = DR1[i1];
                nDR1   = DR1[i1-1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_1
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    kin1 = k - mIntB1;
                    kip1 = k + mIntB1;
                    if(kin1 < LLIM) kin1 = LLIM;
                    if(kip1 > ULIM) kip1 = ULIM;
                    lowB = round((k-kip1)/nDR1) + midn;
                    uppB = round((k-kin1)/nDR1) + midn;
                    for(int k2=kin1; k2<=kip1; ++k2)
                    {//loop for Z_0
                        delta_Z0 = round((k-k2)/nDR1) + midn;
                        f_0[delta_Z0]  = ff1[i1-1][k2];
                        f_2[delta_Z0]  = fL1[i1-1][k2][delta_Z0];
                        delta_t[k2]    = delta_Z0;
                    }//loop for Z_0
                    
                    
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z_2
                        //Nz2 = k + k1 - mIntB0;
                        delta_Z1 = round((k1 - k)/nDR0) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_0
                            delta_Z0 = delta_t[k2];
                            f_1[delta_Z0]  = f_0[delta_Z0]*BesselZero[0][delta_Z0][delta_Z1];
                        }//loop for Z_0
                        
                        
                        ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp1 = ftemp1*D_2; //without renormalization
                        ftemp1 = ftemp1/GRC1[delta_Z1];
                        fL1[i1][k][delta_Z1]= ftemp1;
                        
                        delta_Zt = nBessel0 - delta_Z1;
                        fL1[i1][ngrid-k][delta_Zt] = ftemp1;
                        
                        
                    }//loop for Z_2
                    
                }//loop for Z_1
                
            }//loop for beads
            
            for(short i1=(mp12-1); i1>=0; --i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1+1];
                nDR0   = DR1[i1];
                nDR1   = DR1[i1+1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_0
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //Nz1 = k;
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z1
                        delta_Z0 = round((k1 - k)/nDR0) + midn; //correct
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
                        
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            //delta_Z1 = delta_1[k2];
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR1[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff1[i1+2][k2]*BesselZero[0][delta_Z0][delta_Z1];
                        }//loop for Z_2
                        
                        ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp1 = ftemp1*D_2; //without renormalization
                        ftemp1 = ftemp1/GRC1[delta_Z0];
                        fR1[i1][k][delta_Z0]= ftemp1;
                        
                        delta_Zt = nBessel0 - delta_Z0;
                        fR1[i1][ngrid-k][delta_Zt] = ftemp1;
                    }//loop for Z_1
                }//loop for Z_0
                
            }//loop for beads
            
            for(short i1=0; i1<mp11; ++i1)
            {
                mIntB0 = nIntB[i1];
                nDR0   = DR1[i1];
                D_2    = 0.5/DR1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    lowB = round((kin-k)/nDR0) + midn;
                    uppB = round((kip-k)/nDR0) + midn;
                    for(int k1=kin; k1<=kip; ++k1)
                    {
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        f_1[delta_Z0]  = ff1[i1+1][k1];
                        f_2[delta_Z0]  = fL1[i1][k][delta_Z0]*fR1[i1][k][delta_Z0];
                    }
                    
                    ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s1[i1][k] = ff1[i1][k]*ftemp1*D_2;
                }//loop for Z
            }
            
            ///////////////////////////i1 = mp11////////////////////////////
            mIntBm = nIntB[mp12];
            nDR1   = DR1[mp12];
            D_2    = 0.5/DR1[mp12];;
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                lowB = round((k-kip)/nDR1) + midn;
                uppB = round((k-kin)/nDR1) + midn;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k-k1)/nDR1) + midn;
                    f_1[delta_Z0]  = ff1[mp12][k1];
                    f_2[delta_Z0]  = fL1[mp12][k1][delta_Z0];
                }
                
                ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s1[mp11][k] = ff1[mp11][k]*ftemp1*D_2;
            }
            ///////////////////////////////////////////////////////////////
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    rho1[i2][k] = 0;
                    for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                    {
                        rho1[i2][k]   =  rho1[i2][k] + rho_s1[j2][k];
                    }
                    rho1[i2][ngrid-k] = rho1[i2][k];
                }
            }
            
            
            
            for(short i=0; i<mp11; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL1[i][k];
                    delete [] fR1[i][k];
                }
                delete [] fL1[i];
                delete [] fR1[i];
            }
            delete [] fL1;
            delete [] fR1;
            delete [] nIntB;
            delete [] f_0;
            delete [] f_1;
            delete [] f_2;
            delete [] delta_t;

            
        }//stiff chain
        else
        {
            cerr<<"Error in eulerlagrangeequation.cpp: Polymer model is wrong"<<endl;
            exit(0);
        }
        
        
        
    }
    //loop-2

    //ithresh = 0;
    //loop-3
    if(rhoBM2 > 1E-10)
    {
        mp22 = mp2 - 1;
        if(MODEL[1] == "flexible")//flexible chain
        {
            double** fL2;
            double** fR2;
            fL2  = new double*[mp2]();
            fR2  = new double*[mp2]();
            for(short i=0; i<mp2; ++i)
            {
                fL2[i] = new double[ngrid+1]();
                fR2[i] = new double[ngrid+1]();
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                fL2[0][k]   = 1;
                fR2[mp22][k]= 1;
                
                fL2[0][ngrid-k]   = fL2[0][k];
                fR2[mp22][ngrid-k]= fR2[mp22][k];
            }
            
            for(short i1=1; i1<mp2; ++i1)
            {
                D_2 = 0.5/DL2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DL2[i1];
                    rip = R + DL2[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp2=RombergIntegration(ff2[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    ftemp2=SimpsonIntegration(ff2[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL2[i1][k] = ftemp2*D_2;
                    fL2[i1][ngrid-k] = fL2[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(short i1=(mp2-2); i1>=0; --i1)
            {
                D_2 = 0.5/DR2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DR2[i1];
                    rip = R + DR2[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp2=RombergIntegration(ff2[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    ftemp2=SimpsonIntegration(ff2[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR2[i1][k] = ftemp2*D_2;
                    fR2[i1][ngrid-k] = fR2[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    rho1[i0][k] = 0.0;
                    for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                    {
                        rho_s2[j2][k] = ff2[j2][k]*fL2[j2][k]*fR2[j2][k];
                        rho1[i0][k]   = rho1[i0][k] + rho_s2[j2][k];
                    }
                    rho1[i0][ngrid-k] = rho1[i0][k];
                }
            }
            
            
            for(short i=0; i<mp2; ++i)
            {
                delete [] fL2[i];
                delete [] fR2[i];
            }
            delete [] fL2;
            delete [] fR2;
        }//flexible chain
        else if(MODEL[1] == "rod")//rod chain
        {
            int* N_z;
            double* DLL2;
            double* DRL2;
            double*** fL2;
            double*** fR2;
            double  slope_k,temp0;
            int     Nz2,Nz1,Nz0,Nzn,k2,k3,orient2,orient,ngNz0;
            
            InteLimit** inteBoud;
            
            Dav_2   = Dav_2/nb[1];
            orient2 = round(Dav_2/dr);
            orient  = 2*orient2;
            N_z  = new int[mp2]();
            DLL2 = new double[mp2]();
            DRL2 = new double[mp2]();
            fL2  = new double**[mp2]();
            fR2  = new double**[mp2]();
            
            inteBoud = new InteLimit*[mp2];
            for(short i=0; i<mp2; ++i)
            {
                fL2[i] = new double*[ngrid+1]();
                fR2[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL2[i][k] = new double[orient+1]();
                    fR2[i][k] = new double[orient+1]();
                }
                
                inteBoud[i] = new InteLimit[ngrid+1];
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(int k1=0; k1<=orient; ++k1)
                {
                    fL2[0][k][k1]    = 1.0;
                    fR2[mp22][k][k1] = 1.0;
                    
                    fL2[0][ngrid-k][k1]    = 1.0;
                    fR2[mp22][ngrid-k][k1] = 1.0;
                }
            }
            
            DLL2[0]   = 0;
            DRL2[mp22]= 0;
            for(short i1=1; i1<mp2; ++i1)
            {
                DLL2[i1]      =  DLL2[i1-1] + DL2[i1];
                DRL2[mp22-i1] =  DRL2[mp22-i1+1] + DR2[mp22-i1];
            }
            
            //orientation points to (the last bead - the first bead)
            for(int k=LLIM; k<=ngrid_m; ++k)
            {//loop for grids
                for(int k1=(-orient2); k1<=orient2; ++k1)
                {//loop for orientation
                    k2 = k1 + orient2;
                    k3 = orient2 - k1;
                    
                    slope_k = k1/Dav_2;
                    for(short i1=0; i1<mp2; ++i1)
                    {
                        N_z[i1] = k + round(slope_k*DLL2[i1]);
                        if(N_z[i1] < 0) N_z[i1] = 0;
                        if(N_z[i1] > ngrid) N_z[i1] = ngrid;
                    }
                    
                    Nz0  = N_z[0];
                    Nzn  = N_z[mp22];
                    ngNz0= ngrid - Nz0;
                    if((Nzn >= LLIM) && (Nzn <= ULIM))
                    {
                        temp0 = 1;
                        for(short i2=1; i2<mp2; ++i2)
                        {
                            Nz1    = N_z[i2-1];
                            Nz2    = N_z[i2];
                            
                            temp0 *= ff2[i2-1][Nz1];
                            fL2[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL2[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff2[i2+1][Nz2];
                            fR2[i2][Nz1][k2]       = temp0;
                            if(Nz0 != ngNz0) fR2[i2][ngrid-Nz1][k3] = temp0;
                        }
                        
                        if(k1 < 0)
                        {
                            for(short i2=0; i2<mp2; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kin += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kip += 1;
                            }
                        }
                        else if(k1 > 0)
                        {
                            for(short i2=0; i2<mp2; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kip += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kin += 1;
                            }
                        }
                    }
                    else
                    {
                        for(short i2=1; i2<mp2; ++i2)
                        {
                            Nz1 = N_z[i2];
                            fL2[i2][Nz1][k2]       = 0;
                            fL2[i2][ngrid-Nz1][k3] = 0;
                        }
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1 = N_z[i2];
                            fR2[i2][Nz1][k2]       = 0;
                            fR2[i2][ngrid-Nz1][k3] = 0;
                        }
                        
                    }
                }//loop for orientation
                
            }//loop for grids
            
            
            D_2 = 0.5/Dav_2;
            for(short i1=0; i1<mp2; ++i1)
            {//loop for beads
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for grids
                    kin = orient2 - inteBoud[i1][k].kin;
                    kip = orient2 + inteBoud[i1][k].kip;
                    
                    ftemp2=SimpsonIntegration(fL2[i1][k],fR2[i1][k],orient2,kin,kip);
                    rho_s2[i1][k] = ff2[i1][k]*ftemp2*D_2;
                }//loop for grids
            }//loop for beads
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    rho1[i0][k] = 0.0;
                    for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                    {
                        rho1[i0][k]   = rho1[i0][k] + rho_s2[j2][k];
                    }
                    rho1[i0][ngrid-k] = rho1[i0][k];
                }
            }
            
            
            for(short i=0; i<mp2; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL2[i][k];
                    delete [] fR2[i][k];
                }
                delete [] fL2[i];
                delete [] fR2[i];
                delete [] inteBoud[i];
            }
            delete [] fL2;
            delete [] fR2;
            delete [] DLL2;
            delete [] DRL2;
            delete [] N_z;
            delete [] inteBoud;
        }//rod chain
        else if(MODEL[1] == "semiflexible" || MODEL[1] == "semi-flexible")//stiff chain
        {
            
            int*   nIntB;
            int*   delta_t;
            double*   f_0;
            double*   f_1;
            double*   f_2;
            double*** fL2;
            double*** fR2;
            double  nDR0,nDR1;
            int     nBessel,nBessel0,mp12,delta_Z1,delta_Z0,delta_Zt;
            int     kin1,kip1,midn,mIntB0,mIntB1,mIntBm,uppB,lowB;
            
            mp12    = mp22 - 1;
            nBessel = 2*round(1.0/dr) + 1;
            midn    = nBessel/2;
            nBessel0= midn*2;
            nIntB= new int[mp22]();
            f_0  = new double[nBessel]();
            f_1  = new double[nBessel]();
            f_2  = new double[nBessel]();
            fL2  = new double**[mp22]();
            fR2  = new double**[mp22]();
            delta_t = new int[MAXR1+1];
            for(short i=0; i<mp22; ++i)
            {
                nIntB[i] = round(DR2[i]/dr);
                fL2[i] = new double*[ngrid+1]();
                fR2[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL2[i][k] = new double[nBessel]();
                    fR2[i][k] = new double[nBessel]();
                }
            }
            
            mIntB0 = nIntB[0];
            mIntBm = nIntB[mp12];
            nDR0   = DR2[0];
            nDR1   = DR2[mp12];
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntB0;
                kip = k + mIntB0;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR0) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fL2[0][k][delta_Z0]    = 1.0;
                    fL2[0][ngrid-k][delta_Z1]    = 1.0;
                }
                
                
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR1) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fR2[mp12][k][delta_Z0] = 1.0;
                    fR2[mp12][ngrid-k][delta_Z1] = 1.0;
                }
            }
            
            for(short i1=1; i1<mp22; ++i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1-1];
                nDR0   = DR2[i1];
                nDR1   = DR2[i1-1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_1
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    kin1 = k - mIntB1;
                    kip1 = k + mIntB1;
                    if(kin1 < LLIM) kin1 = LLIM;
                    if(kip1 > ULIM) kip1 = ULIM;
                    
                    lowB = round((k-kip1)/nDR1) + midn;
                    uppB = round((k-kin1)/nDR1) + midn;
                    for(int k2=kin1; k2<=kip1; ++k2)
                    {//loop for Z_0
                        delta_Z0 = round((k-k2)/nDR1) + midn;
                        f_0[delta_Z0]  = ff2[i1-1][k2];
                        f_2[delta_Z0]  = fL2[i1-1][k2][delta_Z0];
                        delta_t[k2]    = delta_Z0;
                    }//loop for Z_0
                    
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z_2
                        //Nz2 = k + k1 - mIntB0;
                        delta_Z1 = round((k1 - k)/nDR0) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_0
                            delta_Z0 = delta_t[k2];
                            f_1[delta_Z0]  = f_0[delta_Z0]*BesselZero[1][delta_Z0][delta_Z1];
                        }//loop for Z_0
                        ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp2 = ftemp2*D_2;//without renormalization
                        ftemp2 = ftemp2/GRC2[delta_Z1];
                        fL2[i1][k][delta_Z1]= ftemp2;
                        
                        delta_Zt = nBessel0 - delta_Z1;
                        fL2[i1][ngrid-k][delta_Zt] = ftemp2;
                    }//loop for Z_2
                }//loop for Z_1
                
            }//loop for beads
            
            
            for(short i1=(mp12-1); i1>=0; --i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1+1];
                nDR0   = DR2[i1];
                nDR1   = DR2[i1+1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_0
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //Nz1 = k;
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z1
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
                        
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR2[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff2[i1+2][k2]*BesselZero[1][delta_Z0][delta_Z1];
                        }//loop for Z_2
                        
                        ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp2 = ftemp2*D_2;//without renormalization
                        ftemp2 = ftemp2/GRC2[delta_Z0];
                        fR2[i1][k][delta_Z0]= ftemp2;
                        
                        delta_Zt = nBessel0 - delta_Z0;
                        fR2[i1][ngrid-k][delta_Zt] = ftemp2;
                    }//loop for Z_1
                }//loop for Z_0
                
            }//loop for beads
            
            
            for(short i1=0; i1<mp22; ++i1)
            {
                mIntB0 = nIntB[i1];
                nDR0   = DR2[i1];
                D_2    = 0.5/DR2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    lowB = round((kin-k)/nDR0) + midn;
                    uppB = round((kip-k)/nDR0) + midn;
                    for(int k1=kin; k1<=kip; ++k1)
                    {
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        f_1[delta_Z0]  = ff2[i1+1][k1];
                        f_2[delta_Z0]  = fL2[i1][k][delta_Z0]*fR2[i1][k][delta_Z0];
                    }
                    
                    ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s2[i1][k] = ff2[i1][k]*ftemp2*D_2;
                }//loop for Z
            }
            
            ///////////////////////////i1 = mp22////////////////////////////
            mIntBm = nIntB[mp12];
            nDR1   = DR2[mp12];
            D_2    = 0.5/DR2[mp12];;
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                lowB = round((k-kip)/nDR1) + midn;
                uppB = round((k-kin)/nDR1) + midn;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k-k1)/nDR1) + midn;
                    f_1[delta_Z0]  = ff2[mp12][k1];
                    f_2[delta_Z0]  = fL2[mp12][k1][delta_Z0];
                }
                
                ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s2[mp22][k] = ff2[mp22][k]*ftemp2*D_2;
            }
            ///////////////////////////////////////////////////////////////
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    rho1[i0][k] = 0.0;
                    for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                    {
                        rho1[i0][k]   = rho1[i0][k] + rho_s2[j2][k];
                    }
                    rho1[i0][ngrid-k] = rho1[i0][k];
                }
            }
            
            
            
            
            
            for(short i=0; i<mp22; ++i)
            {
                for(int k=0; k<=ngrid; ++k)
                {
                    delete [] fL2[i][k];
                    delete [] fR2[i][k];
                }
                delete [] fL2[i];
                delete [] fR2[i];
            }
            delete [] fL2;
            delete [] fR2;
            delete [] nIntB;
            delete [] f_0;
            delete [] f_1;
            delete [] f_2;
            delete [] delta_t;
            
            
        }//stiff chain
        else
        {
            cerr<<"Error in eulerlagrangeequation.cpp: Polymer model is wrong"<<endl;
            exit(0);
        }
    }
    //loop-3
    
    VolumeFractionDFT(eta_t,D,LLI,ULI,etar,rho1);
    
    

  

//////////// normalize the densities such that the block length ratio is kept //////////////////
//    for(short i=0; i<4; ++i)
//    {
//        if(rhoB[i] > 1E-10)
//        {
//            rho_t[i] = 0;
//            for(int k=LLIM; k<ULIM; ++k)
//            {
//                R   = (k + lGau)*dr;
//                rhoD = rho1[i][k] + lGau*(rho1[i][k+1] - rho1[i][k]);
//                rho_t[i] = rho_t[i] + rhoD;
//
//                R   = (k + rGau)*dr;
//                rhoD = rho1[i][k] + rGau*(rho1[i][k+1] - rho1[i][k]);
//                rho_t[i] = rho_t[i] + rhoD;
//            }
//            rho_t[i] = rho_t[i]*coe;
//        }
//
//    }
    
//    for(int k=LLIM; k<=ULIM; ++k)
//    {
//        if((rhoB[0]>1E-10) && (rhoB[1]>1E-10))
//        {
//            rho1[0][k] = rho1[0][k]*rho_t[1]*ma[0]/((mp1-ma[0])*rho_t[0]);
//        }
//
//        if((rhoB[2]>1E-10) && (rhoB[3]>1E-10))
//        {
//            rho1[2][k] = rho1[2][k]*rho_t[3]*ma[1]/((mp2-ma[1])*rho_t[2]);
//        }
//    }

////////////////////////////////////////////////////////////////////////////////////////////////
    delete [] eta;
    delete [] temp;
    delete [] u_im;
    delete [] MB1;
    delete [] MB2;
    delete [] VI;
    
    
    if(rhoBM1 > 1E-10)
    {
        for(short i=0; i<mp1; ++i)
        {
            delete [] ff1[i];
            delete [] rho_s1[i];
        }
        delete [] ff1;
        delete [] rho_s1;
    }
    
    if(rhoBM2 > 1E-10)
    {
        for(short i=0; i<mp2; ++i)
        {
            delete [] ff2[i];
            delete [] rho_s2[i];
        }
        delete [] ff2;
        delete [] rho_s2;
    }
    
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] VanDW[i];
        delete [] DSH[i];
        delete [] DCH[i];
        delete [] DEC[i];
    }
    delete [] VanDW;
    delete [] DSH;
    delete [] DCH;
    delete [] DEC;
    
}




void EulerLagrangeDFT(double gama,double sigma,double f,double eta_t,double& deltaPhi,double err,int* LLI,int* ULI,float* D,
                      double* BB,double* etar,float* Z,double* Psi,double** pairEner,double* mu,double* rhoB,double** Ext,
                      double** ATT,double*** BesselZero,double*** Psi_IJ,double** rho,double** rho1,string* MODEL)
{
    
    double  R,rin,rip;
    double  ftemp1,ftemp2; //,coe;
    double  rhoBM1,rhoBM2,D_2,Dav_1,Dav_2; //thresh
    double* VI;
    double* eta;
    double* temp;
    double* u_im;
    double** DSH;
    double** DCH;
    double** DEC;
    double** VanDW;
    double** ff1;
    double** ff2;
    double** rho_s1;
    double** rho_s2;
    short*   MB1;
    short*   MB2;
    
    
    short mp1,mp2,nblocks,i0,j0,mp11,mp22,hspecies; //ithresh
    int kin,kip,MAXR1,igammar;//,kip0
    
    
    struct InteLimit
    {
        int kin;    //lower limit
        int kip;    //upper limit
        
        InteLimit()
        {
            kin = 0;
            kip = 0;
        }
    };
    
    

    nblocks= nb[0] + nb[1];
    hspecies = nspecies - 1;
    
    MAXR1   = DMAX + ngrid_m;
    if(neutsys == 1)
    {
        igammar = round(1.0/(dr*gama));
        MAXR1   = DMAX + igammar + ngrid_m;
    }
    
    if(MAXR1 > ngrid) MAXR1 = ngrid;
    
    mp1 = mp[0];
    mp2 = mp[1];
    
    //thresh = 1.0E3;
    
    VI   = new double[hspecies]();
    eta  = new double[hspecies]();
    MB1  = new short[nb[0]+1]();
    MB2  = new short[nb[1]+1]();
    temp = new double[nblocks]();
    u_im = new double[ngrid_m+1]();
    
    DCH  = new double*[nspecies]();
    DEC  = new double*[nspecies]();
    DSH  = new double*[nspecies]();
    VanDW= new double*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        VanDW[i]= new double[ngrid_m+1]();
        DCH[i]  = new double[ngrid_m+1]();
        DEC[i]  = new double[ngrid_m+1]();
        DSH[i]  = new double[ngrid_m+1]();
    }
    
    rhoBM1 = 0.0;
    MB1[0] = 0;
    Dav_1  = 0;
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1   = rhoBM1 + rhoB[i];
        MB1[i+1] = MB1[i] + mb[0][i];
        Dav_1   += D[i];
    }
    rhoBM2 = 0.0;
    MB2[0] = 0;
    Dav_2  = 0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        j0        = i - nb[0];
        rhoBM2    = rhoBM2 + rhoB[i];
        MB2[j0+1] = MB2[j0] + mb[1][j0];
        Dav_2    += D[i];
    }
    
    if(rhoBM1 > 1E-10)
    {
        ff1  = new double*[mp1]();
        rho_s1 = new double*[mp1]();
        for(short i=0; i<mp1; ++i)
        {
            ff1[i] = new double[ngrid+1]();
            rho_s1[i] = new double[ngrid+1]();
        }
    }
    
    if(rhoBM2 > 1E-10)
    {
        ff2  = new double*[mp2]();
        rho_s2 = new double*[mp2]();
        for(short i=0; i<mp2; ++i)
        {
            ff2[i] = new double[ngrid+1]();
            rho_s2[i] = new double[ngrid+1]();
        }
    }

    if(neutsys == 1)
    {
        if(f != 0.0) ImageChargePotential(f,Z,rho,u_im);
        
        //mean electric potential from Poisson equation
        if(ngrid_b < ngrid_m) PoissonEquationSingle(Z,rho,Psi);
        if(ngrid_b == ngrid_m)
        {
            PoissonEquationTwo(Z,rho,Psi);
            deltaPhi += (log(C_psi)/ZMax);
            for(int k=0; k<=ngrid_b; ++k)
            {
                Psi[k] = Psi[k] - deltaPhi;
                Psi[ngrid-k] = Psi[k];
            }
        }
        
        //electrostatic correlation from MSA
        DerivElectroCorrel(gama,LLI,ULI,D,Z,rhoB,ATT,rho,DEC);
        //DerivElectroCorrelDirk(MAXR1,gama,LLI,ULI,D,Z,rho,DEC);
        //electrostatic contributions from charge shell
        ChargeShell(LLI,ULI,BB,rho,Psi_IJ,DSH);
        //ChargeShellDirk(LLI,ULI,D,rho,Psi_IJ,DSH);
    }
    
    //Square-well potential
    InhomVanDelWaal(LLI,ULI,D,rho,pairEner,VanDW);
    //hard sphere + hard sphere chain
    DerivHardSphereChain(MAXR1,LLI,ULI,rhoB,D,Z,ATT,rho,DCH);
    //loop-1
    
    for(int k=LLIM; k<=ngrid_m; ++k)
    {
        for(short j=nblocks; j<hspecies; ++j)
        {
            //if((k>=LLI[j])&&(k<=ULI[j])&&(rhoB[j]>1E-10)&&(Ext[j][k]<1.0E10))
            if((Ext[j][k]<1.0E10) && (k<=ngrid_b))
            {

                rho1[j][k] = exp(mu[j]-Ext[j][k]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                 -DCH[j][k]-Psi[k]*Z[j]-Z[j]*Z[j]*u_im[k]);
                
                if(rho1[j][k] > (1000*rhoB[j])) rho1[j][k] = rhoB[j]*1000;
            }
            else if(k > ngrid_b)
            {
                rho1[j][k] = rhoB_S[j];
            }
            else
            {
                rho1[j][k] = 0;
            }
            rho1[j][ngrid-k] = rho1[j][k];
        }
        
        
        if(rhoBM1 > 1E-10)
        {
            for(short j=0; j<nb[0]; ++j)
            {
                temp[j] = 0;
                if((Ext[j][k]<1.0E10) && (k<=ngrid_b))
                {
                    temp[j] = exp(mu[j]-Ext[j][k]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                  -DCH[j][k]-Psi[k]*Z[j]-Z[j]*Z[j]*u_im[k]);
                }
                else if(k > ngrid_b)
                {
                    temp[j] = temp_S[j];
                }
            }
            
            
            for(short i2=0; i2<nb[0]; ++i2)
            {
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    ff1[j2][k] = temp[i2];
                    //if(ff1[j2][k] > thresh) ff1[j2][k] = thresh;
                    ff1[j2][ngrid-k] = ff1[j2][k];
                }
            }

        }
        
        
        if(rhoBM2 > 1E-10)
        {
            for(short j=nb[0]; j<nblocks; ++j)
            {
                temp[j] = 0;
                if((Ext[j][k]<1.0E10) && (k<=ngrid_b))
                {
                    temp[j] = exp(mu[j]-Ext[j][k]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                  -DCH[j][k]-Psi[k]*Z[j]-Z[j]*Z[j]*u_im[k]);
                }
                else if(k > ngrid_b)
                {
                    temp[j] = temp_S[j];
                }
            }
            
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    ff2[j2][k] = temp[i0];
                    //if(ff2[j2][k] > thresh) ff2[j2][k] = thresh;
                    ff2[j2][ngrid-k] = ff2[j2][k];
                }
            }
        }
        
    }
    //loop-1
    //exit(0);
    
    //ithresh = 0;
    //loop-2
    
    if(rhoBM1 > 1E-10)
    {
        mp11 = mp1 - 1;
        if(MODEL[0] == "flexible")//flexible chain
        {
            
            double** fL1;
            double** fR1;
            fL1  = new double*[mp1]();
            fR1  = new double*[mp1]();
            for(short i=0; i<mp1; ++i)
            {
                fL1[i] = new double[ngrid+1]();
                fR1[i] = new double[ngrid+1]();
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                fL1[0][k]   = 1.0;
                fR1[mp11][k]= 1.0;
                
                fL1[0][ngrid-k]   = fL1[0][k];
                fR1[mp11][ngrid-k]= fR1[mp11][k];
            }
            
            for(short i1=1; i1<mp1; ++i1)
            {
                D_2 = 0.5/DL1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DL1[i1];
                    rip = R + DL1[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp1=RombergIntegration(ff1[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    ftemp1=SimpsonIntegration(ff1[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL1[i1][k] = ftemp1*D_2;
                    fL1[i1][ngrid-k] = fL1[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            
            for(short i1=(mp1-2); i1>=0; --i1)
            {
                D_2 = 0.5/DR1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DR1[i1];
                    rip = R + DR1[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp1=RombergIntegration(ff1[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    ftemp1=SimpsonIntegration(ff1[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR1[i1][k] = ftemp1*D_2;
                    fR1[i1][ngrid-k] = fR1[i1][k];
                }
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    if(k<=ngrid_b)
                    {
                        rho1[i2][k] = 0;
                        for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                        {
                            rho_s1[j2][k] =  ff1[j2][k]*fL1[j2][k]*fR1[j2][k];
                            rho1[i2][k]   =  rho1[i2][k] + rho_s1[j2][k];
                        }
                    }
                    else
                    {
                        rho1[i2][k]   = rhoB_S[i2];
                    }
                    
                    rho1[i2][ngrid-k] = rho1[i2][k];
                }
            }
            
            
            
            for(short i=0; i<mp1; ++i)
            {
                delete [] fL1[i];
                delete [] fR1[i];
            }
            delete [] fL1;
            delete [] fR1;
        }//flexible chain
        else if(MODEL[0] == "rod")//rod chain
        {
            
            int* N_z;
            double* DLL1;
            double* DRL1;
            double*** fL1;
            double*** fR1;
            double  slope_k,temp0;
            int     Nz2,Nz1,Nz0,Nzn,k2,k3,orient1,orient,ngNz0;
            
            InteLimit** inteBoud;
            
            Dav_1   = Dav_1/nb[0];
            orient1 = round(Dav_1/dr);
            orient  = 2*orient1;
            N_z  = new int[mp1]();
            DLL1 = new double[mp1]();
            DRL1 = new double[mp1]();
            fL1  = new double**[mp1]();
            fR1  = new double**[mp1]();
            
            inteBoud = new InteLimit*[mp1];
            for(short i=0; i<mp1; ++i)
            {
                fL1[i] = new double*[ngrid+1]();
                fR1[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL1[i][k] = new double[orient+1]();
                    fR1[i][k] = new double[orient+1]();
                }
                
                inteBoud[i] = new InteLimit[ngrid+1];
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(int k1=0; k1<=orient; ++k1)
                {
                    fL1[0][k][k1]    = 1.0;
                    fR1[mp11][k][k1] = 1.0;
                    
                    fL1[0][ngrid-k][k1]    = 1.0;
                    fR1[mp11][ngrid-k][k1] = 1.0;
                }
            }
            
            DLL1[0]   = 0;
            DRL1[mp11]= 0;
            for(short i1=1; i1<mp1; ++i1)
            {
                DLL1[i1]      =  DLL1[i1-1] + DL1[i1];
                DRL1[mp11-i1] =  DRL1[mp11-i1+1] + DR1[mp11-i1];
            }
            
            //orientation points to (the last bead - the first bead)
            for(int k=LLIM; k<=ngrid_m; ++k)
            {//loop for grids
                for(int k1=(-orient1); k1<=orient1; ++k1)
                {//loop for orientation
                    k2 = k1 + orient1;
                    k3 = orient1 - k1;
                    
                    slope_k = k1/Dav_1;
                    for(short i1=0; i1<mp1; ++i1)
                    {
                        N_z[i1] = k + round(slope_k*DLL1[i1]);
                        if(N_z[i1] < 0) N_z[i1] = 0;
                        if(N_z[i1] > ngrid) N_z[i1] = ngrid;
                    }
                    //slope_kR = -slope_kL;
                    //Nzn = k - round(slope_kR*DRL1[0]);
                    
                    Nz0  = N_z[0];
                    Nzn  = N_z[mp11];
                    ngNz0= ngrid - Nz0;
                    if((Nzn >= LLIM) && (Nzn <= ULIM))
                    {
                        temp0 = 1;
                        for(short i2=1; i2<mp1; ++i2)
                        {
                            Nz1    = N_z[i2-1];
                            Nz2    = N_z[i2];
                            
                            temp0 *= ff1[i2-1][Nz1];
                            fL1[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL1[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff1[i2+1][Nz2];
                            fR1[i2][Nz1][k2]       = temp0;
                            if(Nz0 != ngNz0) fR1[i2][ngrid-Nz1][k3] = temp0;
                        }
                        
                        if(k1 < 0)
                        {
                            for(short i2=0; i2<mp1; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kin += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kip += 1;
                            }
                        }
                        else if(k1 > 0)
                        {
                            for(short i2=0; i2<mp1; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kip += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kin += 1;
                            }
                        }
                    }
                    else
                    {
                        for(short i2=1; i2<mp1; ++i2)
                        {
                            Nz1 = N_z[i2];
                            fL1[i2][Nz1][k2]       = 0;
                            fL1[i2][ngrid-Nz1][k3] = 0;
                        }
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1 = N_z[i2];
                            fR1[i2][Nz1][k2]       = 0;
                            fR1[i2][ngrid-Nz1][k3] = 0;
                        }
                        
                    }
                }//loop for orientation
                
            }//loop for grids
            
            
            D_2 = 0.5/Dav_1;
            for(short i1=0; i1<mp1; ++i1)
            {//loop for beads
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for grids
                    kin = orient1 - inteBoud[i1][k].kin;
                    kip = orient1 + inteBoud[i1][k].kip;
                    
                    //ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],kin,kip,kin,kip);
                    ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],orient1,kin,kip);
                    rho_s1[i1][k] = ff1[i1][k]*ftemp1*D_2;
                }//loop for grids
            }//loop for beads
            
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    if(k<=ngrid_b)
                    {
                        rho1[i2][k] = 0;
                        for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                        {
                            rho1[i2][k]   =  rho1[i2][k] + rho_s1[j2][k];
                        }
                    }
                    else
                    {
                        rho1[i2][k]   = rhoB_S[i2];
                    }
                    rho1[i2][ngrid-k] = rho1[i2][k];
                    
                }
            }
            
            
            for(short i=0; i<mp1; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL1[i][k];
                    delete [] fR1[i][k];
                }
                delete [] fL1[i];
                delete [] fR1[i];
                delete [] inteBoud[i];
            }
            delete [] fL1;
            delete [] fR1;
            delete [] DLL1;
            delete [] DRL1;
            delete [] N_z;
            delete [] inteBoud;
        }//rod chain
        else if(MODEL[0] == "semiflexible" || MODEL[0] == "semi-flexible")//stiff chain
        {
            
            int*   nIntB;
            int*   delta_t;
            double*   f_0;
            double*   f_1;
            double*   f_2;
            double*** fL1;
            double*** fR1;
            double  nDR0,nDR1;
            int     nBessel,nBessel0,mp12,delta_Z1,delta_Z0,delta_Zt;
            int     kin1,kip1,midn,mIntB0,mIntB1,mIntBm,uppB,lowB;
            
            mp12    = mp11 - 1;
            nBessel = 2*round(1.0/dr) + 1;
            midn    = nBessel/2;
            nBessel0= midn*2;
            nIntB= new int[mp11]();
            f_0  = new double[nBessel]();
            f_1  = new double[nBessel]();
            f_2  = new double[nBessel]();
            fL1  = new double**[mp11]();
            fR1  = new double**[mp11]();
            delta_t = new int[MAXR1+1];
            for(short i=0; i<mp11; ++i)
            {
                nIntB[i] = round(DR1[i]/dr);
                fL1[i] = new double*[ngrid+1]();
                fR1[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL1[i][k] = new double[nBessel]();
                    fR1[i][k] = new double[nBessel]();
                }
            }
            
            
            mIntB0 = nIntB[0];
            mIntBm = nIntB[mp12];
            nDR0   = DR1[0];
            nDR1   = DR1[mp12];
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntB0;
                kip = k + mIntB0;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR0) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fL1[0][k][delta_Z0]    = 1.0;
                    fL1[0][ngrid-k][delta_Z1]    = 1.0;
                }
                
                
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR1) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fR1[mp12][k][delta_Z0] = 1.0;
                    fR1[mp12][ngrid-k][delta_Z1] = 1.0;
                }
            }
            
            for(short i1=1; i1<mp11; ++i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1-1];
                nDR0   = DR1[i1];
                nDR1   = DR1[i1-1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_1
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    kin1 = k - mIntB1;
                    kip1 = k + mIntB1;
                    if(kin1 < LLIM) kin1 = LLIM;
                    if(kip1 > ULIM) kip1 = ULIM;
                    lowB = round((k-kip1)/nDR1) + midn;
                    uppB = round((k-kin1)/nDR1) + midn;
                    for(int k2=kin1; k2<=kip1; ++k2)
                    {//loop for Z_0
                        delta_Z0 = round((k-k2)/nDR1) + midn;
                        f_0[delta_Z0]  = ff1[i1-1][k2];
                        f_2[delta_Z0]  = fL1[i1-1][k2][delta_Z0];
                        delta_t[k2]    = delta_Z0;
                    }//loop for Z_0
                    
                    
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z_2
                        //Nz2 = k + k1 - mIntB0;
                        delta_Z1 = round((k1 - k)/nDR0) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_0
                            delta_Z0 = delta_t[k2];
                            f_1[delta_Z0]  = f_0[delta_Z0]*BesselZero[0][delta_Z0][delta_Z1];
                        }//loop for Z_0
                        
                        
                        ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp1 = ftemp1*D_2; //without renormalization
                        ftemp1 = ftemp1/GRC1[delta_Z1];
                        fL1[i1][k][delta_Z1]= ftemp1;
                        
                        delta_Zt = nBessel0 - delta_Z1;
                        fL1[i1][ngrid-k][delta_Zt] = ftemp1;
                        
                    }//loop for Z_2
                   
                }//loop for Z_1
                
            }//loop for beads
            
            for(short i1=(mp12-1); i1>=0; --i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1+1];
                nDR0   = DR1[i1];
                nDR1   = DR1[i1+1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_0
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //Nz1 = k;
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z1
                        delta_Z0 = round((k1 - k)/nDR0) + midn; //correct
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
            
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            //delta_Z1 = delta_1[k2];
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR1[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff1[i1+2][k2]*BesselZero[0][delta_Z0][delta_Z1];
                        }//loop for Z_2
                        
                        ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp1 = ftemp1*D_2; //without renormalization
                        ftemp1 = ftemp1/GRC1[delta_Z0];
                        fR1[i1][k][delta_Z0]= ftemp1;
                        
                        delta_Zt = nBessel0 - delta_Z0;
                        fR1[i1][ngrid-k][delta_Zt] = ftemp1;
                    }//loop for Z_1
                }//loop for Z_0
                
            }//loop for beads
            
            for(short i1=0; i1<mp11; ++i1)
            {
                mIntB0 = nIntB[i1];
                nDR0   = DR1[i1];
                D_2    = 0.5/DR1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    lowB = round((kin-k)/nDR0) + midn;
                    uppB = round((kip-k)/nDR0) + midn;
                    for(int k1=kin; k1<=kip; ++k1)
                    {
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        f_1[delta_Z0]  = ff1[i1+1][k1];
                        f_2[delta_Z0]  = fL1[i1][k][delta_Z0]*fR1[i1][k][delta_Z0];
                    }
                    
                    ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s1[i1][k] = ff1[i1][k]*ftemp1*D_2;
                }//loop for Z
            }
            
            ///////////////////////////i1 = mp11////////////////////////////
            mIntBm = nIntB[mp12];
            nDR1   = DR1[mp12];
            D_2    = 0.5/DR1[mp12];;
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                lowB = round((k-kip)/nDR1) + midn;
                uppB = round((k-kin)/nDR1) + midn;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k-k1)/nDR1) + midn;
                    f_1[delta_Z0]  = ff1[mp12][k1];
                    f_2[delta_Z0]  = fL1[mp12][k1][delta_Z0];
                }
                
                ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s1[mp11][k] = ff1[mp11][k]*ftemp1*D_2;
            }
            ///////////////////////////////////////////////////////////////
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    if(k<=ngrid_b)
                    {
                        rho1[i2][k] = 0;
                        for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                        {
                            rho1[i2][k]   =  rho1[i2][k] + rho_s1[j2][k];
                        }
                    }
                    else
                    {
                        rho1[i2][k]   = rhoB_S[i2];
                    }
                    
                    rho1[i2][ngrid-k] = rho1[i2][k];
                }
            }
            
            
            
            for(short i=0; i<mp11; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL1[i][k];
                    delete [] fR1[i][k];
                }
                delete [] fL1[i];
                delete [] fR1[i];
            }
            delete [] fL1;
            delete [] fR1;
            delete [] nIntB;
            delete [] f_0;
            delete [] f_1;
            delete [] f_2;
            delete [] delta_t;
            
            
        }//stiff chain
        else
        {
            cerr<<"Error in eulerlagrangeequation.cpp: Polymer model is wrong"<<endl;
            exit(0);
        }
        
        
        
    }
    //loop-2
    
    
    //ithresh = 0;
    //loop-3
    if(rhoBM2 > 1E-10)
    {
        mp22 = mp2 - 1;
        if(MODEL[1] == "flexible")//flexible chain
        {
            double** fL2;
            double** fR2;
            fL2  = new double*[mp2]();
            fR2  = new double*[mp2]();
            for(short i=0; i<mp2; ++i)
            {
                fL2[i] = new double[ngrid+1]();
                fR2[i] = new double[ngrid+1]();
            }
        
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                fL2[0][k]   = 1;
                fR2[mp22][k]= 1;
                
                fL2[0][ngrid-k]   = fL2[0][k];
                fR2[mp22][ngrid-k]= fR2[mp22][k];
            }
            
            for(short i1=1; i1<mp2; ++i1)
            {
                D_2 = 0.5/DL2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DL2[i1];
                    rip = R + DL2[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp2=RombergIntegration(ff2[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    ftemp2=SimpsonIntegration(ff2[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL2[i1][k] = ftemp2*D_2;
                    fL2[i1][ngrid-k] = fL2[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(short i1=(mp2-2); i1>=0; --i1)
            {
                D_2 = 0.5/DR2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DR2[i1];
                    rip = R + DR2[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp2=RombergIntegration(ff2[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    ftemp2=SimpsonIntegration(ff2[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR2[i1][k] = ftemp2*D_2;
                    fR2[i1][ngrid-k] = fR2[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    if(k<=ngrid_b)
                    {
                        rho1[i0][k] = 0.0;
                        for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                        {
                            rho_s2[j2][k] = ff2[j2][k]*fL2[j2][k]*fR2[j2][k];
                            rho1[i0][k]   = rho1[i0][k] + rho_s2[j2][k];
                        }
                    }
                    else
                    {
                        rho1[i0][k]   = rhoB_S[i0];
                    }
                    
                    rho1[i0][ngrid-k] = rho1[i0][k];
                }
            }
            
            for(short i=0; i<mp2; ++i)
            {
                delete [] fL2[i];
                delete [] fR2[i];
            }
            delete [] fL2;
            delete [] fR2;
        }//flexible chain
        else if(MODEL[1] == "rod")//rod chain
        {
            int* N_z;
            double* DLL2;
            double* DRL2;
            double*** fL2;
            double*** fR2;
            double  slope_k,temp0;
            int     Nz2,Nz1,Nz0,Nzn,k2,k3,orient2,orient,ngNz0;
            
            InteLimit** inteBoud;
            
            Dav_2   = Dav_2/nb[1];
            orient2 = round(Dav_2/dr);
            orient  = 2*orient2;
            N_z  = new int[mp2]();
            DLL2 = new double[mp2]();
            DRL2 = new double[mp2]();
            fL2  = new double**[mp2]();
            fR2  = new double**[mp2]();
            
            inteBoud = new InteLimit*[mp2];
            for(short i=0; i<mp2; ++i)
            {
                fL2[i] = new double*[ngrid+1]();
                fR2[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL2[i][k] = new double[orient+1]();
                    fR2[i][k] = new double[orient+1]();
                }
                
                inteBoud[i] = new InteLimit[ngrid+1];
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(int k1=0; k1<=orient; ++k1)
                {
                    fL2[0][k][k1]    = 1.0;
                    fR2[mp22][k][k1] = 1.0;
                    
                    fL2[0][ngrid-k][k1]    = 1.0;
                    fR2[mp22][ngrid-k][k1] = 1.0;
                }
            }
            
            DLL2[0]   = 0;
            DRL2[mp22]= 0;
            for(short i1=1; i1<mp2; ++i1)
            {
                DLL2[i1]      =  DLL2[i1-1] + DL2[i1];
                DRL2[mp22-i1] =  DRL2[mp22-i1+1] + DR2[mp22-i1];
            }
            
            //orientation points to (the last bead - the first bead)
            for(int k=LLIM; k<=ngrid_m; ++k)
            {//loop for grids
                for(int k1=(-orient2); k1<=orient2; ++k1)
                {//loop for orientation
                    k2 = k1 + orient2;
                    k3 = orient2 - k1;
                    
                    slope_k = k1/Dav_2;
                    for(short i1=0; i1<mp2; ++i1)
                    {
                        N_z[i1] = k + round(slope_k*DLL2[i1]);
                        if(N_z[i1] < 0) N_z[i1] = 0;
                        if(N_z[i1] > ngrid) N_z[i1] = ngrid;
                    }
                    //slope_kR = -slope_kL;
                    //Nzn = k - round(slope_kR*DRL1[0]);
                    
                    Nz0  = N_z[0];
                    Nzn  = N_z[mp22];
                    ngNz0= ngrid - Nz0;
                    if((Nzn >= LLIM) && (Nzn <= ULIM))
                    {
                        temp0 = 1;
                        for(short i2=1; i2<mp2; ++i2)
                        {
                            Nz1    = N_z[i2-1];
                            Nz2    = N_z[i2];
                            
                            temp0 *= ff2[i2-1][Nz1];
                            fL2[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL2[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff2[i2+1][Nz2];
                            fR2[i2][Nz1][k2]       = temp0;
                            if(Nz0 != ngNz0) fR2[i2][ngrid-Nz1][k3] = temp0;
                        }
                        
                        if(k1 < 0)
                        {
                            for(short i2=0; i2<mp2; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kin += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kip += 1;
                            }
                        }
                        else if(k1 > 0)
                        {
                            for(short i2=0; i2<mp2; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kip += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kin += 1;
                            }
                        }
                    }
                    else
                    {
                        for(short i2=1; i2<mp2; ++i2)
                        {
                            Nz1 = N_z[i2];
                            fL2[i2][Nz1][k2]       = 0;
                            fL2[i2][ngrid-Nz1][k3] = 0;
                        }
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1 = N_z[i2];
                            fR2[i2][Nz1][k2]       = 0;
                            fR2[i2][ngrid-Nz1][k3] = 0;
                        }
                        
                    }
                }//loop for orientation
                
            }//loop for grids
            
            
            D_2 = 0.5/Dav_2;
            for(short i1=0; i1<mp2; ++i1)
            {//loop for beads
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for grids
                    kin = orient2 - inteBoud[i1][k].kin;
                    kip = orient2 + inteBoud[i1][k].kip;
                    
                    ftemp2=SimpsonIntegration(fL2[i1][k],fR2[i1][k],orient2,kin,kip);
                    rho_s2[i1][k] = ff2[i1][k]*ftemp2*D_2;
                }//loop for grids
            }//loop for beads
            
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    if(k<=ngrid_b)
                    {
                        rho1[i0][k] = 0.0;
                        for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                        {
                            rho1[i0][k]   = rho1[i0][k] + rho_s2[j2][k];
                        }
                    }
                    else
                    {
                        rho1[i0][k]   = rhoB_S[i0];
                    }
                    
                    rho1[i0][ngrid-k] = rho1[i0][k];
                }
            }
            
            for(short i=0; i<mp2; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL2[i][k];
                    delete [] fR2[i][k];
                }
                delete [] fL2[i];
                delete [] fR2[i];
                delete [] inteBoud[i];
            }
            delete [] fL2;
            delete [] fR2;
            delete [] DLL2;
            delete [] DRL2;
            delete [] N_z;
            delete [] inteBoud;
        }//rod chain
        else if(MODEL[1] == "semiflexible" || MODEL[1] == "semi-flexible")//stiff chain
        {
            
            int*   nIntB;
            int*   delta_t;
            double*   f_0;
            double*   f_1;
            double*   f_2;
            double*** fL2;
            double*** fR2;
            double  nDR0,nDR1;
            int     nBessel,nBessel0,mp12,delta_Z1,delta_Z0,delta_Zt;
            int     kin1,kip1,midn,mIntB0,mIntB1,mIntBm,uppB,lowB;
            
            mp12    = mp22 - 1;
            nBessel = 2*round(1.0/dr) + 1;
            midn    = nBessel/2;
            nBessel0= midn*2;
            nIntB= new int[mp22]();
            f_0  = new double[nBessel]();
            f_1  = new double[nBessel]();
            f_2  = new double[nBessel]();
            fL2  = new double**[mp22]();
            fR2  = new double**[mp22]();
            delta_t = new int[MAXR1+1];
            for(short i=0; i<mp22; ++i)
            {
                nIntB[i] = round(DR2[i]/dr);
                fL2[i] = new double*[ngrid+1]();
                fR2[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL2[i][k] = new double[nBessel]();
                    fR2[i][k] = new double[nBessel]();
                }
            }
            
            mIntB0 = nIntB[0];
            mIntBm = nIntB[mp12];
            nDR0   = DR2[0];
            nDR1   = DR2[mp12];
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntB0;
                kip = k + mIntB0;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR0) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fL2[0][k][delta_Z0]    = 1.0;
                    fL2[0][ngrid-k][delta_Z1]    = 1.0;
                }
                
                
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR1) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fR2[mp12][k][delta_Z0] = 1.0;
                    fR2[mp12][ngrid-k][delta_Z1] = 1.0;
                }
            }
            
            for(short i1=1; i1<mp22; ++i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1-1];
                nDR0   = DR2[i1];
                nDR1   = DR2[i1-1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_1
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    kin1 = k - mIntB1;
                    kip1 = k + mIntB1;
                    if(kin1 < LLIM) kin1 = LLIM;
                    if(kip1 > ULIM) kip1 = ULIM;
                    
                    lowB = round((k-kip1)/nDR1) + midn;
                    uppB = round((k-kin1)/nDR1) + midn;
                    for(int k2=kin1; k2<=kip1; ++k2)
                    {//loop for Z_0
                        delta_Z0 = round((k-k2)/nDR1) + midn;
                        f_0[delta_Z0]  = ff2[i1-1][k2];
                        f_2[delta_Z0]  = fL2[i1-1][k2][delta_Z0];
                        delta_t[k2]    = delta_Z0;
                    }//loop for Z_0
                    
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z_2
                        //Nz2 = k + k1 - mIntB0;
                        delta_Z1 = round((k1 - k)/nDR0) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_0
                            delta_Z0 = delta_t[k2];
                            f_1[delta_Z0]  = f_0[delta_Z0]*BesselZero[1][delta_Z0][delta_Z1];
                        }//loop for Z_0
                        ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp2 = ftemp2*D_2;//without renormalization
                        ftemp2 = ftemp2/GRC2[delta_Z1];
                        fL2[i1][k][delta_Z1]= ftemp2;
                        
                        delta_Zt = nBessel0 - delta_Z1;
                        fL2[i1][ngrid-k][delta_Zt] = ftemp2;
                    }//loop for Z_2
                }//loop for Z_1
                
            }//loop for beads
            
            
            for(short i1=(mp12-1); i1>=0; --i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1+1];
                nDR0   = DR2[i1];
                nDR1   = DR2[i1+1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_0
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //Nz1 = k;
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z1
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
                        
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR2[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff2[i1+2][k2]*BesselZero[1][delta_Z0][delta_Z1];
                        }//loop for Z_2
                        
                        ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp2 = ftemp2*D_2;//without renormalization
                        ftemp2 = ftemp2/GRC2[delta_Z0];
                        fR2[i1][k][delta_Z0]= ftemp2;
                        
                        delta_Zt = nBessel0 - delta_Z0;
                        fR2[i1][ngrid-k][delta_Zt] = ftemp2;
                    }//loop for Z_1
                }//loop for Z_0
                
            }//loop for beads
            
            
            for(short i1=0; i1<mp22; ++i1)
            {
                mIntB0 = nIntB[i1];
                nDR0   = DR2[i1];
                D_2    = 0.5/DR2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    lowB = round((kin-k)/nDR0) + midn;
                    uppB = round((kip-k)/nDR0) + midn;
                    for(int k1=kin; k1<=kip; ++k1)
                    {
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        f_1[delta_Z0]  = ff2[i1+1][k1];
                        f_2[delta_Z0]  = fL2[i1][k][delta_Z0]*fR2[i1][k][delta_Z0];
                    }
                    
                    ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s2[i1][k] = ff2[i1][k]*ftemp2*D_2;
                }//loop for Z
            }
            
            ///////////////////////////i1 = mp22////////////////////////////
            mIntBm = nIntB[mp12];
            nDR1   = DR2[mp12];
            D_2    = 0.5/DR2[mp12];;
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                lowB = round((k-kip)/nDR1) + midn;
                uppB = round((k-kin)/nDR1) + midn;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k-k1)/nDR1) + midn;
                    f_1[delta_Z0]  = ff2[mp12][k1];
                    f_2[delta_Z0]  = fL2[mp12][k1][delta_Z0];
                }
                
                ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s2[mp22][k] = ff2[mp22][k]*ftemp2*D_2;
            }
            ///////////////////////////////////////////////////////////////
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    if(k<=ngrid_b)
                    {
                        rho1[i0][k] = 0.0;
                        for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                        {
                            rho1[i0][k]   = rho1[i0][k] + rho_s2[j2][k];
                        }
                    }
                    else
                    {
                        rho1[i0][k]   = rhoB_S[i0];
                    }
                    rho1[i0][ngrid-k] = rho1[i0][k];
                }
            }
                
                
                
            
            
            for(short i=0; i<mp22; ++i)
            {
                for(int k=0; k<=ngrid; ++k)
                {
                    delete [] fL2[i][k];
                    delete [] fR2[i][k];
                }
                delete [] fL2[i];
                delete [] fR2[i];
            }
            delete [] fL2;
            delete [] fR2;
            delete [] nIntB;
            delete [] f_0;
            delete [] f_1;
            delete [] f_2;
            delete [] delta_t;
            
            
        }//stiff chain
        
        else
        {
            cerr<<"Error in eulerlagrangeequation.cpp: Polymer model is wrong"<<endl;
            exit(0);
        }
    }
    //loop-3
    
    
    ///Here we only consider single surface. If two surface, we have to change this part/////
    //double temp01;
    //for(short i=0; i<hspecies; ++i)
    //{
    //    temp01 = rhoB_S[i]/rho1[i][ngrid_b];
    //    for(int k=LLIM; k<=ngrid_b; ++k)
    //    {
    //        rho1[i][k] = rho1[i][k]*temp01;
    //        rho1[i][ngrid-k] = rho1[i][k];
    //    }
    //}
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    if(neutsys == 1)
    {
       Renorm_Newton_Downhill(sigma,err,rhoBM1,rhoBM2,Z,rho1);
       //if(ngrid_b <ngrid_m) RenormDensity_Newton(sigma,deltaPhi,err,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,BesselZero,MODEL);
       //if(ngrid_b==ngrid_m) Renorm_Newton_Downhill(sigma,err,rhoBM1,rhoBM2,Z,rho1);
       //if(ngrid_b==ngrid_m) RenormDensity_Newton(sigma,err,rhoBM1,rhoBM2,Z,rho1);
    }
    
    VolumeFractionDFT(eta_t,D,LLI,ULI,etar,rho1);


    delete [] eta;
    delete [] temp;
    delete [] u_im;
    delete [] MB1;
    delete [] MB2;
    delete [] VI;
    
    if(rhoBM1 > 1E-10)
    {
        for(short i=0; i<mp1; ++i)
        {
            delete [] ff1[i];
            delete [] rho_s1[i];
        }
        delete [] ff1;
        delete [] rho_s1;
    }
    
    if(rhoBM2 > 1E-10)
    {
        for(short i=0; i<mp2; ++i)
        {
            delete [] ff2[i];
            delete [] rho_s2[i];
        }
        delete [] ff2;
        delete [] rho_s2;
    }
    
    
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] VanDW[i];
        delete [] DSH[i];
        delete [] DCH[i];
        delete [] DEC[i];
    }
    delete [] VanDW;
    delete [] DSH;
    delete [] DCH;
    delete [] DEC;
    
}




void EulerLagrangeDFT(double gama,int* LLI,int* ULI,float* D,double* BB,float* Z,double** pairEner,double* mu,
                      double* rhoB,double** ATT,double*** BesselZero,double*** Psi_IJ,string* MODEL)
{
    
    double  R,rin,rip;
    double  ftemp1,ftemp2; //,coe;
    double  rhoBM1,rhoBM2,D_2,Dav_1,Dav_2; //thresh
    double* VI;
    double* temp;
    double** rho;
    double** DSH;
    double** DCH;
    double** DEC;
    double** VanDW;
    double** ff1;
    double** ff2;
    double** rho_s1;
    double** rho_s2;
    short*   MB1;
    short*   MB2;
    
    
    short mp1,mp2,nblocks,i0,j0,mp11,mp22,hspecies; //ithresh
    int kin,kip,MAXR1,igammar;//,kip0
    
    
    struct InteLimit
    {
        int kin;    //lower limit
        int kip;    //upper limit
        
        InteLimit()
        {
            kin = 0;
            kip = 0;
        }
    };
    
    

    nblocks= nb[0] + nb[1];
    hspecies = nspecies - 1;
    
    MAXR1   = DMAX + ngrid_m;
    if(neutsys == 1)
    {
        igammar = round(1.0/(dr*gama));
        MAXR1   = DMAX + igammar + ngrid_m;
    }
    if(MAXR1 > ngrid) MAXR1 = ngrid;
    
    mp1 = mp[0];
    mp2 = mp[1];
    
    //thresh = 1.0E3;
    
    MB1  = new short[nb[0]+1]();
    MB2  = new short[nb[1]+1]();
    temp = new double[nblocks]();
    VI   = new double[hspecies]();
    
    rho  = new double*[nspecies]();
    for(short i=0; i<nspecies; ++i)
    {
        rho[i] = new double[ngrid+1]();
    }
    
    DCH  = new double*[nspecies]();
    DEC  = new double*[nspecies]();
    DSH  = new double*[nspecies]();
    VanDW= new double*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        VanDW[i]= new double[ngrid_m+1]();
        DCH[i]  = new double[ngrid_m+1]();
        DEC[i]  = new double[ngrid_m+1]();
        DSH[i]  = new double[ngrid_m+1]();
    }

    
    rhoBM1 = 0.0;
    MB1[0] = 0;
    Dav_1  = 0;
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1   = rhoBM1 + rhoB[i];
        MB1[i+1] = MB1[i] + mb[0][i];
        Dav_1   += D[i];
    }
    rhoBM2 = 0.0;
    MB2[0] = 0;
    Dav_2  = 0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        j0        = i - nb[0];
        rhoBM2    = rhoBM2 + rhoB[i];
        MB2[j0+1] = MB2[j0] + mb[1][j0];
        Dav_2    += D[i];
    }
    
    if(rhoBM1 > 1E-10)
    {
        ff1  = new double*[mp1]();
        rho_s1 = new double*[mp1]();
        for(short i=0; i<mp1; ++i)
        {
            ff1[i] = new double[ngrid+1]();
            rho_s1[i] = new double[ngrid+1]();
        }
    }
    
    if(rhoBM2 > 1E-10)
    {
        ff2  = new double*[mp2]();
        rho_s2 = new double*[mp2]();
        for(short i=0; i<mp2; ++i)
        {
            ff2[i] = new double[ngrid+1]();
            rho_s2[i] = new double[ngrid+1]();
        }
    }
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        for(short i=0; i<nspecies; ++i)
        {
            rho[i][k] = rhoB[i];
            if((k<LLI[i]) || (k>ULI[i])) rho[i][k] = 0;
            
            rho[i][ngrid-k] = rho[i][k];
        }
    }
    
    if(neutsys == 1)
    {
        //electrostatic correlation from MSA
        DerivElectroCorrel(gama,LLI,ULI,D,Z,rhoB,ATT,rho,DEC);
        //DerivElectroCorrelDirk(MAXR1,gama,LLI,ULI,D,Z,rho,DEC);
        
        //electrostatic contributions from charge shell
        ChargeShell(LLI,ULI,BB,rho,Psi_IJ,DSH);
        //ChargeShellDirk(LLI,ULI,D,rho,Psi_IJ,DSH);
    }
    
    //Square-well potential
    InhomVanDelWaal(LLI,ULI,D,rho,pairEner,VanDW);
    
    //hard sphere + hard sphere chain
    DerivHardSphereChain(MAXR1,LLI,ULI,rhoB,D,Z,ATT,rho,DCH);
    
    
    for(short j=nblocks; j<hspecies; ++j)
    {
        rhoB_S[j] = exp(mu[j]-VanDW[j][ngrid_m]-DEC[j][ngrid_m]-DSH[j][ngrid_m]
                        -DCH[j][ngrid_m]);
    }
    //loop-1
    for(int k=LLIM; k<=ngrid_m; ++k)
    {
        if(rhoBM1 > 1E-10)
        {
            for(short j=0; j<nb[0]; ++j)
            {
                temp[j] = 0;
                if(k >= LLI[j])
                {
                    temp[j]  = exp(mu[j]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                   -DCH[j][k]);
                }
                if(k==ngrid_m)  temp_S[j]= temp[j];
            }
            
            
            for(short i2=0; i2<nb[0]; ++i2)
            {
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    ff1[j2][k] = temp[i2];
                    ff1[j2][ngrid-k] = ff1[j2][k];
                }
            }
            
        }
        
        
        if(rhoBM2 > 1E-10)
        {
            for(short j=nb[0]; j<nblocks; ++j)
            {
                temp[j] = 0;
                if(k >= LLI[j])
                {
                    temp[j]  = exp(mu[j]-VanDW[j][k]-DEC[j][k]-DSH[j][k]
                                   -DCH[j][k]);
                }
                if(k==ngrid_m)  temp_S[j]= temp[j];
            }
            
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    ff2[j2][k] = temp[i0];
                    ff2[j2][ngrid-k] = ff2[j2][k];
                }
            }
        }
        
    }
    //loop-1
    
    
    //ithresh = 0;
    //loop-2
    
    if(rhoBM1 > 1E-10)
    {
        mp11 = mp1 - 1;
        if(MODEL[0] == "flexible")//flexible chain
        {
            
            double** fL1;
            double** fR1;
            fL1  = new double*[mp1]();
            fR1  = new double*[mp1]();
            for(short i=0; i<mp1; ++i)
            {
                fL1[i] = new double[ngrid+1]();
                fR1[i] = new double[ngrid+1]();
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                fL1[0][k]   = 1.0;
                fR1[mp11][k]= 1.0;
                
                fL1[0][ngrid-k]   = fL1[0][k];
                fR1[mp11][ngrid-k]= fR1[mp11][k];
            }
            
            for(short i1=1; i1<mp1; ++i1)
            {
                D_2 = 0.5/DL1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DL1[i1];
                    rip = R + DL1[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp1=RombergIntegration(ff1[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    ftemp1=SimpsonIntegration(ff1[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL1[i1][k] = ftemp1*D_2;
                    fL1[i1][ngrid-k] = fL1[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            
            for(short i1=(mp1-2); i1>=0; --i1)
            {
                D_2 = 0.5/DR1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DR1[i1];
                    rip = R + DR1[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp1=RombergIntegration(ff1[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    ftemp1=SimpsonIntegration(ff1[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR1[i1][k] = ftemp1*D_2;
                    fR1[i1][ngrid-k] = fR1[i1][k];
                }
            }
            
            for(short i2=0; i2<nb[0]; ++i2)
            {
                rhoB_S[i2] = 0;
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    rho_s1[j2][ngrid_m] =  ff1[j2][ngrid_m]*fL1[j2][ngrid_m]*fR1[j2][ngrid_m];
                    rhoB_S[i2]   =   rhoB_S[i2] + rho_s1[j2][ngrid_m];
                }
            }
            
            
            for(short i=0; i<mp1; ++i)
            {
                delete [] fL1[i];
                delete [] fR1[i];
            }
            delete [] fL1;
            delete [] fR1;
        }//flexible chain
        else if(MODEL[0] == "rod")//rod chain
        {
            
            int* N_z;
            double* DLL1;
            double* DRL1;
            double*** fL1;
            double*** fR1;
            double  slope_k,temp0;
            int     Nz2,Nz1,Nz0,Nzn,k2,k3,orient1,orient,ngNz0;
            
            InteLimit** inteBoud;
            
            Dav_1   = Dav_1/nb[0];
            orient1 = round(Dav_1/dr);
            orient  = 2*orient1;
            N_z  = new int[mp1]();
            DLL1 = new double[mp1]();
            DRL1 = new double[mp1]();
            fL1  = new double**[mp1]();
            fR1  = new double**[mp1]();
            
            inteBoud = new InteLimit*[mp1];
            for(short i=0; i<mp1; ++i)
            {
                fL1[i] = new double*[ngrid+1]();
                fR1[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL1[i][k] = new double[orient+1]();
                    fR1[i][k] = new double[orient+1]();
                }
                
                inteBoud[i] = new InteLimit[ngrid+1];
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(int k1=0; k1<=orient; ++k1)
                {
                    fL1[0][k][k1]    = 1.0;
                    fR1[mp11][k][k1] = 1.0;
                    
                    fL1[0][ngrid-k][k1]    = 1.0;
                    fR1[mp11][ngrid-k][k1] = 1.0;
                }
            }
            
            DLL1[0]   = 0;
            DRL1[mp11]= 0;
            for(short i1=1; i1<mp1; ++i1)
            {
                DLL1[i1]      =  DLL1[i1-1] + DL1[i1];
                DRL1[mp11-i1] =  DRL1[mp11-i1+1] + DR1[mp11-i1];
            }
            
            //orientation points to (the last bead - the first bead)
            for(int k=LLIM; k<=ngrid_m; ++k)
            {//loop for grids
                for(int k1=(-orient1); k1<=orient1; ++k1)
                {//loop for orientation
                    k2 = k1 + orient1;
                    k3 = orient1 - k1;
                    
                    slope_k = k1/Dav_1;
                    for(short i1=0; i1<mp1; ++i1)
                    {
                        N_z[i1] = k + round(slope_k*DLL1[i1]);
                        if(N_z[i1] < 0) N_z[i1] = 0;
                        if(N_z[i1] > ngrid) N_z[i1] = ngrid;
                    }
                    //slope_kR = -slope_kL;
                    //Nzn = k - round(slope_kR*DRL1[0]);
                    
                    Nz0  = N_z[0];
                    Nzn  = N_z[mp11];
                    ngNz0= ngrid - Nz0;
                    if((Nzn >= LLIM) && (Nzn <= ULIM))
                    {
                        temp0 = 1;
                        for(short i2=1; i2<mp1; ++i2)
                        {
                            Nz1    = N_z[i2-1];
                            Nz2    = N_z[i2];
                            
                            temp0 *= ff1[i2-1][Nz1];
                            fL1[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL1[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff1[i2+1][Nz2];
                            fR1[i2][Nz1][k2]       = temp0;
                            if(Nz0 != ngNz0) fR1[i2][ngrid-Nz1][k3] = temp0;
                        }
                        
                        if(k1 < 0)
                        {
                            for(short i2=0; i2<mp1; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kin += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kip += 1;
                            }
                        }
                        else if(k1 > 0)
                        {
                            for(short i2=0; i2<mp1; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kip += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kin += 1;
                            }
                        }
                    }
                    else
                    {
                        for(short i2=1; i2<mp1; ++i2)
                        {
                            Nz1 = N_z[i2];
                            fL1[i2][Nz1][k2]       = 0;
                            fL1[i2][ngrid-Nz1][k3] = 0;
                        }
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1 = N_z[i2];
                            fR1[i2][Nz1][k2]       = 0;
                            fR1[i2][ngrid-Nz1][k3] = 0;
                        }
                        
                    }
                }//loop for orientation
                
            }//loop for grids
            
            
            D_2 = 0.5/Dav_1;
            for(short i1=0; i1<mp1; ++i1)
            {//loop for beads
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for grids
                    kin = orient1 - inteBoud[i1][k].kin;
                    kip = orient1 + inteBoud[i1][k].kip;
                    
                    //ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],kin,kip,kin,kip);
                    ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],orient1,kin,kip);
                    rho_s1[i1][k] = ff1[i1][k]*ftemp1*D_2;
                }//loop for grids
            }//loop for beads
            
            
            
            for(short i2=0; i2<nb[0]; ++i2)
            {
                rhoB_S[i2] = 0;
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    rhoB_S[i2]   =  rhoB_S[i2] + rho_s1[j2][ngrid_m];
                }
                
            }
            
            
            for(short i=0; i<mp1; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL1[i][k];
                    delete [] fR1[i][k];
                }
                delete [] fL1[i];
                delete [] fR1[i];
                delete [] inteBoud[i];
            }
            delete [] fL1;
            delete [] fR1;
            delete [] DLL1;
            delete [] DRL1;
            delete [] N_z;
            delete [] inteBoud;
        }//rod chain
        else if(MODEL[0] == "semiflexible" || MODEL[0] == "semi-flexible")//stiff chain
        {
            
            int*   nIntB;
            int*   delta_t;
            double*   f_0;
            double*   f_1;
            double*   f_2;
            double*** fL1;
            double*** fR1;
            double  nDR0,nDR1;
            int     nBessel,nBessel0,mp12,delta_Z1,delta_Z0,delta_Zt;
            int     kin1,kip1,midn,mIntB0,mIntB1,mIntBm,uppB,lowB;
            
            mp12    = mp11 - 1;
            nBessel = 2*round(1.0/dr) + 1;
            midn    = nBessel/2;
            nBessel0= midn*2;
            nIntB= new int[mp11]();
            f_0  = new double[nBessel]();
            f_1  = new double[nBessel]();
            f_2  = new double[nBessel]();
            fL1  = new double**[mp11]();
            fR1  = new double**[mp11]();
            delta_t = new int[MAXR1+1];
            for(short i=0; i<mp11; ++i)
            {
                nIntB[i] = round(DR1[i]/dr);
                fL1[i] = new double*[ngrid+1]();
                fR1[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL1[i][k] = new double[nBessel]();
                    fR1[i][k] = new double[nBessel]();
                }
            }
            
            
            mIntB0 = nIntB[0];
            mIntBm = nIntB[mp12];
            nDR0   = DR1[0];
            nDR1   = DR1[mp12];
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntB0;
                kip = k + mIntB0;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR0) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fL1[0][k][delta_Z0]    = 1.0;
                    fL1[0][ngrid-k][delta_Z1]    = 1.0;
                }
                
                
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR1) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fR1[mp12][k][delta_Z0] = 1.0;
                    fR1[mp12][ngrid-k][delta_Z1] = 1.0;
                }
            }
            
            for(short i1=1; i1<mp11; ++i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1-1];
                nDR0   = DR1[i1];
                nDR1   = DR1[i1-1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_1
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    kin1 = k - mIntB1;
                    kip1 = k + mIntB1;
                    if(kin1 < LLIM) kin1 = LLIM;
                    if(kip1 > ULIM) kip1 = ULIM;
                    lowB = round((k-kip1)/nDR1) + midn;
                    uppB = round((k-kin1)/nDR1) + midn;
                    for(int k2=kin1; k2<=kip1; ++k2)
                    {//loop for Z_0
                        delta_Z0 = round((k-k2)/nDR1) + midn;
                        f_0[delta_Z0]  = ff1[i1-1][k2];
                        f_2[delta_Z0]  = fL1[i1-1][k2][delta_Z0];
                        delta_t[k2]    = delta_Z0;
                    }//loop for Z_0
                    
                    
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z_2
                        //Nz2 = k + k1 - mIntB0;
                        delta_Z1 = round((k1 - k)/nDR0) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_0
                            delta_Z0 = delta_t[k2];
                            f_1[delta_Z0]  = f_0[delta_Z0]*BesselZero[0][delta_Z0][delta_Z1];
                        }//loop for Z_0
                        
                        
                        ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp1 = ftemp1*D_2; //without renormalization
                        ftemp1 = ftemp1/GRC1[delta_Z1];
                        fL1[i1][k][delta_Z1]= ftemp1;
                        
                        delta_Zt = nBessel0 - delta_Z1;
                        fL1[i1][ngrid-k][delta_Zt] = ftemp1;
                        
                    }//loop for Z_2
                    
                }//loop for Z_1
                
            }//loop for beads
            
            for(short i1=(mp12-1); i1>=0; --i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1+1];
                nDR0   = DR1[i1];
                nDR1   = DR1[i1+1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_0
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //Nz1 = k;
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z1
                        delta_Z0 = round((k1 - k)/nDR0) + midn; //correct
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
                        
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            //delta_Z1 = delta_1[k2];
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR1[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff1[i1+2][k2]*BesselZero[0][delta_Z0][delta_Z1];
                        }//loop for Z_2
                        
                        ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp1 = ftemp1*D_2; //without renormalization
                        ftemp1 = ftemp1/GRC1[delta_Z0];
                        fR1[i1][k][delta_Z0]= ftemp1;
                        
                        delta_Zt = nBessel0 - delta_Z0;
                        fR1[i1][ngrid-k][delta_Zt] = ftemp1;
                    }//loop for Z_1
                }//loop for Z_0
                
            }//loop for beads
            
            for(short i1=0; i1<mp11; ++i1)
            {
                mIntB0 = nIntB[i1];
                nDR0   = DR1[i1];
                D_2    = 0.5/DR1[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    lowB = round((kin-k)/nDR0) + midn;
                    uppB = round((kip-k)/nDR0) + midn;
                    for(int k1=kin; k1<=kip; ++k1)
                    {
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        f_1[delta_Z0]  = ff1[i1+1][k1];
                        f_2[delta_Z0]  = fL1[i1][k][delta_Z0]*fR1[i1][k][delta_Z0];
                    }
                    
                    ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s1[i1][k] = ff1[i1][k]*ftemp1*D_2;
                }//loop for Z
            }
            
            ///////////////////////////i1 = mp11////////////////////////////
            mIntBm = nIntB[mp12];
            nDR1   = DR1[mp12];
            D_2    = 0.5/DR1[mp12];;
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                lowB = round((k-kip)/nDR1) + midn;
                uppB = round((k-kin)/nDR1) + midn;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k-k1)/nDR1) + midn;
                    f_1[delta_Z0]  = ff1[mp12][k1];
                    f_2[delta_Z0]  = fL1[mp12][k1][delta_Z0];
                }
                
                ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s1[mp11][k] = ff1[mp11][k]*ftemp1*D_2;
            }
            
            
            for(short i2=0; i2<nb[0]; ++i2)
            {
                rhoB_S[i2] = 0;
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    rhoB_S[i2]  =  rhoB_S[i2] + rho_s1[j2][ngrid_m];
                }
            }
            
            
            
            for(short i=0; i<mp11; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL1[i][k];
                    delete [] fR1[i][k];
                }
                delete [] fL1[i];
                delete [] fR1[i];
            }
            delete [] fL1;
            delete [] fR1;
            delete [] nIntB;
            delete [] f_0;
            delete [] f_1;
            delete [] f_2;
            delete [] delta_t;
            
            
        }//stiff chain
        else
        {
            cerr<<"Error in eulerlagrangeequation.cpp: Polymer model is wrong"<<endl;
            exit(0);
        }
        
        
        
    }
    //loop-2
    
    
    //ithresh = 0;
    //loop-3
    if(rhoBM2 > 1E-10)
    {
        mp22 = mp2 - 1;
        if(MODEL[1] == "flexible")//flexible chain
        {
            double** fL2;
            double** fR2;
            fL2  = new double*[mp2]();
            fR2  = new double*[mp2]();
            for(short i=0; i<mp2; ++i)
            {
                fL2[i] = new double[ngrid+1]();
                fR2[i] = new double[ngrid+1]();
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                fL2[0][k]   = 1;
                fR2[mp22][k]= 1;
                
                fL2[0][ngrid-k]   = fL2[0][k];
                fR2[mp22][ngrid-k]= fR2[mp22][k];
            }
            
            for(short i1=1; i1<mp2; ++i1)
            {
                D_2 = 0.5/DL2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DL2[i1];
                    rip = R + DL2[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp2=RombergIntegration(ff2[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    ftemp2=SimpsonIntegration(ff2[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL2[i1][k] = ftemp2*D_2;
                    fL2[i1][ngrid-k] = fL2[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(short i1=(mp2-2); i1>=0; --i1)
            {
                D_2 = 0.5/DR2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {
                    R   = k*dr;
                    rin = R - DR2[i1];
                    rip = R + DR2[i1];
                    kin = round(rin/dr);
                    kip = round(rip/dr);
                    
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //ftemp2=RombergIntegration(ff2[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    ftemp2=SimpsonIntegration(ff2[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR2[i1][k] = ftemp2*D_2;
                    fR2[i1][ngrid-k] = fR2[i1][k];
                }
                //if(ithresh == 1) break;
            }
            
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                rhoB_S[i0] = 0.0;
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    rho_s2[j2][ngrid_m] = ff2[j2][ngrid_m]*fL2[j2][ngrid_m]*fR2[j2][ngrid_m];
                    rhoB_S[i0]   = rhoB_S[i0] + rho_s2[j2][ngrid_m];
                }
            }
            
            
            for(short i=0; i<mp2; ++i)
            {
                delete [] fL2[i];
                delete [] fR2[i];
            }
            delete [] fL2;
            delete [] fR2;
        }//flexible chain
        else if(MODEL[1] == "rod")//rod chain
        {
            int* N_z;
            double* DLL2;
            double* DRL2;
            double*** fL2;
            double*** fR2;
            double  slope_k,temp0;
            int     Nz2,Nz1,Nz0,Nzn,k2,k3,orient2,orient,ngNz0;
            
            InteLimit** inteBoud;
            
            Dav_2   = Dav_2/nb[1];
            orient2 = round(Dav_2/dr);
            orient  = 2*orient2;
            N_z  = new int[mp2]();
            DLL2 = new double[mp2]();
            DRL2 = new double[mp2]();
            fL2  = new double**[mp2]();
            fR2  = new double**[mp2]();
            
            inteBoud = new InteLimit*[mp2];
            for(short i=0; i<mp2; ++i)
            {
                fL2[i] = new double*[ngrid+1]();
                fR2[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL2[i][k] = new double[orient+1]();
                    fR2[i][k] = new double[orient+1]();
                }
                
                inteBoud[i] = new InteLimit[ngrid+1];
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(int k1=0; k1<=orient; ++k1)
                {
                    fL2[0][k][k1]    = 1.0;
                    fR2[mp22][k][k1] = 1.0;
                    
                    fL2[0][ngrid-k][k1]    = 1.0;
                    fR2[mp22][ngrid-k][k1] = 1.0;
                }
            }
            
            DLL2[0]   = 0;
            DRL2[mp22]= 0;
            for(short i1=1; i1<mp2; ++i1)
            {
                DLL2[i1]      =  DLL2[i1-1] + DL2[i1];
                DRL2[mp22-i1] =  DRL2[mp22-i1+1] + DR2[mp22-i1];
            }
            
            //orientation points to (the last bead - the first bead)
            for(int k=LLIM; k<=ngrid_m; ++k)
            {//loop for grids
                for(int k1=(-orient2); k1<=orient2; ++k1)
                {//loop for orientation
                    k2 = k1 + orient2;
                    k3 = orient2 - k1;
                    
                    slope_k = k1/Dav_2;
                    for(short i1=0; i1<mp2; ++i1)
                    {
                        N_z[i1] = k + round(slope_k*DLL2[i1]);
                        if(N_z[i1] < 0) N_z[i1] = 0;
                        if(N_z[i1] > ngrid) N_z[i1] = ngrid;
                    }
                    //slope_kR = -slope_kL;
                    //Nzn = k - round(slope_kR*DRL1[0]);
                    
                    Nz0  = N_z[0];
                    Nzn  = N_z[mp22];
                    ngNz0= ngrid - Nz0;
                    if((Nzn >= LLIM) && (Nzn <= ULIM))
                    {
                        temp0 = 1;
                        for(short i2=1; i2<mp2; ++i2)
                        {
                            Nz1    = N_z[i2-1];
                            Nz2    = N_z[i2];
                            
                            temp0 *= ff2[i2-1][Nz1];
                            fL2[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL2[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff2[i2+1][Nz2];
                            fR2[i2][Nz1][k2]       = temp0;
                            if(Nz0 != ngNz0) fR2[i2][ngrid-Nz1][k3] = temp0;
                        }
                        
                        if(k1 < 0)
                        {
                            for(short i2=0; i2<mp2; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kin += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kip += 1;
                            }
                        }
                        else if(k1 > 0)
                        {
                            for(short i2=0; i2<mp2; ++i2)
                            {
                                Nz1 = N_z[i2];
                                inteBoud[i2][Nz1].kip += 1;
                                if(Nz0 != ngNz0) inteBoud[i2][ngrid - Nz1].kin += 1;
                            }
                        }
                    }
                    else
                    {
                        for(short i2=1; i2<mp2; ++i2)
                        {
                            Nz1 = N_z[i2];
                            fL2[i2][Nz1][k2]       = 0;
                            fL2[i2][ngrid-Nz1][k3] = 0;
                        }
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1 = N_z[i2];
                            fR2[i2][Nz1][k2]       = 0;
                            fR2[i2][ngrid-Nz1][k3] = 0;
                        }
                        
                    }
                }//loop for orientation
                
            }//loop for grids
            
            
            D_2 = 0.5/Dav_2;
            for(short i1=0; i1<mp2; ++i1)
            {//loop for beads
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for grids
                    kin = orient2 - inteBoud[i1][k].kin;
                    kip = orient2 + inteBoud[i1][k].kip;
                    
                    ftemp2=SimpsonIntegration(fL2[i1][k],fR2[i1][k],orient2,kin,kip);
                    rho_s2[i1][k] = ff2[i1][k]*ftemp2*D_2;
                }//loop for grids
            }//loop for beads
            
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                rhoB_S[i0] = 0.0;
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    rhoB_S[i0]   = rhoB_S[i0] + rho_s2[j2][ngrid_m];
                }
            }
            
            for(short i=0; i<mp2; ++i)
            {
                for(int k=0; k<(ngrid+1); ++k)
                {
                    delete [] fL2[i][k];
                    delete [] fR2[i][k];
                }
                delete [] fL2[i];
                delete [] fR2[i];
                delete [] inteBoud[i];
            }
            delete [] fL2;
            delete [] fR2;
            delete [] DLL2;
            delete [] DRL2;
            delete [] N_z;
            delete [] inteBoud;
        }//rod chain
        else if(MODEL[1] == "semiflexible" || MODEL[1] == "semi-flexible")//stiff chain
        {
            
            int*   nIntB;
            int*   delta_t;
            double*   f_0;
            double*   f_1;
            double*   f_2;
            double*** fL2;
            double*** fR2;
            double  nDR0,nDR1;
            int     nBessel,nBessel0,mp12,delta_Z1,delta_Z0,delta_Zt;
            int     kin1,kip1,midn,mIntB0,mIntB1,mIntBm,uppB,lowB;
            
            mp12    = mp22 - 1;
            nBessel = 2*round(1.0/dr) + 1;
            midn    = nBessel/2;
            nBessel0= midn*2;
            nIntB= new int[mp22]();
            f_0  = new double[nBessel]();
            f_1  = new double[nBessel]();
            f_2  = new double[nBessel]();
            fL2  = new double**[mp22]();
            fR2  = new double**[mp22]();
            delta_t = new int[MAXR1+1];
            for(short i=0; i<mp22; ++i)
            {
                nIntB[i] = round(DR2[i]/dr);
                fL2[i] = new double*[ngrid+1]();
                fR2[i] = new double*[ngrid+1]();
                for(int k=0; k<(ngrid+1); ++k)
                {
                    fL2[i][k] = new double[nBessel]();
                    fR2[i][k] = new double[nBessel]();
                }
            }
            
            mIntB0 = nIntB[0];
            mIntBm = nIntB[mp12];
            nDR0   = DR2[0];
            nDR1   = DR2[mp12];
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntB0;
                kip = k + mIntB0;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR0) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fL2[0][k][delta_Z0]    = 1.0;
                    fL2[0][ngrid-k][delta_Z1]    = 1.0;
                }
                
                
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k1-k)/nDR1) + midn;
                    delta_Z1 = nBessel0 - delta_Z0;
                    fR2[mp12][k][delta_Z0] = 1.0;
                    fR2[mp12][ngrid-k][delta_Z1] = 1.0;
                }
            }
            
            for(short i1=1; i1<mp22; ++i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1-1];
                nDR0   = DR2[i1];
                nDR1   = DR2[i1-1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_1
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    kin1 = k - mIntB1;
                    kip1 = k + mIntB1;
                    if(kin1 < LLIM) kin1 = LLIM;
                    if(kip1 > ULIM) kip1 = ULIM;
                    
                    lowB = round((k-kip1)/nDR1) + midn;
                    uppB = round((k-kin1)/nDR1) + midn;
                    for(int k2=kin1; k2<=kip1; ++k2)
                    {//loop for Z_0
                        delta_Z0 = round((k-k2)/nDR1) + midn;
                        f_0[delta_Z0]  = ff2[i1-1][k2];
                        f_2[delta_Z0]  = fL2[i1-1][k2][delta_Z0];
                        delta_t[k2]    = delta_Z0;
                    }//loop for Z_0
                    
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z_2
                        //Nz2 = k + k1 - mIntB0;
                        delta_Z1 = round((k1 - k)/nDR0) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_0
                            delta_Z0 = delta_t[k2];
                            f_1[delta_Z0]  = f_0[delta_Z0]*BesselZero[1][delta_Z0][delta_Z1];
                        }//loop for Z_0
                        ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp2 = ftemp2*D_2;//without renormalization
                        ftemp2 = ftemp2/GRC2[delta_Z1];
                        fL2[i1][k][delta_Z1]= ftemp2;
                        
                        delta_Zt = nBessel0 - delta_Z1;
                        fL2[i1][ngrid-k][delta_Zt] = ftemp2;
                    }//loop for Z_2
                }//loop for Z_1
                
            }//loop for beads
            
            
            for(short i1=(mp12-1); i1>=0; --i1)
            {//loop for beads
                mIntB0 = nIntB[i1];
                mIntB1 = nIntB[i1+1];
                nDR0   = DR2[i1];
                nDR1   = DR2[i1+1];
                //D_2    = 0.5/nDR1;
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z_0
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    
                    //Nz1 = k;
                    for(int k1=kin; k1<=kip; ++k1)
                    {//loop for Z1
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
                        
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR2[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff2[i1+2][k2]*BesselZero[1][delta_Z0][delta_Z1];
                        }//loop for Z_2
                        
                        ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                        //ftemp2 = ftemp2*D_2;//without renormalization
                        ftemp2 = ftemp2/GRC2[delta_Z0];
                        fR2[i1][k][delta_Z0]= ftemp2;
                        
                        delta_Zt = nBessel0 - delta_Z0;
                        fR2[i1][ngrid-k][delta_Zt] = ftemp2;
                    }//loop for Z_1
                }//loop for Z_0
                
            }//loop for beads
            
            
            for(short i1=0; i1<mp22; ++i1)
            {
                mIntB0 = nIntB[i1];
                nDR0   = DR2[i1];
                D_2    = 0.5/DR2[i1];
                for(int k=LLIM; k<=ngrid_m; ++k)
                {//loop for Z
                    kin = k - mIntB0;
                    kip = k + mIntB0;
                    if(kin < LLIM) kin = LLIM;
                    if(kip > ULIM) kip = ULIM;
                    lowB = round((kin-k)/nDR0) + midn;
                    uppB = round((kip-k)/nDR0) + midn;
                    for(int k1=kin; k1<=kip; ++k1)
                    {
                        delta_Z0 = round((k1 - k)/nDR0) + midn;
                        f_1[delta_Z0]  = ff2[i1+1][k1];
                        f_2[delta_Z0]  = fL2[i1][k][delta_Z0]*fR2[i1][k][delta_Z0];
                    }
                    
                    ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s2[i1][k] = ff2[i1][k]*ftemp2*D_2;
                }//loop for Z
            }
            
            ///////////////////////////i1 = mp22////////////////////////////
            mIntBm = nIntB[mp12];
            nDR1   = DR2[mp12];
            D_2    = 0.5/DR2[mp12];;
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                kin = k - mIntBm;
                kip = k + mIntBm;
                if(kin < LLIM) kin = LLIM;
                if(kip > ULIM) kip = ULIM;
                lowB = round((k-kip)/nDR1) + midn;
                uppB = round((k-kin)/nDR1) + midn;
                for(int k1=kin; k1<=kip; ++k1)
                {
                    delta_Z0 = round((k-k1)/nDR1) + midn;
                    f_1[delta_Z0]  = ff2[mp12][k1];
                    f_2[delta_Z0]  = fL2[mp12][k1][delta_Z0];
                }
                
                ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s2[mp22][k] = ff2[mp22][k]*ftemp2*D_2;
            }
            ///////////////////////////////////////////////////////////////
            
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                rhoB_S[i0] = 0.0;
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    rhoB_S[i0]  = rhoB_S[i0] + rho_s2[j2][ngrid_m];
                }
            }
            
            
            
            
            
            for(short i=0; i<mp22; ++i)
            {
                for(int k=0; k<=ngrid; ++k)
                {
                    delete [] fL2[i][k];
                    delete [] fR2[i][k];
                }
                delete [] fL2[i];
                delete [] fR2[i];
            }
            delete [] fL2;
            delete [] fR2;
            delete [] nIntB;
            delete [] f_0;
            delete [] f_1;
            delete [] f_2;
            delete [] delta_t;
            
            
        }//stiff chain
        
        else
        {
            cerr<<"Error in eulerlagrangeequation.cpp: Polymer model is wrong"<<endl;
            exit(0);
        }
    }
    //loop-3
    
    delete [] temp;
    delete [] MB1;
    delete [] MB2;
    delete [] VI;
    
    if(rhoBM1 > 1E-10)
    {
        for(short i=0; i<mp1; ++i)
        {
            delete [] ff1[i];
            delete [] rho_s1[i];
        }
        delete [] ff1;
        delete [] rho_s1;
    }
    
    if(rhoBM2 > 1E-10)
    {
        for(short i=0; i<mp2; ++i)
        {
            delete [] ff2[i];
            delete [] rho_s2[i];
        }
        delete [] ff2;
        delete [] rho_s2;
    }
    
    for(short i=0; i<nspecies; ++i)
    {
        delete [] rho[i];
    }
    for(short i=0; i<hspecies; ++i)
    {
        delete [] VanDW[i];
    }
    delete [] rho;
    delete [] VanDW;
    
}
