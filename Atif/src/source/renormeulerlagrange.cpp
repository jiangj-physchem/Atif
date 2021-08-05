//*******************solve Euler Lagrange***************//
#include "renormeulerlagrange.h"
//#include "rombergintegration.h"
#include "simpsonintegration.h"
#include "constantnum.h"

extern double dr;
extern int ngrid; //the number of grids
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int ngrid_b;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern int ULIM; //the maxium uper limit of intergral
extern int LLIM; //the minimum lower limit of intergral
extern int DMAX;
extern short nspecies;
extern double* DL1;
extern double* DR1;
extern double* DL2;
extern double* DR2;
extern double* GRC1;
extern double* GRC2;

double RenormEulerLagrange(double ratio,double rhoBM1,double rhoBM2,float* D,short* MB1,short* MB2,float* Z,
                   double** ff1,double** ff2,double** rho1,double** rho2,double*** BesselZero,string* MODEL)
{
    double  R,rin,rip,temp,D_2,Dav_1,Dav_2;
    double  ftemp1,ftemp2,fsigma,size_b,dr_size_b,R_L; //,coe
    double* Psi_delta;
    
    double** rho_s1;
    double** rho_s2;
    double** ff11;
    double** ff22;
    
    short mp1,mp2,nblocks,i0,mp11,mp22,hspecies;
    int kin,kip,MAXR0,MAXR1;//kip0;
    
    hspecies= nspecies - 1;
    
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
    
    //coe     = dr*0.5;
    nblocks = nb[0] + nb[1];
    
    MAXR0   = DMAX/2 + ngrid_m;
    MAXR1   = DMAX + ngrid_m;
    
    size_b    = ngrid_b*dr;
    //if(ngrid%2 != 0) size_mid = ngrid_m*dr + 0.5*dr;
    dr_size_b = dr/size_b;
    
    mp1 = mp[0];
    mp2 = mp[1];
    
    Psi_delta = new double[ngrid_m+1]();
    
    Dav_1  = 0;
    for(short i=0; i<nb[0]; ++i)
    {
        Dav_1 += D[i];
    }
    Dav_2  = 0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        Dav_2 += D[i];
    }
    
    
    if(rhoBM1 > 1E-10)
    {
        ff11 = new double*[mp1]();
        rho_s1 = new double*[mp1]();
        for(short i=0; i<mp1; ++i)
        {
            ff11[i]= new double[ngrid+1]();
            rho_s1[i] = new double[ngrid+1]();
        }
    }
    
    if(rhoBM2 > 1E-10)
    {
        ff22 = new double*[mp2]();
        rho_s2 = new double*[mp2]();
        for(short i=0; i<mp2; ++i)
        {
            ff22[i]= new double[ngrid+1]();
            rho_s2[i] = new double[ngrid+1]();
        }
    }
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        if(k<=ngrid_b)
        {
            //R_L = (k*dr)/size_b;
            R_L = k*dr_size_b;
            Psi_delta[k] =  ratio*R_L - ratio;
        }
        else
        {
            Psi_delta[k] = 0;
        }
        
    }
    
    //loop-1
    for(int k=LLIM; k<=ngrid_b; ++k)
    {
        for(short j=nblocks; j<hspecies; ++j)
        {
            rho2[j][k] = rho1[j][k]*exp(Psi_delta[k]*Z[j]);
        }
        
        
        
        if(rhoBM1 > 1E-10)
        {
            for(short i2=0; i2<nb[0]; ++i2)
            {
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    ff11[j2][k] =  ff1[j2][k]*exp(Psi_delta[k]*Z[i2]);
                    ff11[j2][ngrid-k] = ff11[j2][k];
                }
            }
        }
        
        
        if(rhoBM2 > 1E-10)
        {
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    ff22[j2][k] = ff2[j2][k]*exp(Psi_delta[k]*Z[i0]);
                    ff22[j2][ngrid-k] = ff22[j2][k];
                }
            }
        }
        
    }
    //loop-1
    
    
    //loop-1
    for(int k=(ngrid_b+1); k<=ngrid_m; ++k)
    {
        for(short j=nblocks; j<hspecies; ++j)
        {
            rho2[j][k] = rho1[j][k];
        }
        
        
        
        if(rhoBM1 > 1E-10)
        {
            for(short i2=0; i2<nb[0]; ++i2)
            {
                for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                {
                    ff11[j2][k] =  ff1[j2][k];
                    ff11[j2][ngrid-k] = ff11[j2][k];
                }
            }
        }
        
        
        if(rhoBM2 > 1E-10)
        {
            for(short i2=0; i2<nb[1]; ++i2)
            {
                i0 = i2 + nb[0];
                for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                {
                    ff22[j2][k] = ff2[j2][k];
                    ff22[j2][ngrid-k] = ff22[j2][k];
                }
            }
        }
        
    }
    //loop-1
    
    

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
                    ftemp1=SimpsonIntegration(ff11[i1-1],fL1[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL1[i1][k] = ftemp1*D_2;
                    fL1[i1][ngrid-k] = fL1[i1][k];
                }
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
                    
                    ftemp1=SimpsonIntegration(ff11[i1+1],fR1[i1+1],LLIM,ULIM,kin,kip);
                    
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
                        rho2[i2][k] = 0;
                        for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                        {
                            
                            rho_s1[j2][k] =  ff11[j2][k]*fL1[j2][k]*fR1[j2][k];
                            rho2[i2][k]   =  rho2[i2][k] + rho_s1[j2][k];
                        }
                    }
                    else
                    {
                        rho2[i2][k] = rho1[i2][k];
                    }
                    
                    
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
                            
                            temp0 *= ff11[i2-1][Nz1];
                            fL1[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL1[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp11-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff11[i2+1][Nz2];
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
                    
                    ftemp1=SimpsonIntegration(fL1[i1][k],fR1[i1][k],orient1,kin,kip);
                    rho_s1[i1][k] = ff11[i1][k]*ftemp1*D_2;
                }//loop for grids
            }//loop for beads
            
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    if(k<=ngrid_b)
                    {
                        rho2[i2][k] = 0;
                        for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                        {
                            rho2[i2][k]   =  rho2[i2][k] + rho_s1[j2][k];
                        }
                    }
                    else
                    {
                        rho2[i2][k] = rho1[i2][k];
                    }
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
                        f_0[delta_Z0]  = ff11[i1-1][k2];
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
                        delta_Z0 = round((k1 - k)/nDR0) + midn;//correct
                        
                        kin1 = k1 - mIntB1;
                        kip1 = k1 + mIntB1;
                        if(kin1 < LLIM) kin1 = LLIM;
                        if(kip1 > ULIM) kip1 = ULIM;
                        
                        lowB = round((kin1-k1)/nDR1) + midn;
                        uppB = round((kip1-k1)/nDR1) + midn;
                        for(int k2=kin1; k2<=kip1; ++k2)
                        {//loop for Z_2
                            delta_Z1 = round((k2-k1)/nDR1) + midn;
                            
                            f_2[delta_Z1] = fR1[i1+1][k1][delta_Z1];
                            f_1[delta_Z1] = ff11[i1+2][k2]*BesselZero[0][delta_Z0][delta_Z1];
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
                        f_1[delta_Z0]  = ff11[i1+1][k1];
                        f_2[delta_Z0]  = fL1[i1][k][delta_Z0]*fR1[i1][k][delta_Z0];
                    }
                    
                    ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s1[i1][k] = ff11[i1][k]*ftemp1*D_2;
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
                    f_1[delta_Z0]  = ff11[mp12][k1];
                    f_2[delta_Z0]  = fL1[mp12][k1][delta_Z0];
                }
                
                ftemp1 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s1[mp11][k] = ff11[mp11][k]*ftemp1*D_2;
            }
            ///////////////////////////////////////////////////////////////
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[0]; ++i2)
                {
                    if(k<=ngrid_b)
                    {
                        rho2[i2][k] = 0;
                        for(short j2=MB1[i2]; j2<MB1[i2+1]; ++j2)
                        {
                            rho2[i2][k]   =  rho2[i2][k] + rho_s1[j2][k];
                        }
                    }
                    else
                    {
                        rho2[i2][k] = rho1[i2][k];
                    }
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
                    
                    ftemp2=SimpsonIntegration(ff22[i1-1],fL2[i1-1],LLIM,ULIM,kin,kip);
                    
                    fL2[i1][k] = ftemp2*D_2;
                    fL2[i1][ngrid-k] = fL2[i1][k];
                }
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
                    
                    ftemp2=SimpsonIntegration(ff22[i1+1],fR2[i1+1],LLIM,ULIM,kin,kip);
                    
                    fR2[i1][k] = ftemp2*D_2;
                    fR2[i1][ngrid-k] = fR2[i1][k];
                }
            }
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    if(k<=ngrid_b)
                    {
                        rho2[i0][k] = 0.0;
                        for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                        {
                            rho_s2[j2][k] = ff22[j2][k]*fL2[j2][k]*fR2[j2][k];
                            rho2[i0][k]   = rho2[i0][k] + rho_s2[j2][k];
                        }
                    }
                    else
                    {
                        rho2[i0][k] = rho1[i0][k];
                    }
                    
                    
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
                            
                            temp0 *= ff22[i2-1][Nz1];
                            fL2[i2][Nz2][k2]       = temp0;
                            if(Nz0 != ngNz0) fL2[i2][ngrid-Nz2][k3] = temp0;
                        }
                        
                        temp0 = 1;
                        for(short i2=(mp22-1); i2>=0; --i2)
                        {
                            Nz1    = N_z[i2];
                            Nz2    = N_z[i2+1];
                            
                            temp0 *= ff22[i2+1][Nz2];
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
                    rho_s2[i1][k] = ff22[i1][k]*ftemp2*D_2;
                }//loop for grids
            }//loop for beads
            
            
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    if(k<=ngrid_b)
                    {
                        rho2[i0][k] = 0.0;
                        for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                        {
                            rho2[i0][k]   = rho2[i0][k] + rho_s2[j2][k];
                        }
                    }
                    else
                    {
                        rho2[i0][k] = rho1[i0][k];
                    }
                    
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
                D_2    = 0.5/nDR1;
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
                        f_0[delta_Z0]  = ff22[i1-1][k2];
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
                D_2    = 0.5/nDR1;
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
                            f_1[delta_Z1] = ff22[i1+2][k2]*BesselZero[1][delta_Z0][delta_Z1];
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
                        f_1[delta_Z0]  = ff22[i1+1][k1];
                        f_2[delta_Z0]  = fL2[i1][k][delta_Z0]*fR2[i1][k][delta_Z0];
                    }
                    
                    ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                    rho_s2[i1][k] = ff22[i1][k]*ftemp2*D_2;
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
                    f_1[delta_Z0]  = ff22[mp12][k1];
                    f_2[delta_Z0]  = fL2[mp12][k1][delta_Z0];
                }
                
                ftemp2 = SimpsonIntegration(f_1,f_2,lowB,uppB,lowB,uppB);
                rho_s2[mp22][k] = ff22[mp22][k]*ftemp2*D_2;
            }
            ///////////////////////////////////////////////////////////////
            
            for(int k=LLIM; k<=ngrid_m; ++k)
            {
                for(short i2=0; i2<nb[1]; ++i2)
                {
                    i0 = i2 + nb[0];
                    if(k<=ngrid_b)
                    {
                        rho2[i0][k] = 0.0;
                        for(short j2=MB2[i2]; j2<MB2[i2+1]; ++j2)
                        {
                            rho2[i0][k]   = rho2[i0][k] + rho_s2[j2][k];
                        }
                    }
                    else
                    {
                        rho2[i0][k] = rho1[i0][k];
                    }
                    
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
    
    
    fsigma = 0;
    for(short i=0; i<hspecies; ++i)
    {
        if(Z[i] != 0)
        {
            //temp    = RombergIntegration(rho2[i],0,ngridm1,0,ngrid_m);
            temp    = SimpsonIntegration(rho2[i],0,ngrid_m,0,ngrid_b);
            fsigma += temp*Z[i];
        }
    }
    
    if(rhoBM1 > 1E-10)
    {
        for(short i=0; i<mp1; ++i)
        {
            delete [] ff11[i];
            delete [] rho_s1[i];
        }
        delete [] ff11;
        delete [] rho_s1;
    }
    
    if(rhoBM2 > 1E-10)
    {
        for(short i=0; i<mp2; ++i)
        {
            delete [] ff22[i];
            delete [] rho_s2[i];
        }
        delete [] ff22;
        delete [] rho_s2;
    }
    delete [] Psi_delta;
    
    
    return fsigma;
    
}
