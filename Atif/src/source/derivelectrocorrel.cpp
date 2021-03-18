//**calculating inhomogeneous chemical potential for electrostatic correlation**//
#include "clibrary.h"
#include "weighteddensity.h"
#include "simpsonintegration.h"
#include "constantnum.h"
#include "calculategama.h"
#include "derivelectrocorrel.h"

extern double dr;
extern double BJ;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int DMAX;
extern short  nspecies;
//using namespace std;

//////////////////////////////////////// 1, R, S, V ///////////////////////////////////////////////////
void DerivElectroCorrel(double gammab,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                        double** ATT,double** rho,double** DES)
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
    double*  Z_eff;
    double*  gama_s_j;
    double*  GammaS;
    double*  GammaS1;
    double*  GammaS2;
    double*  GammaS3;
    double** DDC;
    double** H;
    double** NI;
    double** DD1;
    double** DD2;
    double** sIJ;
    double** Y;
    double**  D_Z_eff0;
    double*** DY1;
    double*** DY2;
    double*** DDGg;
    double   temp1,temp2,temp3,N3,Pi2,Pi4,gammar,gamma_k,BJGA,PI_BJ,Chi,gamma_s,DDH1,DDH2;
    double   f_gq,temp4,temp5,temp6,temp7,temp8,d_gama_s,d_chi,d_f_gq,DDC1,DDC2;
    double   temp9,temp10,temp11,temp12,AA1,BB1,temp_d_gs,temp_d_f,temp_d_chi;
    short    hspecies,nblocks,p_max,i1;
    
    double   R,rin,rip,rhoBM1,rhoBM2,AN0M1,AN0M2,D_Z_eff1,tempm1,tempm2;
    int      kin,kip,igammar,k1,n_m_gama,n_gama,RMAX,b_max;
    
    hspecies= nspecies - 1;
    nblocks = nb[0] + nb[1];
    RMAX    = DMAX/2;
    
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
    GammaS2 = new double[hspecies]();
    GammaS3 = new double[hspecies]();
    Z_eff   = new double[nblocks]();
    gama_s_j= new double[nblocks]();
    DDH   = new double[4]();
    


    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    DD2    = new double*[hspecies]();
    DDC    = new double*[4]();
    
    DDGg    = new double**[2]();
    DDGg[0] = new double*[4]();
    DDGg[1] = new double*[4]();
    D_Z_eff0= new double*[4]();
    for(short i=0; i<4; ++i)
    {
        H[i]  = new double[hspecies]();
        
        DDC[i]= new double[hspecies]();
        
        D_Z_eff0[i]= new double[hspecies]();
        DDGg[0][i] = new double[hspecies]();
        DDGg[1][i] = new double[hspecies]();
        
    }
    
    
    sIJ     = new double*[nblocks]();
    Y       = new double*[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
        sIJ[i]  = new double[nblocks]();
    }
    
    DY1  = new double**[nblocks]();
    DY2  = new double**[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        DY1[i]  = new double*[4]();
        DY2[i]  = new double*[4]();
        for(short j=0; j<4; ++j)
        {
            DY1[i][j]  = new double[hspecies]();
            DY2[i][j]  = new double[hspecies]();
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
    
    
    
    
    Pi2    = 2.0*Pi;
    Pi4    = 2.0*Pi2;
    PI_BJ  = Pi*BJ;
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    
    b_max  = RMAX + igammar;
    
    n_m_gama = ngrid_m + b_max;
    n_gama   = ngrid + b_max*2;
    for(short i=0; i<hspecies; ++i)
    {
        DD1[i] = new double[n_gama+1]();
        DD2[i] = new double[n_gama+1]();
        
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
    for(int k=0; k<=n_m_gama; ++k)
    {
        k1 = k - b_max;
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
        temp4 = 0;
        temp5 = 0;
        temp6 = 0;
        temp7 = 0;
        temp8 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gamma_k*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            GammaS2[i] = GammaS1[i]*GammaS1[i];
            GammaS3[i] = GammaS2[i]*GammaS1[i];
            
            temp1 += (Z[i]*N1I[i]*GammaS1[i]);
            temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
            temp3 += (N3I[i]*GammaS1[i]);
            
            temp4 += (Z[i]*N2I[i]*GammaS2[i]);
            temp5 += (GammaS[i]*N3I[i]*GammaS2[i]);
            
            temp6 += (D[i]*N3I[i]*GammaS2[i]);
            temp7 += (D[i]*N3I[i]*GammaS[i]*GammaS3[i]);
            temp8 += (D[i]*Z[i]*N2I[i]*GammaS3[i]);
        }
        Chi     = 1.0/(1-N3+3*temp3); //eq.(15.16)
        gamma_s = temp2*Chi;
        
        for(short i=0; i<nblocks; ++i)
        {
            gama_s_j[i] = gamma_s*D[i]*0.5;
            
            Z_eff[i] = (GammaS1[i]*(Z[i]-gama_s_j[i]*GammaS[i])); //eq.(15.12)
        }
        
        
        f_gq = 3*gamma_s*Chi*temp5 - Chi*temp4 - gamma_s; //eq.(15.19)
        d_gama_s = f_gq/gamma_k; //eq.(15.24)
        d_chi    = 3*Chi*Chi*temp6; //eq.(15.25)
        d_f_gq   = 3*(gamma_s*d_chi+Chi*d_gama_s)*temp5 + 3*gamma_s*Chi*temp6
                   - 6*gamma_s*Chi*temp7 - d_chi*temp4 + 2*Chi*temp8 - d_gama_s;  //eq.(15.23)
        
        if(N3 < 1E-20)
        {
            for(short j1=0; j1<hspecies; ++j1)
            {
                DD1[j1][k] = 0;
                DD2[j1][k] = 0;
                
                DD1[j1][n_gama-k] = DD1[j1][k];
                DD2[j1][n_gama-k] = DD2[j1][k];
            }
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
            
            temp9  = 0;
            temp10 = 0;
            temp11 = 0;
            temp12 = 0;
            for(short i=0; i<hspecies; ++i)
            {
                temp9 += (D[i]*N0I[i]*Z[i]*Z[i]*GammaS3[i]);
                temp10+= (D[i]*N1I[i]*Z[i]*GammaS3[i]);
                temp11+= (D[i]*N1I[i]*Z[i]*GammaS2[i]);
                temp12+= (N1I[i]*Z[i]*GammaS2[i]);
            }
            BB1 = 2*gamma_k + PI_BJ*(2*temp9+2*gamma_s*temp10-temp1*d_f_gq+f_gq*temp11-d_gama_s*temp12);
            
            for(short i=0; i<hspecies; ++i)
            {
                
                ///////////////electrostatic correlation for small ions /////////////////////
                DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
                DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
                DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
                DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
                
                
                DDH1   = BJGA*(DDH[2]*H[2][i] + DDH[1]*RI[i]*H[1][i] + RI2[i]*DDH[0]);
                DDH2   = BJGA*DDH[3];
                ///////////////electrostatic correlation for small ions /////////////////////
                
                
                ///////////////electrostatic correlation for chain connectivity /////////////
                AA1 = PI_BJ*Z[i]*Z[i]*GammaS2[i];
                DDGg[0][0][i]  = AA1/BB1;              //Eq.15.26 in our note
                DDGg[1][0][i]  = d_gama_s*DDGg[0][0][i];    //Eq.15.22 in our note
                
                AA1 = PI_BJ*Z[i]*(gamma_s*GammaS2[i]+f_gq*GammaS1[i]);
                DDGg[0][1][i]  = AA1/BB1;              //Eq.15.32 in our note
                DDGg[1][1][i]  = d_gama_s*DDGg[0][1][i];    //Eq.15.31 in our note
                
                temp_d_gs= Z[i]*Chi*GammaS1[i]/GammaS[i]; //Eq.15.36 in our note
                temp_d_f = 3*Chi*temp_d_gs*temp5 - Z[i]*Chi*GammaS2[i] - temp_d_gs;  //Eq.15.39 in our note
                AA1 = PI_BJ*(temp1*temp_d_f + temp12*temp_d_gs);
                DDGg[0][2][i]  = AA1/BB1;              //Eq.15.40 in our note
                DDGg[1][2][i]  = temp_d_gs + d_gama_s*DDGg[0][2][i]; //Eq.15.38 in our note
                
                temp_d_gs = Chi*gamma_s*GammaS1[i]*(GammaS[i]-2);
                temp_d_chi= Chi*Chi*GammaS1[i]*(GammaS[i]-2);
                temp_d_f  = 3*(Chi*temp_d_gs+gamma_s*temp_d_chi)*temp5 + 3*Chi*gamma_s*GammaS2[i]*GammaS[i]
                            - temp_d_chi*temp4 - temp_d_gs;
                
                AA1 = PI_BJ*(temp1*temp_d_f + temp12*temp_d_gs);
                DDGg[0][3][i]  = AA1/BB1;              //Eq.15.49 in our note
                DDGg[1][3][i]  = temp_d_gs + d_gama_s*DDGg[0][3][i]; //Eq.15.46 in our note
                ///////////////electrostatic correlation for chain connectivity /////////////
                
                DD1[i][k] = DDH1;
                DD2[i][k] = DDH2;
                
                DD1[i][n_gama-k] = DDH1;
                DD2[i][n_gama-k] = DDH2;
            }
            
            
            
            ///////////////electrostatic correlation for chain connectivity -start /////////////
            //Start: chemical potential for chain connectivity
            //Here Y is log(Y) actually
            for(short p=0; p<p_max; ++p)
            {
                i1 = p*nb[1];
                for(short i=(p*nb[0]); i<(nb[0]+i1); ++i)
                {
                    
                    Y[i][i] = BJ*(Z_eff[i]*Z_eff[i] - Z[i]*Z[i])/sIJ[i][i]; //eq.(15.2)
                    
                    for(short i0=0; i0<4; ++i0)
                    {
                        for(short j0=0; j0<hspecies; ++j0)
                        {
                            D_Z_eff0[i0][j0] = -D[i]*((gama_s_j[i]+Z_eff[i])*GammaS1[i]*DDGg[0][i0][j0]
                                                      + GammaS1[i]*GammaS[i]*0.5*DDGg[1][i0][j0]);//eq.(15.17)
                            DY1[i][i0][j0]   = (2*BJ*Z_eff[i]*D_Z_eff0[i0][j0]/sIJ[i][i]);//eq.(15.6)
                        }
                    }
                    
                    if(i < (nb[0]+i1-1))
                    {
                        Y[i][i+1] = BJ*(Z_eff[i]*Z_eff[i+1] - Z[i]*Z[i+1])/sIJ[i][i+1];//eq.(15.2)
                        
                        for(short i0=0; i0<4; ++i0)
                        {
                            for(short j0=0; j0<hspecies; ++j0)
                            {
                                D_Z_eff1 = -D[i+1]*((gama_s_j[i+1]+Z_eff[i+1])*GammaS1[i+1]*DDGg[0][i0][j0]
                                                    + GammaS1[i+1]*GammaS[i+1]*0.5*DDGg[1][i0][j0]);//eq.(15.17)
                                DY2[i][i0][j0]   = (BJ*(Z_eff[i]*D_Z_eff1 + Z_eff[i+1]*D_Z_eff0[i0][j0])/sIJ[i][i+1]);//eq.(15.6)
                            }
                        }
                        
                    }
                }
            }
            
            
            //The derivation of excess free energy density for chain connectivity
            for(short i=0; i<4; ++i)
            {
                for(short j=0; j<hspecies; ++j)
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
                
                for(short i0=0; i0<4; ++i0)
                {
                    for(short j=0; j<hspecies; ++j)
                    {
                        tempm2 = 0;
                        for(short i=0; i<nb[0]; ++i)
                        {
                            tempm2 += (ATT[i][i]*DY1[i][i0][j]);
                            if(i < (nb[0]-1)) tempm2 += (ATT[i][i+1]*DY2[i][i0][j]);
                        }
                        
                        DDC[i0][j] += AN0M1*tempm2;
                    }
                }
                
                tempm1 = 0;
                for(short i=0; i<nb[0]; ++i)
                {
                    tempm1 += (ATT[i][i]*Y[i][i]);
                    if(i < (nb[0]-1)) tempm1 += (ATT[i][i+1]*Y[i][i+1]);
                }
                
                for(short i=0; i<nb[0]; ++i)
                {
                    DDC[0][i] += tempm1;
                }
            }
            
            
            if(rhoBM2 > 1.0e-12)
            {
                AN0M2 = 0.0;
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    AN0M2 = AN0M2 + N0I[i];
                }
                
                for(short i0=0; i0<4; ++i0)
                {
                    for(short j=0; j<hspecies; ++j)
                    {
                        tempm2 = 0;
                        for(short i=nb[0]; i<nblocks; ++i)
                        {
                            tempm2 += (ATT[i][i]*DY1[i][i0][j]);
                            if(i < (nblocks-1)) tempm2 += (ATT[i][i+1]*DY2[i][i0][j]);
                        }
                        
                        DDC[i0][j] += AN0M2*tempm2;
                    }
                }
                
                tempm1 = 0;
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    tempm1 += (ATT[i][i]*Y[i][i]);
                    if(i < (nblocks-1)) tempm1 += (ATT[i][i+1]*Y[i][i+1]);
                }
                
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    DDC[0][i] += tempm1;
                }
            }
            ///////////////electrostatic correlation for chain connectivity -end   /////////////
            
            if((rhoBM1 + rhoBM2) > 1.0e-12)
            {
                for(short j=0; j<hspecies; ++j)
                {
                    DDC1 = DDC[2][j]*H[2][j] + DDC[1][j]*RI[j]*H[1][j] + RI2[j]*DDC[0][j];
                    DDC2 = DDC[3][j];
                    
                    DD1[j][k] += DDC1;
                    DD2[j][k] += DDC2;
                    
                    DD1[j][n_gama-k] = DD1[j][k];
                    DD2[j][n_gama-k] = DD2[j][k];
                }
            }
            
            
            
            ///////////////electrostatic correlation for small ions /////////////////////
            
            
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //time_start = clock();
    
    for(int i=b_max; i<=n_m_gama; ++i)
    {
        R = i*dr;
        k1= i - b_max;
        for(short j=0; j<hspecies; ++j)
        {
            rin = R - B[j];
            rip = R + B[j];
            kin = round(rin/dr);
            kip = round(rip/dr);
            if(kin < 0) kin = 0;
            if(kip > n_gama) kip = n_gama;
                        
            temp1 = SimpsonIntegration(DD1[j],0,n_gama,kin,kip);
            temp2 = SimpsonIntegration(DD2[j],0,n_gama,kin,kip,B2[j],i);
            
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
    delete [] GammaS2;
    delete [] GammaS3;
    delete [] gama_s_j;
    delete [] Z_eff;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
        delete [] DDC[i];
        
        delete [] DDGg[0][i];
        delete [] DDGg[1][i];
        delete [] D_Z_eff0[i];
    }
    delete [] H;
    delete [] DDC;
    delete [] DDGg[0];
    delete [] DDGg[1];
    delete [] DDGg;
    delete [] D_Z_eff0;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] DD2[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] DD2;
    delete [] NI;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] Y[i];
        delete [] sIJ[i];
    }
    delete [] Y;
    delete [] sIJ;
    
    for(short i=0; i<nblocks; ++i)
    {
        for(short j=0; j<4; ++j)
        {
            delete [] DY1[i][j];
            delete [] DY2[i][j];
        }
        delete [] DY1[i];
        delete [] DY2[i];
    }
    delete [] DY1;
    delete [] DY2;   

}


void DerivElectroCorrel(double* rhoB,double gammab,float* D,float* Z,double** ATT,
                        double* DESB,double& ES_EN)
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
    double*  Z_eff;
    double*  gama_s_j;
    double*  GammaS;
    double*  GammaS1;
    double*  GammaS2;
    double*  GammaS3;
    double** DDC;
    double** NI;
    double*  DDB1;
    double*  DDB2;
    double** H;
    double** DD1;
    double** DD2;
    double** sIJ;
    double** Y;
    double**  D_Z_eff0;
    double*** DY1;
    double*** DY2;
    double*** DDGg;
    
    double   temp1,temp2,temp3,temp4,N3,Pi2,Pi4,gammar,gamma_k,BJGA,gamma_s,Chi;
    double   f_gq,temp5,temp6,temp7,temp8,d_gama_s,d_chi,d_f_gq,DDC1,DDC2,PI_BJ;
    double   temp9,temp10,temp11,temp12,temp13,AA1,BB1,temp_d_gs,temp_d_f,temp_d_chi,En_C;
    short    hspecies,nblocks,p_max,i1;
    
    double   rin,rip,rhoBM1,rhoBM2,AN0M1,AN0M2,D_Z_eff1,tempm1,tempm2;
    int      kin,kip,RMAX,BMAX,BMAX1,igammar;
    
    hspecies= nspecies - 1;
    nblocks = nb[0] + nb[1];
    
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
    GammaS2 = new double[hspecies]();
    GammaS3 = new double[hspecies]();
    Z_eff   = new double[nblocks]();
    gama_s_j= new double[nblocks]();
    
    DDH   = new double[4]();
    
    H      = new double*[4]();
    NI     = new double*[hspecies]();
    DD1    = new double*[hspecies]();
    DD2    = new double*[hspecies]();
    DDC    = new double*[4]();
    
    DDGg    = new double**[2]();
    DDGg[0] = new double*[4]();
    DDGg[1] = new double*[4]();
    D_Z_eff0= new double*[4]();
    
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
        
        DDC[i]= new double[hspecies]();
        
        D_Z_eff0[i]= new double[hspecies]();
        DDGg[0][i] = new double[hspecies]();
        DDGg[1][i] = new double[hspecies]();
    }
    
    sIJ     = new double*[nblocks]();
    Y       = new double*[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
        sIJ[i]  = new double[nblocks]();
    }
    
    DY1  = new double**[nblocks]();
    DY2  = new double**[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        DY1[i]  = new double*[4]();
        DY2[i]  = new double*[4]();
        for(short j=0; j<4; ++j)
        {
            DY1[i][j]  = new double[hspecies]();
            DY2[i][j]  = new double[hspecies]();
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
    
    
    PI_BJ  = Pi*BJ;
    
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
    temp5 = 0;
    temp6 = 0;
    temp7 = 0;
    temp8 = 0;
    temp13= 0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gamma_k*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        GammaS2[i] = GammaS1[i]*GammaS1[i];
        GammaS3[i] = GammaS2[i]*GammaS1[i];
        
        temp1 += (Z[i]*N1I[i]*GammaS1[i]);
        temp2 += (Z[i]*N2I[i]*GammaS1[i]/GammaS[i]);
        temp3 += (N3I[i]*GammaS1[i]);
        temp13+= (Z[i]*Z[i]*N0I[i]*GammaS1[i]);
        
        temp4 += (Z[i]*N2I[i]*GammaS2[i]);
        temp5 += (GammaS[i]*N3I[i]*GammaS2[i]);
        
        temp6 += (D[i]*N3I[i]*GammaS2[i]);
        temp7 += (D[i]*N3I[i]*GammaS[i]*GammaS3[i]);
        temp8 += (D[i]*Z[i]*N2I[i]*GammaS3[i]);
    }
    Chi     = 1.0/(1-N3+3*temp3);
    gamma_s = temp2*Chi;
    
    for(short i=0; i<nblocks; ++i)
    {
        gama_s_j[i] = gamma_s*D[i]*0.5;
        
        Z_eff[i] = (GammaS1[i]*(Z[i]-gama_s_j[i]*GammaS[i]));
    }
    
    
    f_gq = 3*gamma_s*Chi*temp5 - Chi*temp4 - gamma_s;
    d_gama_s = f_gq/gamma_k;
    d_chi    = 3*Chi*Chi*temp6;
    d_f_gq   = 3*(gamma_s*d_chi+Chi*d_gama_s)*temp5 + 3*gamma_s*Chi*temp6
               - 6*gamma_s*Chi*temp7 - d_chi*temp4 + 2*Chi*temp8 - d_gama_s;
    
    
    ES_EN = 0;
    En_C  = 0;
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
        
        temp9  = 0;
        temp10 = 0;
        temp11 = 0;
        temp12 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            temp9 += (D[i]*N0I[i]*Z[i]*Z[i]*GammaS3[i]);
            temp10+= (D[i]*N1I[i]*Z[i]*GammaS3[i]);
            temp11+= (D[i]*N1I[i]*Z[i]*GammaS2[i]);
            temp12+= (N1I[i]*Z[i]*GammaS2[i]);
        }
        BB1 = 2*gamma_k + PI_BJ*(2*temp9+2*gamma_s*temp10-temp1*d_f_gq+f_gq*temp11-d_gama_s*temp12);
        

        for(short i=0; i<hspecies; ++i)
        {
            DDH[0] = -(GammaS1[i]*Z[i]*Z[i]);
            DDH[1] = -(GammaS1[i]*Z[i]*gamma_s);
            DDH[2] = -((temp1*Chi*Z[i]*GammaS1[i])/GammaS[i]);
            DDH[3] = (gamma_s*temp1*Chi*GammaS1[i]*(2-GammaS[i]));
            
            
            DDB1[i] = BJGA*(DDH[2]*H[2][i] + DDH[1]*RI[i]*H[1][i] + RI2[i]*DDH[0]);
            DDB2[i] = BJGA*DDH[3];
            
            ///////////////electrostatic correlation for chain connectivity /////////////
            AA1 = PI_BJ*Z[i]*Z[i]*GammaS2[i];
            DDGg[0][0][i]  = AA1/BB1;              //Eq.15.26 in our note
            DDGg[1][0][i]  = d_gama_s*DDGg[0][0][i];    //Eq.15.22 in our note
            
            AA1 = PI_BJ*Z[i]*(gamma_s*GammaS2[i]+f_gq*GammaS1[i]);
            DDGg[0][1][i]  = AA1/BB1;              //Eq.15.32 in our note
            DDGg[1][1][i]  = d_gama_s*DDGg[0][1][i];    //Eq.15.31 in our note
            
            temp_d_gs= Z[i]*Chi*GammaS1[i]/GammaS[i]; //Eq.15.36 in our note
            temp_d_f = 3*Chi*temp_d_gs*temp5 - Z[i]*Chi*GammaS2[i] - temp_d_gs;  //Eq.15.39 in our note
            AA1 = PI_BJ*(temp1*temp_d_f + temp12*temp_d_gs);
            DDGg[0][2][i]  = AA1/BB1;              //Eq.15.40 in our note
            DDGg[1][2][i]  = temp_d_gs + d_gama_s*DDGg[0][2][i]; //Eq.15.38 in our note
            
            temp_d_gs = Chi*gamma_s*GammaS1[i]*(GammaS[i]-2);
            temp_d_chi= Chi*Chi*GammaS1[i]*(GammaS[i]-2);
            temp_d_f  = 3*(Chi*temp_d_gs+gamma_s*temp_d_chi)*temp5 + 3*Chi*gamma_s*GammaS2[i]*GammaS[i]
                        - temp_d_chi*temp4 - temp_d_gs;
            
            AA1 = PI_BJ*(temp1*temp_d_f + temp12*temp_d_gs);
            DDGg[0][3][i]  = AA1/BB1;              //Eq.15.49 in our note
            DDGg[1][3][i]  = temp_d_gs + d_gama_s*DDGg[0][3][i]; //Eq.15.46 in our note
            
            ///////////////electrostatic correlation for chain connectivity /////////////
        }
        
        
        ///////////////electrostatic correlation for chain connectivity -start /////////////
        //Start: chemical potential for chain connectivity
        //Here Y is log(Y) actually
        for(short p=0; p<p_max; ++p)
        {
            i1 = p*nb[1];
            for(short i=(p*nb[0]); i<(nb[0]+i1); ++i)
            {
                
                Y[i][i] = BJ*(Z_eff[i]*Z_eff[i] - Z[i]*Z[i])/sIJ[i][i];
                
                for(short i0=0; i0<4; ++i0)
                {
                    for(short j0=0; j0<hspecies; ++j0)
                    {
                        D_Z_eff0[i0][j0] = -D[i]*((gama_s_j[i]+Z_eff[i])*GammaS1[i]*DDGg[0][i0][j0] +
                                                  GammaS1[i]*GammaS[i]*0.5*DDGg[1][i0][j0]);
                        DY1[i][i0][j0]   = (2*BJ*Z_eff[i]*D_Z_eff0[i0][j0]/sIJ[i][i]);
                    }
                }
                
                if(i < (nb[0]+i1-1))
                {
                    Y[i][i+1] = BJ*(Z_eff[i]*Z_eff[i+1] - Z[i]*Z[i+1])/sIJ[i][i+1];
                    
                    for(short i0=0; i0<4; ++i0)
                    {
                        for(short j0=0; j0<hspecies; ++j0)
                        {
                            D_Z_eff1 = -D[i+1]*((gama_s_j[i+1]+Z_eff[i+1])*GammaS1[i+1]*DDGg[0][i0][j0] +
                                                GammaS1[i+1]*GammaS[i+1]*0.5*DDGg[1][i0][j0]);
                            DY2[i][i0][j0]   = (BJ*(Z_eff[i]*D_Z_eff1 + Z_eff[i+1]*D_Z_eff0[i0][j0])/sIJ[i][i+1]);
                        }
                    }
                }
            }
        }
        
        
        //The derivation of excess free energy density for chain connectivity
        for(short i=0; i<4; ++i)
        {
            for(short j=0; j<hspecies; ++j)
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
            
            for(short i0=0; i0<4; ++i0)
            {
                for(short j=0; j<hspecies; ++j)
                {
                    tempm2 = 0;
                    for(short i=0; i<nb[0]; ++i)
                    {
                        tempm2 += (ATT[i][i]*DY1[i][i0][j]);
                        if(i < (nb[0]-1)) tempm2 += (ATT[i][i+1]*DY2[i][i0][j]);
                    }
                    
                    DDC[i0][j] += AN0M1*tempm2;
                }
            }
            
            tempm1 = 0;
            for(short i=0; i<nb[0]; ++i)
            {
                tempm1 += (ATT[i][i]*Y[i][i]);
                if(i < (nb[0]-1)) tempm1 += (ATT[i][i+1]*Y[i][i+1]);
            }
            
            for(short i=0; i<nb[0]; ++i)
            {
                DDC[0][i] += tempm1;
            }
            
            En_C += (tempm1*AN0M1);
        }
        
        
        if(rhoBM2 > 1.0e-12)
        {
            AN0M2 = 0.0;
            for(short i=nb[0]; i<nblocks; ++i)
            {
                AN0M2 = AN0M2 + N0I[i];
            }
            
            for(short i0=0; i0<4; ++i0)
            {
                for(short j=0; j<hspecies; ++j)
                {
                    tempm2 = 0;
                    for(short i=nb[0]; i<nblocks; ++i)
                    {
                        tempm2 += (ATT[i][i]*DY1[i][i0][j]);
                        if(i < (nblocks-1)) tempm2 += (ATT[i][i+1]*DY2[i][i0][j]);
                    }
                    
                    DDC[i0][j] += AN0M2*tempm2;
                }
            }
            
            tempm1 = 0;
            for(short i=nb[0]; i<nblocks; ++i)
            {
                tempm1 += (ATT[i][i]*Y[i][i]);
                if(i < (nblocks-1)) tempm1 += (ATT[i][i+1]*Y[i][i+1]);
            }
            
            for(short i=nb[0]; i<nblocks; ++i)
            {
                DDC[0][i] += tempm1;
            }
            
            En_C += (tempm1*AN0M2);
        }
        ///////////////electrostatic correlation for chain connectivity -end   /////////////
        
        
        if((rhoBM1 + rhoBM2) > 1.0e-12)
        {
            for(short j=0; j<hspecies; ++j)
            {
                DDC1 = DDC[2][j]*H[2][j] + DDC[1][j]*RI[j]*H[1][j] + RI2[j]*DDC[0][j];
                DDC2 = DDC[3][j];
                
                DDB1[j] += DDC1;
                DDB2[j] += DDC2;
            }
        }
        
        
        ES_EN = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp13+gamma_s*temp1)*BJGA + En_C;
        
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
        
        DESB[j] = temp1 + temp2;
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
    delete [] GammaS2;
    delete [] GammaS3;
    delete [] gama_s_j;
    delete [] Z_eff;
    delete [] DDB1;
    delete [] DDB2;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
        delete [] DDC[i];
        
        delete [] DDGg[0][i];
        delete [] DDGg[1][i];
        delete [] D_Z_eff0[i];
    }
    delete [] H;
    delete [] DDC;
    delete [] DDGg[0];
    delete [] DDGg[1];
    delete [] DDGg;
    delete [] D_Z_eff0;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] DD1[i];
        delete [] DD2[i];
        delete [] NI[i];
    }
    delete [] DD1;
    delete [] DD2;
    delete [] NI;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] Y[i];
        delete [] sIJ[i];
    }
    delete [] Y;
    delete [] sIJ;
    
    for(short i=0; i<nblocks; ++i)
    {
        for(short j=0; j<4; ++j)
        {
            delete [] DY1[i][j];
            delete [] DY2[i][j];
        }
        delete [] DY1[i];
        delete [] DY2[i];
    }
    delete [] DY1;
    delete [] DY2;

}
