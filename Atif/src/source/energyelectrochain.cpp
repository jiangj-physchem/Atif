//**calculating inhomogeneous chemical potential for electrostatic correlation**//
#include "clibrary.h"
#include "weighteddensity.h"
#include "simpsonintegration.h"
#include "constantnum.h"
#include "calculategama.h"
#include "energyelectrochain.h"

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
void EnergyElectroChain(double gammab,int* LLI,int* ULI,float* D,float* Z,double* rhoB,
                        double** ATT,double** rho,double& EN_EC)
{
    //LLI: the lower limit of intergral: LLI[i]=round(D[i]*0.5/dr)
    //ULI: the upper limit of intergral: ULI[i]=round((size-D[i]*0.5)/dr)
    double*  f_EN;
    double*  B;
    double*  N0I;
    double*  N1I;
    double*  N2I;
    double*  N3I;
    double*  Z_eff;
    double*  GammaS;
    double*  GammaS1;
    double** H;
    double** NI;
    double** sIJ;
    double** Y;
    double   temp1,temp2,temp3,temp4,N3,gammar,gamma_k,BJGA,Chi,gamma_s;

    short    hspecies,nblocks,p_max,i1;
    
    double   rhoBM1,rhoBM2,AN0M1,AN0M2,tempm1,En_C,En_E;
    int      igammar,k1,n_m_gama,n_gama,RMAX,b_max;
    
    hspecies= nspecies - 1;
    nblocks = nb[0] + nb[1];
    RMAX    = DMAX/2;
    

    N0I   = new double[hspecies]();
    N1I   = new double[hspecies]();
    N2I   = new double[hspecies]();
    N3I   = new double[hspecies]();
    B     = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    Z_eff   = new double[nblocks]();

    NI     = new double*[hspecies]();
    H      = new double*[4]();
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    
    sIJ     = new double*[nblocks]();
    Y       = new double*[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        Y[i]    = new double[nblocks]();
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
    
    
    gammar = 0.5/gammab;
    igammar= round(gammar/dr);
    
    b_max  = RMAX + igammar;
    
    n_m_gama = ngrid_m + b_max;
    n_gama   = ngrid + b_max*2;
    
    f_EN = new double[n_gama+1]();
    for(short i=0; i<hspecies; ++i)
    {
        NI[i]  = new double[4]();
        
        
        B[i] = 0.5*D[i] + igammar*dr;
        temp1   = 0.5*D[i]/B[i];
        H[0][i] = 1.0;
        H[1][i] = temp1;
        H[2][i] = temp1*temp1;
        H[3][i] = H[2][i]*temp1;
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
        
        for(short i=0; i<nblocks; ++i)
        {
            Z_eff[i] = (GammaS1[i]*(Z[i]-gamma_s*D[i]*0.5*GammaS[i]));
        }
        
        En_E = gamma_k*gamma_k*gamma_k/(3*Pi) - (temp4+gamma_s*temp1)*BJGA;
    
        if(N3 < 1E-20)
        {
            f_EN[k] = 0;
        }
        else //small loop
        {
            //Start: derivation of the hard sphere contribution
            if(N3 > 0.99) N3 = 0.99;
        
            ///////////////electrostatic correlation for chain connectivity -start /////////////
            //Start: chemical potential for chain connectivity
            //Here Y is log(Y) actually
            for(short p=0; p<p_max; ++p)
            {
                i1 = p*nb[1];
                for(short i=(p*nb[0]); i<(nb[0]+i1); ++i)
                {
                    
                    Y[i][i] = BJ*(Z_eff[i]*Z_eff[i] - Z[i]*Z[i])/sIJ[i][i];
                    
                    
                    if(i < (nb[0]+i1-1))
                    {
                        Y[i][i+1] = BJ*(Z_eff[i]*Z_eff[i+1] - Z[i]*Z[i+1])/sIJ[i][i+1];
                    }
                }
            }
            
            
            //The derivation of excess free energy density for chain connectivity
            En_C = 0;
            if(rhoBM1 > 1.0e-12)
            {
                AN0M1 = 0.0;
                for(short i=0; i<nb[0]; ++i)
                {
                    AN0M1 = AN0M1 + N0I[i];
                }
                
                
                tempm1 = 0;
                for(short i=0; i<nb[0]; ++i)
                {
                    tempm1 += (ATT[i][i]*Y[i][i]);
                    if(i < (nb[0]-1)) tempm1 += (ATT[i][i+1]*Y[i][i+1]);
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
                
                tempm1 = 0;
                for(short i=nb[0]; i<nblocks; ++i)
                {
                    tempm1 += (ATT[i][i]*Y[i][i]);
                    if(i < (nblocks-1)) tempm1 += (ATT[i][i+1]*Y[i][i+1]);
                }
                
                
                En_C += (tempm1*AN0M2);
            }
            ///////////////electrostatic correlation for chain connectivity -end   /////////////
            
            f_EN[k] = En_E + En_C;
            
            f_EN[n_gama-k] = f_EN[k];
            
        }//small loop
        
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    EN_EC = SimpsonIntegration(f_EN,0,n_gama,0,n_gama);
    
    
    delete [] N0I;
    delete [] N1I;
    delete [] N2I;
    delete [] N3I;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] Z_eff;
    delete [] f_EN;
    delete [] B;
    
    for(short i=0; i<4; ++i)
    {
        delete [] H[i];
    }
    delete [] H;
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] NI[i];
    }
    delete [] NI;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] Y[i];
        delete [] sIJ[i];
    }
    delete [] Y;
    delete [] sIJ;

}



