//***calculating the screening parameter in MSA***//
#include "clibrary.h"
#include "constantnum.h"
#include "calculategama.h"

extern double errTol;
extern short  nspecies;
extern double BJ;

double CalculateGama(float* Z,float* D,double** QI_k)
{
    int    iter;
    short  N,hspecies;
    double K_Debye,DEB0,eta,PiBJ,error,iter_MAX,D_av;
    double gama0,gama1,gama2,gama_temp,zzg0,zzg1;
    double Chi,gamma_s,temp1,temp2,temp3,temp4,temp5,temp6,temp7;
    double *QI0_k, *QI1_k, *QI2_k, *QI3_k, *GammaS, *GammaS1, *GammaS2;
    
    iter_MAX = 1E5;
    error= errTol*0.01;
    PiBJ = Pi*BJ;
    
    hspecies= nspecies - 1;
    
    QI0_k   = new double[hspecies]();
    QI1_k   = new double[hspecies]();
    QI2_k   = new double[hspecies]();
    QI3_k   = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    GammaS2 = new double[hspecies]();
    
    
    for(short j=0; j<hspecies; ++j)
    {
        QI0_k[j] = QI_k[j][0];
        QI1_k[j] = QI_k[j][1];
        QI2_k[j] = QI_k[j][2];
        QI3_k[j] = QI_k[j][3];
    }
    
    D_av = 0.0;
    N    = 0;
    for(short j=0; j<hspecies; ++j)
    {
        if((Z[j] != 0) && (QI0_k[j] != 0))
        {
            D_av += D[j];
            ++N;
        }
    }
    K_Debye = 0;
    for(short j=0; j<hspecies; ++j)
    {
        if(Z[j] != 0) K_Debye += (QI0_k[j]*Z[j]*Z[j]);
    }
    K_Debye= sqrt(PiBJ*K_Debye);
    DEB0   = 0;
    
    if(N != 0)
    {
        D_av = D_av/N;
        DEB0 = (sqrt(0.25+K_Debye*D_av)-0.5)/D_av;
    }
    
    eta = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        eta += QI3_k[i];
    }
    
    
    if(DEB0 == 0)
    {
        delete [] QI0_k;
        delete [] QI1_k;
        delete [] QI2_k;
        delete [] QI3_k;
        delete [] GammaS;
        delete [] GammaS1;
        delete [] GammaS2;
        
        return DEB0;
    }
    
    
    gama0= DEB0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gama0*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        GammaS2[i] = GammaS1[i]*GammaS1[i];
    }
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    temp5 = 0;
    temp6 = 0;
    temp7 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        temp1 += (Z[i]*Z[i]*QI0_k[i]*GammaS2[i]);
        temp2 += (Z[i]*QI1_k[i]*GammaS2[i]);
        temp3 += (Z[i]*QI1_k[i]*GammaS1[i]);
        temp4 += (QI3_k[i]*GammaS2[i]*GammaS[i]);
        temp5 += (Z[i]*QI2_k[i]*GammaS2[i]*(2+1.0/GammaS[i]));
        temp6 += (Z[i]*QI2_k[i]*GammaS1[i]/GammaS[i]);
        temp7 += (QI3_k[i]*GammaS1[i]);
    }
    
    Chi     = 1.0/(1-eta+3*temp7);
    gamma_s = temp6*Chi;
    
    
    gama1 = PiBJ*(temp1+gamma_s*temp2+temp3*Chi*(3*gamma_s*temp4-temp5));
    gama1 = sqrt(gama1);
    zzg0  = gama1;
    
    iter = 0;
    while(fabs(gama1-gama0) > (error*DEB0))
    {
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gama1*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            GammaS2[i] = GammaS1[i]*GammaS1[i];
        }
        
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        temp4 = 0;
        temp5 = 0;
        temp6 = 0;
        temp7 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            temp1 += (Z[i]*Z[i]*QI0_k[i]*GammaS2[i]);
            temp2 += (Z[i]*QI1_k[i]*GammaS2[i]);
            temp3 += (Z[i]*QI1_k[i]*GammaS1[i]);
            temp4 += (QI3_k[i]*GammaS2[i]*GammaS[i]);
            temp5 += (Z[i]*QI2_k[i]*GammaS2[i]*(2+1.0/GammaS[i]));
            temp6 += (Z[i]*QI2_k[i]*GammaS1[i]/GammaS[i]);
            temp7 += (QI3_k[i]*GammaS1[i]);
        }
        
        Chi     = 1.0/(1-eta+3*temp7);
        gamma_s = temp6*Chi;
        
        
        gama2 = PiBJ*(temp1+gamma_s*temp2+temp3*Chi*(3*gamma_s*temp4-temp5));
        gama2 = sqrt(gama2);
        zzg1  = gama2;
        
        gama_temp = gama1 + (gama1-gama0)*(gama1-zzg1)/(gama0-zzg0-gama1+zzg1);
        
        gama0 = gama1;
        gama1 = gama_temp;
        zzg0  = zzg1;
        
        ++iter;
        if(iter > iter_MAX)
        {
            std::cerr<<"Exceed the iterative maximum in calculategamma.cpp"<<std::endl;
            exit(0);
        }
        
    }
    
    
    delete [] QI0_k;
    delete [] QI1_k;
    delete [] QI2_k;
    delete [] QI3_k;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] GammaS2;
    
    return gama1;
}



double CalculateGama(float* Z,float* D,double** H,double* QI_k)
{
    int    iter;
    short  N,hspecies;
    double K_Debye,DEB0,eta,PiBJ,error,iter_MAX,D_av;
    double gama0,gama1,gama2,gama_temp,zzg0,zzg1;
    double Chi,gamma_s,temp1,temp2,temp3,temp4,temp5,temp6,temp7;
    double *QI0_k, *QI1_k, *QI2_k, *QI3_k, *GammaS, *GammaS1, *GammaS2;
    
    iter_MAX = 1E5;
    error= errTol*0.01;
    PiBJ = Pi*BJ;
    
    hspecies= nspecies - 1;
    
    QI0_k   = new double[hspecies]();
    QI1_k   = new double[hspecies]();
    QI2_k   = new double[hspecies]();
    QI3_k   = new double[hspecies]();
    GammaS  = new double[hspecies]();
    GammaS1 = new double[hspecies]();
    GammaS2 = new double[hspecies]();
    
    
    for(short j=0; j<hspecies; ++j)
    {
        QI0_k[j] = QI_k[j];
        QI1_k[j] = QI_k[j]*H[0][j];
        QI2_k[j] = QI_k[j]*H[1][j];
        QI3_k[j] = QI_k[j]*H[2][j];
    }
    
    D_av = 0.0;
    N    = 0;
    for(short j=0; j<hspecies; ++j)
    {
        if((Z[j] != 0) && (QI0_k[j] != 0))
        {
            D_av += D[j];
            ++N;
        }
    }
    K_Debye = 0;
    for(short j=0; j<hspecies; ++j)
    {
        if(Z[j] != 0) K_Debye += (QI0_k[j]*Z[j]*Z[j]);
    }
    K_Debye= sqrt(PiBJ*K_Debye);
    DEB0   = 0;
    
    if(N != 0)
    {
        D_av = D_av/N;
        DEB0 = (sqrt(0.25+K_Debye*D_av)-0.5)/D_av;
    }
    
    eta = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        eta += QI3_k[i];
    }
    
    
    if(DEB0 == 0)
    {
        delete [] QI0_k;
        delete [] QI1_k;
        delete [] QI2_k;
        delete [] QI3_k;
        delete [] GammaS;
        delete [] GammaS1;
        delete [] GammaS2;
        
        return DEB0;
    }
    
    
    gama0= DEB0;
    for(short i=0; i<hspecies; ++i)
    {
        GammaS[i]  = gama0*D[i];
        GammaS1[i] = 1.0/(1 + GammaS[i]);
        GammaS2[i] = GammaS1[i]*GammaS1[i];
    }
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    temp5 = 0;
    temp6 = 0;
    temp7 = 0;
    for(short i=0; i<hspecies; ++i)
    {
        temp1 += (Z[i]*Z[i]*QI0_k[i]*GammaS2[i]);
        temp2 += (Z[i]*QI1_k[i]*GammaS2[i]);
        temp3 += (Z[i]*QI1_k[i]*GammaS1[i]);
        temp4 += (QI3_k[i]*GammaS2[i]*GammaS[i]);
        temp5 += (Z[i]*QI2_k[i]*GammaS2[i]*(2+1.0/GammaS[i]));
        temp6 += (Z[i]*QI2_k[i]*GammaS1[i]/GammaS[i]);
        temp7 += (QI3_k[i]*GammaS1[i]);
    }
    
    Chi     = 1.0/(1-eta+3*temp7);
    gamma_s = temp6*Chi;
    
    
    gama1 = PiBJ*(temp1+gamma_s*temp2+temp3*Chi*(3*gamma_s*temp4-temp5));
    gama1 = sqrt(gama1);
    zzg0  = gama1;
    
    iter = 0;
    while(fabs(gama1-gama0) > (error*DEB0))
    {
        for(short i=0; i<hspecies; ++i)
        {
            GammaS[i]  = gama1*D[i];
            GammaS1[i] = 1.0/(1 + GammaS[i]);
            GammaS2[i] = GammaS1[i]*GammaS1[i];
        }
        
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        temp4 = 0;
        temp5 = 0;
        temp6 = 0;
        temp7 = 0;
        for(short i=0; i<hspecies; ++i)
        {
            temp1 += (Z[i]*Z[i]*QI0_k[i]*GammaS2[i]);
            temp2 += (Z[i]*QI1_k[i]*GammaS2[i]);
            temp3 += (Z[i]*QI1_k[i]*GammaS1[i]);
            temp4 += (QI3_k[i]*GammaS2[i]*GammaS[i]);
            temp5 += (Z[i]*QI2_k[i]*GammaS2[i]*(2+1.0/GammaS[i]));
            temp6 += (Z[i]*QI2_k[i]*GammaS1[i]/GammaS[i]);
            temp7 += (QI3_k[i]*GammaS1[i]);
        }
        
        Chi     = 1.0/(1-eta+3*temp7);
        gamma_s = temp6*Chi;
        
        
        gama2 = PiBJ*(temp1+gamma_s*temp2+temp3*Chi*(3*gamma_s*temp4-temp5));
        gama2 = sqrt(gama2);
        zzg1  = gama2;
        
        gama_temp = gama1 + (gama1-gama0)*(gama1-zzg1)/(gama0-zzg0-gama1+zzg1);
        
        gama0 = gama1;
        gama1 = gama_temp;
        zzg0  = zzg1;
        
        ++iter;
        if(iter > iter_MAX)
        {
            std::cerr<<"Exceed the iterative maximum in calculategamma.cpp"<<std::endl;
            exit(0);
        }
        
    }
    
    
    delete [] QI0_k;
    delete [] QI1_k;
    delete [] QI2_k;
    delete [] QI3_k;
    delete [] GammaS;
    delete [] GammaS1;
    delete [] GammaS2;
    
    return gama1;
}

