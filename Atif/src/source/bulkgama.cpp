//***********calculating bulk gama****************//
#include "clibrary.h"
#include "bulkgama.h"
#include "constantnum.h"
extern double errTol;
extern double BJ; //The Bjerrum length
extern short  nspecies;

void BulkGama(double* rhoB,float* D,float* Z,double& gama)
{
    double*  X;
    int      iter;
    short    hspecies;
    double   DAV,i_D,kapa0,C,gamai,error;
    double   gama0,gama2,gama_temp,zzg0,zzg1;
    double   temp1,temp2,temp3,iter_MAX,rhot;
    
    iter_MAX = 1.0E5;
    error= errTol*0.01;
    
    hspecies = nspecies - 1;
    
    
    X      = new double[hspecies]();
    
    
    //calculate the bulk gama
    C = 0.0;
    rhot= 0;
    DAV = 0;
    i_D = 0;
    for(short i=0; i<hspecies; ++i)
    {
        rhot= rhot + Z[i]*Z[i]*rhoB[i];
        C = C + Pi*rhoB[i]*D[i]*D[i]*D[i]/6;
        if((Z[i] != 0) && (rhoB[i] != 0))
        {
            i_D = i_D + 1;
            DAV = DAV + D[i];
        }
    }
    C = Pi*0.5/(1.0-C);
    
    kapa0= sqrt(4*Pi*BJ*rhot);
    if(i_D != 0) DAV = DAV/i_D;
    
    gamai = 0.5*kapa0;
    if(i_D != 0) gamai = (sqrt(2*kapa0*DAV+1)-1)/(2*DAV);
    gama0 = gamai;
    
    temp1 = 0.0;
    temp2 = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        temp3 = gama0*D[i] + 1;
        temp1 = temp1 + rhoB[i]*D[i]*Z[i]/temp3;
        temp2 = temp2 + rhoB[i]*D[i]*D[i]*D[i]/temp3;
    }
    temp2 = C*temp2;
    
    for(short i=0; i<hspecies; ++i)
    {
        temp3 = gama0*D[i] + 1;
        X[i] = Z[i]/temp3 - C*D[i]*D[i]*temp1/(temp3*(1.0+temp2));// X[i] is the effective valency
    }
    
    gama = 0.0;
    for(short i=0; i<hspecies; ++i)
    {
        gama = gama + rhoB[i]*X[i]*X[i];
    }
    
    gama = sqrt(Pi*BJ*gama);
    zzg0  = gama;
    //gama  = (1-mixf)*gama + mixf*gama0;
    
    iter  = 0;
    while(fabs(gama-gama0) > (error*gamai))
    {
        temp1 = 0.0;
        temp2 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            temp3 = gama*D[i] + 1;
            temp1 = temp1 + rhoB[i]*D[i]*Z[i]/temp3;
            temp2 = temp2 + rhoB[i]*D[i]*D[i]*D[i]/temp3;
        }
        temp2 = C*temp2;
        
        for(short i=0; i<hspecies; ++i)
        {
            temp3 = gama*D[i] + 1;
            X[i] = Z[i]/temp3 - C*D[i]*D[i]*temp1/(temp3*(1.0+temp2));// X[i] is the effective valency
        }
        
        gama2 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            gama2 = gama2 + rhoB[i]*X[i]*X[i];
        }
        
        gama2 = sqrt(Pi*BJ*gama2);
        zzg1  = gama2;
        //gama  = (1-mixf)*gama + mixf*gama0;
        
        gama_temp = gama + (gama-gama0)*(gama-zzg1)/(gama0-zzg0-gama+zzg1);
        
        gama0 = gama;
        gama  = gama_temp;
        zzg0  = zzg1;
        
        ++iter;
        
        if(iter > iter_MAX)
        {
            std::cerr<<"Exceed the iterative maximum bulkpotentialsolvent code: 1"<<std::endl;
            exit(0);
        }
        
    }
    
    delete [] X;
}



/*
void BulkGama(double* rhoB,float* D,float* Z,double& gama)
{
    double*  X;
    int      iter;
    short    hspecies;
    double   DAV,i_D,kapa0,C,gama0,mixf,gamai;
    double   temp1,temp2,temp3,iter_MAX,rhot;
    
    mixf= 0.1;
    iter_MAX = 1.0E5;
    
    hspecies = nspecies - 1;
    
    
    X      = new double[hspecies]();
    
    
    //calculate the bulk gama
    C = 0.0;
    rhot= 0;
    DAV = 0;
    i_D = 0;
    for(short i=0; i<hspecies; ++i)
    {
        rhot= rhot + Z[i]*Z[i]*rhoB[i];
        C = C + Pi*rhoB[i]*D[i]*D[i]*D[i]/6;
        if((Z[i] != 0) && (rhoB[i] != 0))
        {
            i_D = i_D + 1;
            DAV = DAV + D[i];
        }
    }
    C = Pi*0.5/(1.0-C);
    
    kapa0= sqrt(4*Pi*BJ*rhot);
    if(i_D != 0) DAV = DAV/i_D;
    
    gamai = 0.5*kapa0;
    if(i_D != 0) gamai = (sqrt(2*kapa0*DAV+1)-1)/(2*DAV);
    gama  = gamai;
    iter  = 0;
    do
    {
        temp1 = 0.0;
        temp2 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            temp3 = gama*D[i] + 1;
            temp1 = temp1 + rhoB[i]*D[i]*Z[i]/temp3;
            temp2 = temp2 + rhoB[i]*D[i]*D[i]*D[i]/temp3;
        }
        temp2 = C*temp2;
        
        for(short i=0; i<hspecies; ++i)
        {
            temp3 = gama*D[i] + 1;
            X[i] = Z[i]/temp3 - C*D[i]*D[i]*temp1/(temp3*(1.0+temp2));// X[i] is the effective valency
        }
        
        gama0 = 0.0;
        for(short i=0; i<hspecies; ++i)
        {
            gama0 = gama0 + rhoB[i]*X[i]*X[i];
        }
        
        gama0 = sqrt(Pi*BJ*gama0);
        gama  = (1-mixf)*gama + mixf*gama0;
        
        ++iter;
        
        if(iter > iter_MAX)
        {
            std::cerr<<"Exceed the iterative maximum bulkpotentialsolvent code: 1"<<std::endl;
            exit(0);
        }
        
        std::cout<<"old  iter= "<<iter<<std::endl;
        
    }
    while(fabs(gama0-gama) > (errTol*gamai));
    
    delete [] X;
}
 */
