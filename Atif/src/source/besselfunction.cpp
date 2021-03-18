//*******************for charge Bessel Function***************//
#include "clibrary.h"
#include "constantnum.h"
#include "besselfunction.h"

//using namespace std;
extern double dr;


void BesselFunction(int nBessel,double* epsilon_b,double*** BesselZero)
{
    int    midn,n_11,nz_1,nz_2;
    double temp1,temp2,temp3,temp4,temp5,epsilon,dr2;
    
    dr2  = dr*dr;
    midn = nBessel/2;
    n_11 = round(1.0/dr2);
    for(short p=0; p<2; ++p)
    {
        epsilon = epsilon_b[p];
        for(int i=0; i<nBessel; ++i)
        {
            nz_1 = i-midn;
            temp1= dr*epsilon*sqrt(n_11-nz_1*nz_1);
            for(int j=0; j<nBessel; ++j)
            {
                nz_2 = j-midn;
                temp2= dr*sqrt(n_11-nz_2*nz_2);
                temp3= temp1*temp2;
                temp4= exp(epsilon*dr2*(nz_1*nz_2-n_11));
                temp5= Bessel_Zero(temp3);
                
                BesselZero[p][i][j] = temp4*temp5;
            }
            
            
        }
    }
}


double Bessel_Zero(double xx)
{
    int    k;
    double temp_x,temp_2,temp_f,results0,temp;
    double Error_T,error,Max_k;
    
    results0 = 1;
    if(xx == 0) return results0;
    
    Error_T = 1.0E-12;
    Max_k   = 1.0E5;
    error   = 1;
    temp_x  = xx;
    temp_2  = 2;
    temp_f  = 1;
    temp    = temp_x/(temp_2*temp_f);
    k       = 1;
    results0= temp*temp + 1;
    do
    {
        ++k;
        temp_x *= xx;
        temp_2 *= 2;
        temp_f *= k;
        
        temp    = temp_x/(temp_2*temp_f);
        error   = temp*temp;
        results0 += error;
        
    }while((error > Error_T) && (k < Max_k));
    
    if(k >= Max_k)
    {
        std::cerr<<"we cannot get converge value for Bessel function"<<std::endl;
        exit(0);
    }
    
    return results0;
}

double Bessel_Zero_Int(double xx)
{
    int     n_pi,n_pi0;
    double  dr_pi,temp0,temp1,results0,r_pi;
    double* exp_cos;
    
    dr_pi = 1.0E-6;
    n_pi  = round(Pi/dr_pi);
    
    exp_cos = new double[n_pi+1]();
    
    for(int i=0; i<=n_pi; ++i)
    {
        r_pi  = i*dr_pi;
        temp0 = cos(r_pi);
        temp1 = temp0*xx;
        exp_cos[i] = exp(temp1);
    }
    
    n_pi0 = n_pi;
    if((n_pi+1)%2 == 0) n_pi0 = n_pi -1;
    
    results0  = 0;
    results0 += (exp_cos[0] + exp_cos[n_pi0]);
    results0 += (exp_cos[n_pi0-1]*4.0);
    for(int k=1; k<(n_pi0-2); k=(k+2))
    {
        results0 += (exp_cos[k]*4.0);
        results0 += (exp_cos[k+1]*2.0);
    }
    results0 = results0*dr_pi/3.0;
    
    if(n_pi0 != n_pi) results0 += ((exp_cos[n_pi0] + exp_cos[n_pi])*0.5*dr_pi);
    results0 = results0/Pi;
    
    delete [] exp_cos;
    return results0;
}
