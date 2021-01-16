//*******************for charge neutrality***************//
#include "clibrary.h"
#include "constantnum.h"
//#include "rombergintegration.h"
#include "simpsonintegration.h"
#include "renormeulerlagrange.h"
#include "renormalization.h"

using namespace std;

extern double dr;
extern int ngrid; //the number of grids
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int ngrid_b;
extern int LLIM;
extern int ULIM;
extern int DMAX;
extern short* mp;//monomer # on copolymer
extern short* nb;//# of blocks in copolymer i;
extern short** mb;//# of monomers in each block
extern short  nspecies;



void RenormDensity(double sigma,double& deltaPhi,double err,double rhoBM1,double rhoBM2,float* D,
                   short* MB1,short* MB2,float* Z,double** ff1,double** ff2,double** rho1,
                   double*** BesselZero,string* MODEL)
{
    
    int     ratio_iter,ratio_err0,MAXR0;
    double  ratio_err,ratio_MAX;
    double  ratio,temp;//,ratio_t;
    double  ratio_l,ratio_r,dphi_L;
    double  f_l,f_r,f;
    double** rho2;
    short   hspecies;
    
    hspecies= nspecies - 1;
    
    f      = 0;
    for(short i=0; i<hspecies; ++i)
    {
        if(Z[i] != 0)
        {
            temp = SimpsonIntegration(rho1[i],0,ngrid_m,0,ngrid_b);
            f   += temp*Z[i];
        }
    }
    f = f + sigma;
    ratio = 0;
    
    if(fabs(f) < 1E-12)
    {
        deltaPhi = ratio;
    }
    else//start renormalize
    {
        ratio_MAX = 1E10;
        ratio_err = 1E-12;
        if(err > 1E-2)
        {
            ratio_err = 1E-10;
        }
        else if(err > 1E-3)
        {
            ratio_err = 1E-11;
        }
        ratio_err0= ratio_err*0.1;
        MAXR0     = DMAX/2 + ngrid_m;
        
        dphi_L    = 0.001;
        
        rho2      = new double*[hspecies]();
        for(short i=0; i<hspecies; ++i)
        {
            rho2[i] = new double[MAXR0+2]();
        }
        
        
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        //advance renormalization
        //renormalize the surface charge density: sigma
        if(f < 0)
        {
            do
            {
                ratio_r = ratio;
                ratio   = ratio_r - dphi_L;
                f=RenormEulerLagrange(ratio,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
                f += sigma;
            }
            while(f < -ratio_err0);
            ratio_l = ratio;
        }
        else
        {
            do
            {
                ratio_l = ratio;
                ratio   = ratio_l + dphi_L;
                f=RenormEulerLagrange(ratio,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
                f += sigma;
            }
            while(f > ratio_err0);
            ratio_r = ratio;
        }
        
        //cout<<"f="<<f<<" ratio_l="<<ratio_l<<" ratio_r="<<ratio_r<<endl;
        
        
        ratio_iter= 0;
        ratio     = (ratio_l + ratio_r)*0.5;
        do
        {
            ++ratio_iter;
            //ratio_t = ratio;
            
            f_l=RenormEulerLagrange(ratio_l,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
            f_r=RenormEulerLagrange(ratio_r,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
            f  =RenormEulerLagrange(ratio,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
            f_l += sigma;
            f_r += sigma;
            f   += sigma;
            
            
            //cout<<"ratio_iter="<<ratio_iter<<" f="<<f<<" ratio_l="<<ratio_l<<" ratio_r="<<ratio_r<<endl;
            if(fabs(f_l) < ratio_err)
            {
                ratio = ratio_l;
                ratio_iter = 0;
                break;
            }
            
            if(fabs(f_r) < ratio_err)
            {
                ratio = ratio_r;
                ratio_iter = 0;
                break;
            }
            
            if(fabs(f) < ratio_err)
            {
                ratio_iter = 0;
                break;
            }
            
            if(f*f_l > 0)
            {
                ratio_l = ratio;
                ratio   = (ratio_l + ratio_r)*0.5;
            }
            else
            {
                ratio_r = ratio;
                ratio   = (ratio_r + ratio_l)*0.5;
            }
        }
        while(ratio_iter<ratio_MAX);
        
        if(ratio_iter >= ratio_MAX)
        {
            cerr<<"we cannot renormalize the charge neutrality for bulk"<<endl;
            exit(0);
            
        }
        
        
        
        //ratio = (sqrt(sigma*sigma+sigmaP*sigmaN)-sigma)/sigmaP;
        for(int k=LLIM; k<=ngrid_b; ++k)
        {
            for(short i=0; i<hspecies; ++i)
            {
                rho1[i][k]= rho2[i][k];
                rho1[i][ngrid-k]= rho1[i][k];
            }
        }
        
        deltaPhi = ratio;
        
        for(short i=0; i<hspecies; ++i)
        {
            delete [] rho2[i];
        }
        delete [] rho2;
        
    }//end renormalize
    
    
    
    
    //charge neutrality
    //////////////////////////////////////////////////////////////////////////////////////////////////
}




void RenormDensity_Newton(double sigma,double& deltaPhi,double err,double rhoBM1,double rhoBM2,float* D,
                   short* MB1,short* MB2,float* Z,double** ff1,double** ff2,double** rho1,
                   double*** BesselZero,string* MODEL)
{
    int     ratio_iter,MAXR0;
    double  ratio_err,ratio_MAX,ratio,temp,dphi_L;
    double  f_0,f_1,f_2,ratio_0,ratio_1,ratio_2,dphi_i,temp_1,temp_2,ratio_temp,f_temp;
    double** rho2;
    short   hspecies;
    
    hspecies= nspecies - 1;
    
    f_0      = 0;
    for(short i=0; i<hspecies; ++i)
    {
        if(Z[i] != 0)
        {
            temp = SimpsonIntegration(rho1[i],0,ngrid_m,0,ngrid_b);
            f_0   += temp*Z[i];
        }
    }
    ratio_0 = 0;
    f_0 = f_0 + sigma;
    
    if(fabs(f_0) < 1E-12)
    {
        deltaPhi = ratio_0;
    }
    else//start renormalize
    {
        ratio_MAX = 1E6;
        ratio_err = 1E-14;
        if(err > 1E-2)
        {
            ratio_err = 1E-12;
        }
        else if(err > 1E-3)
        {
            ratio_err = 1E-13;
        }
        
        MAXR0     = DMAX/2 + ngrid_m;
        dphi_L    = 0.01;
        
        rho2      = new double*[hspecies]();
        for(short i=0; i<hspecies; ++i)
        {
            rho2[i] = new double[MAXR0+2]();
        }
        
        
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        //advance renormalization
        //renormalize the surface charge density: sigma
        f_1 = f_0;
        
        //cout<<"test-0 "<<"f_0="<<f_0<<endl;
        if(f_0 < 0)
        {
            dphi_i = dphi_L;
            do
            {
                ratio_1   = ratio_0 - dphi_i;
                f_1=RenormEulerLagrange(ratio_1,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
                f_1 += sigma;
                dphi_i *= 0.5;
                //cout<<"test-1 "<<"f_1="<<f_1<<endl;
            }
            while(f_1 <= f_0);
        }
        else
        {
            dphi_i = dphi_L;
            do
            {
                ratio_1   = ratio_0 + dphi_i;
                f_1=RenormEulerLagrange(ratio_1,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
                f_1 += sigma;
                dphi_i *= 0.5;
                //cout<<"test-2 "<<"f_1="<<f_1<<endl;
            }
            while(f_1 >= f_0);
        }
        
        /////////////////------------------initialization---------///////////////////////
        ratio_2 = ratio_1 - f_1*(ratio_1-ratio_0)/(f_1-f_0);
        f_2=RenormEulerLagrange(ratio_2,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
        f_2 += sigma;
        
        //cout<<"test-3 "<<"f_2="<<f_2<<endl;
        ratio_iter= 0;
        do
        {
            ++ratio_iter;
            //ratio_t = ratio;
            
            temp_2 = (ratio_2-ratio_1)/(f_2-f_1);
            temp_1 = (ratio_1-ratio_0)/(f_1-f_0);
            
            ratio_temp = ratio_2 - temp_2*f_2 + (temp_2-temp_1)*f_2*f_1/(f_2-f_0);
            //cout<<"test-4 "<<"ratio_temp="<<ratio_temp<<endl;
            //cout<<"test-0-5   ratio_temp="<<ratio_temp<<" ratio_iter="<<ratio_iter<<" f_2="<<f_2<<" f_0="<<f_0<<endl;
            f_temp=RenormEulerLagrange(ratio_temp,rhoBM1,rhoBM2,D,MB1,MB2,Z,ff1,ff2,rho1,rho2,BesselZero,MODEL);
            //cout<<"test-5 "<<"f_temp="<<f_temp<<endl;
            f_temp += sigma;
            //cout<<"test-6 "<<"f_temp="<<f_temp<<endl;
            ratio_0 = ratio_1;
            ratio_1 = ratio_2;
            ratio_2 = ratio_temp;
            
            f_0 = f_1;
            f_1 = f_2;
            f_2 = f_temp;
            
            
            
            //cout<<"ratio_iter="<<ratio_iter<<" f="<<f<<" ratio_l="<<ratio_l<<" ratio_r="<<ratio_r<<endl;
            if(fabs(f_2) < ratio_err)
            {
                ratio = ratio_2;
                ratio_iter = 0;
                break;
            }
        }
        while(ratio_iter<ratio_MAX);
        
        if(ratio_iter >= ratio_MAX)
        {
            cerr<<"we cannot renormalize the charge neutrality for bulk"<<endl;
            exit(0);
            
        }
        
        //ratio = (sqrt(sigma*sigma+sigmaP*sigmaN)-sigma)/sigmaP;
        for(int k=LLIM; k<=ngrid_b; ++k)
        {
            for(short i=0; i<hspecies; ++i)
            {
                rho1[i][k]= rho2[i][k];
                rho1[i][ngrid-k]= rho1[i][k];
            }
        }
        
        deltaPhi = ratio;
        
        for(short i=0; i<hspecies; ++i)
        {
            delete [] rho2[i];
        }
        delete [] rho2;
        
    }//end renormalize
    //charge neutrality
    //////////////////////////////////////////////////////////////////////////////////////////////////
}
