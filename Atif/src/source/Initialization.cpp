//*******************solve Euler Lagrange***************//

#include "Initialization.h"
#include "eulerlagrangeequation.h"
#include "eulerlagrangeequationDFT.h"
#include "getnumbersbyline.h"
//#include "poissonequation.h"
#include "fileoutput.h"

extern int ngrid; //the number of grids
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern short nspecies;
extern double errTol;


void Initialization(double gama,double f,double eta_t,int* LLI, int* ULI,float* D,double* BB,double* etar,float* Z,
                    double* rhoB,double* mu,double** Ext,double** ATT,double*** BesselZero,double*** Psi_IJ,double** rho,
                    double** pairEner,string* MODEL,string* MethG)
{
    int     iter;
    //int     thresh_copy,thresh_in,thresh_out;
    //short   copy_yes;
    double  mixCoe,maxItera,err,mix_f,temp,deltaPhi;//,lunit;
    //double  err_temp1,err_temp2,mix_f_temp1,mix_f_temp2;
    double*   Psi;
    double**  rho1;
    short   hspecies;
    //double**  rho_temp1;
    //double**  rho_temp2;
    
    mixCoe = 0.001;
    maxItera  = 1.0E6;
    //lunit     = 5.0E-10;
    hspecies = nspecies - 1;
    Psi  = new double[ngrid+1]();
    rho1 = new double*[nspecies]();
    //rho_temp1 = new double*[nspecies]();
    //rho_temp2 = new double*[nspecies]();
    for(short i=0; i<nspecies; ++i)
    {
        rho1[i] = new double[ngrid+1]();
        //rho_temp1[i] = new double[ngrid+1]();
        //rho_temp2[i] = new double[ngrid+1]();
    }
    
    for(int k=0; k<=ngrid; ++k)
    {
        for(short i=0; i<nspecies; ++i)
        {
            rho[i][k] = rhoB[i];
            if((k<LLI[i]) || (k>ULI[i])) rho[i][k] = 0;
        }
    }
    
    iter    = 0;
    mix_f   = mixCoe;
    deltaPhi= 0;
    err     = 0;
    do
    {
        if(MethG[0] == "scft")
        {
            EulerLagrange(deltaPhi,f,eta_t,LLI,ULI,D,Z,etar,Psi,pairEner,mu,rhoB,Ext,BesselZero,
                          rho,rho1,MODEL);
        }
        if(MethG[0] == "dft")
        {
            EulerLagrangeDFT(gama,deltaPhi,f,eta_t,LLI,ULI,D,BB,Z,etar,Psi,pairEner,mu,rhoB,Ext,ATT,
                             BesselZero,Psi_IJ,rho,rho1,MODEL);
        }
        
        
        err = 0;
        for(short i=0; i<hspecies; ++i)
        {
            if(rhoB[i] > 1E-10)
            {
                for(int k=LLI[i]; k<=ngrid_m; ++k)
                {
                    temp = fabs(rho1[i][k]-rho[i][k])/rhoB[i];
                    if(err < temp) err = temp;
                }
            }
        }
        
        for(short i=0; i<nspecies; ++i)
        {
            for(int k=LLI[i]; k<=ngrid_m; ++k)
            {
                rho[i][k] = (1-mix_f)*rho[i][k] + mix_f*rho1[i][k];
                rho[i][ngrid-k] = rho[i][k];
            }
        }
        
        ++iter;
        
    }while((err >= errTol) && (iter <= maxItera));
    
    //make sure the iteration convergence
    for(short i=0; i<nspecies; ++i)
    {
        delete [] rho1[i];
        //delete [] rho_temp1[i];
        //delete [] rho_temp2[i];
    }
    delete [] rho1;
    //delete [] rho_temp1;
    //delete [] rho_temp2;
    delete [] Psi;
    
}



void Initialization(double& sigma_t,double& dcharge,int& cstep,int& iter,double** rho,ifstream& inFile)
{
    int     ngrid_old,ngrid_old2,m_old2,kL,kR,ngridk;
    double* FirstLine;
    double* InputRho;
    
    FirstLine = new double[5]();
    InputRho  = new double[nspecies]();
    GetNumbersByLine(inFile,FirstLine);
    sigma_t      = FirstLine[0];
    dcharge      = FirstLine[1];
    cstep        = (int) FirstLine[2];
    iter         = FirstLine[3];
    ngrid_old    = FirstLine[4];
    //inFile>>psi_delta>>C_psi>>iter>>ngrid_old>>Ener_Save[0];
    
    if(ngrid == ngrid_old)
    {
        for(int k=0; k<=ngrid_m; ++k)
        {
            //inFile>>rho[0][k]>>rho[1][k]>>rho[2][k]>>rho[3][k]>>rho[4][k]>>rho[5][k]>>rho[6][k]>>rho[7][k]>>rho[8][k];
            GetNumbersByLine(inFile,InputRho);
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][k] = InputRho[i];
                rho[i][ngrid-k] = rho[i][k];
            }
        }
        
    }
    else
    {//222222
        ngrid_old2 = ngrid_old/2;
        m_old2     = ngrid_old%2;
        for(int k=0; k<=ngrid_old2; ++k)
        {
            //inFile>>rho[0][k]>>rho[1][k]>>rho[2][k]>>rho[3][k]>>rho[4][k]>>rho[5][k]>>rho[6][k]>>rho[7][k]>>rho[8][k];
            GetNumbersByLine(inFile,InputRho);
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][k] = InputRho[i];
                rho[i][ngrid-k] = rho[i][k];
            }
        }

        
        kL = ngrid_old2 + 1;
        kR = ngrid-ngrid_old2-m_old2;
        if(m_old2 == 0)
        {
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][kR] = rho[i][kL-1];
            }
            kR= kR - 1;
        }
        
        if(kR >= kL)//111111
        {
            ngridk = (kL+kR)/2;
            if((kR-kL+1)%2 == 0)
            {
                for(int k=kL; k<=ngridk; ++k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k-1] - rho[i][k-2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                }
                
                for(int k=kR; k>=(ngridk+1); --k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k+1] - rho[i][k+2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                    
                }
                
            }
            else
            {
                for(int k=kL; k<ngridk; ++k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k-1] - rho[i][k-2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                }
                
                for(int k=kR; k>ngridk; --k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k+1] - rho[i][k+2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                    
                }
                
                for(short i=0; i<nspecies; ++i)
                {
                    rho[i][ngridk] = rho[i][ngridk+1] + rho[i][ngridk-1] -
                    0.5*(rho[i][ngridk+2]+rho[i][ngridk-2]);
                    if(rho[i][ngridk] < 0) rho[i][ngridk] = 0.0;
                }
                
            }
        }//111111
    }//2222222
    
    delete [] FirstLine;
    delete [] InputRho;
}


void Initialization(int* LLI, int* ULI,double* rhoB, double** rho)
{
    for(int k=0; k<=ngrid; ++k)
    {
        for(short i=0; i<nspecies; ++i)
        {
            rho[i][k] = rhoB[i];
            if((k<LLI[i]) || (k>ULI[i])) rho[i][k] = 0;
        }
    }
}


void Initialization(double& sigma_t,double& dcharge,int& cstep,int& iter,int* LLI, int* ULI,double* rhoB,
                    double** rho,ifstream& inFile)
{
    
    double* FirstLine;
    
    FirstLine = new double[5]();
    GetNumbersByLine(inFile,FirstLine);
    sigma_t      = FirstLine[0];
    dcharge      = FirstLine[1];
    cstep        = (int) FirstLine[2];
    iter         = FirstLine[3];
    
    
    for(int k=0; k<=ngrid; ++k)
    {
        for(short i=0; i<nspecies; ++i)
        {
            rho[i][k] = rhoB[i];
            if((k<LLI[i]) || (k>ULI[i])) rho[i][k] = 0;
        }
    }
    
    
    delete [] FirstLine;
}

/*
void Initialization(double& sigma_t,double& dcharge,int& cstep,int& iter,double** rho,ifstream& inFile)
{
    int     ngrid_old,ngrid_old2,m_old2,kL,kR,ngridk;
    double* FirstLine;
    double* InputRho;
    
    FirstLine = new double[5]();
    InputRho  = new double[nspecies]();
    GetNumbersByLine(inFile,FirstLine);
    sigma_t      = FirstLine[0];
    dcharge      = FirstLine[1];
    cstep        = (int) FirstLine[2];
    iter         = FirstLine[3];
    ngrid_old    = FirstLine[4];
    //inFile>>psi_delta>>C_psi>>iter>>ngrid_old>>Ener_Save[0];
    
    if(ngrid == ngrid_old)
    {
        for(int k=0; k<=ngrid; ++k)
        {
            //inFile>>rho[0][k]>>rho[1][k]>>rho[2][k]>>rho[3][k]>>rho[4][k]>>rho[5][k]>>rho[6][k]>>rho[7][k]>>rho[8][k];
            GetNumbersByLine(inFile,InputRho);
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][k] = InputRho[i];
            }
        }
        
    }
    else
    {//222222
        ngrid_old2 = ngrid_old/2;
        m_old2     = ngrid_old%2;
        for(int k=0; k<=ngrid_old2; ++k)
        {
            //inFile>>rho[0][k]>>rho[1][k]>>rho[2][k]>>rho[3][k]>>rho[4][k]>>rho[5][k]>>rho[6][k]>>rho[7][k]>>rho[8][k];
            GetNumbersByLine(inFile,InputRho);
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][k] = InputRho[i];
            }
        }
        for(int k=(ngrid-ngrid_old2+1-m_old2); k<=ngrid; ++k)
        {
            //inFile>>rho[0][k]>>rho[1][k]>>rho[2][k]>>rho[3][k]>>rho[4][k]>>rho[5][k]>>rho[6][k]>>rho[7][k]>>rho[8][k];
            GetNumbersByLine(inFile,InputRho);
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][k] = InputRho[i];
            }
        }
        
        kL = ngrid_old2 + 1;
        kR = ngrid-ngrid_old2-m_old2;
        if(m_old2 == 0)
        {
            for(short i=0; i<nspecies; ++i)
            {
                rho[i][kR] = rho[i][kL-1];
            }
            kR= kR - 1;
        }
        
        if(kR >= kL)//111111
        {
            ngridk = (kL+kR)/2;
            if((kR-kL+1)%2 == 0)
            {
                for(int k=kL; k<=ngridk; ++k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k-1] - rho[i][k-2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                }
                
                for(int k=kR; k>=(ngridk+1); --k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k+1] - rho[i][k+2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                    
                }
                
            }
            else
            {
                for(int k=kL; k<ngridk; ++k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k-1] - rho[i][k-2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                }
                
                for(int k=kR; k>ngridk; --k)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        rho[i][k] = 2.0*rho[i][k+1] - rho[i][k+2];
                        if(rho[i][k] < 0) rho[i][k] = 0.0;
                    }
                    
                }
                
                for(short i=0; i<nspecies; ++i)
                {
                    rho[i][ngridk] = rho[i][ngridk+1] + rho[i][ngridk-1] -
                    0.5*(rho[i][ngridk+2]+rho[i][ngridk-2]);
                    if(rho[i][ngridk] < 0) rho[i][ngridk] = 0.0;
                }
                
            }
        }//111111
    }//2222222
    
    delete [] FirstLine;
    delete [] InputRho;
}
 */


