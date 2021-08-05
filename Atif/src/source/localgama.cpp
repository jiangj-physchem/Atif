//**calculating local screening parameter Gama(r) in funtionalized MSA**//
#include "clibrary.h"
#include "weighteddensity.h"
#include "localgama.h"
#include "constantnum.h"
#include "calculategama.h"


extern double dr;
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern short  nspecies;

///////////////////////////////////////////////////////////////////////////////////////////////////////
// In this subroutine, we consider two symmetric surface
void LocalGama(double gammab,int* LLI,int* ULI,float* D,float* Z,double** rho,double* gamar)
{
    double*  B;
    double** H;
    double** NI;

    double   temp1;
    short    hspecies;
    int      i_gama,k1,n_m_gama,n_gama;
    
    hspecies= nspecies - 1;
    
    B     = new double[hspecies]();


    H      = new double*[4]();
    NI     = new double*[hspecies]();
    for(short i=0; i<4; ++i)
    {
        H[i] = new double[hspecies]();
    }
    
    i_gama = round(0.5/(gammab*dr));
    n_m_gama = ngrid_m + i_gama;
    n_gama   = ngrid + i_gama*2;
    for(short i=0; i<hspecies; ++i)
    {
        NI[i]  = new double[4]();
        
        B[i] = 0.5*D[i] + i_gama*dr;

        temp1   = 0.5*D[i]/B[i];
        
        H[0][i] = 1.0;
        H[1][i] = temp1;
        H[2][i] = temp1*temp1;
        H[3][i] = H[2][i]*temp1;
    }
    
    //big loop
    for(int k=0; k<=n_m_gama; ++k)
    {
        k1 = k - i_gama;
        WeightedDensity(k1,LLI,ULI,B,H,rho,NI);
        
        gamar[k]=CalculateGama(Z,D,NI);
        
        gamar[n_gama-k] = gamar[k];
    }//big loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

}

