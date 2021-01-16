//***********calculating external potential****************//
#include "clibrary.h"
#include "constantnum.h"
#include "externalpotential.h"
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern int LLIM; //the minimum lower limit of intergral
extern short nspecies;
extern double dr;

void ExternalPotential(float* uWall,int** depthLR,double** Ext)
{
    short hspecies;
    
    hspecies= nspecies - 1;
    for(int k=LLIM; k<=ngrid_m; ++k)
    {
        for(short i=0; i<hspecies; ++i)
        {
            Ext[i][k] = 0;
            if(k<depthLR[0][i] || k>depthLR[2][i])
            {
                Ext[i][k] += 1E20;
            }
            
            if(k>=depthLR[0][i] && k<=depthLR[1][i])
            {
                Ext[i][k] += uWall[i];
            }
            
            if(k>=depthLR[3][i] && k<=depthLR[2][i])
            {
                Ext[i][k] += uWall[i];
            }
        }
    }
}




void ExternalPotential(double alpha,float* D,float* uWall,int** depthLR,double** Ext)
{
    short hspecies;
    //int depthL1,depthL2,depthR1,depthR2;
    double DR3L,DR3R,coe,RL,RR,D3,R_star,E_cut,R_cut,R_min,E_min,R_minr,E_min0;
    
    hspecies= nspecies - 1;
    R_min  = pow(0.4,1.0/6.0);
    R_minr = 1.0/R_min;
    //E_min  = 2.0*Pi*(2.0*pow(R_minr,9.0)/45.0-pow(R_minr,3.0)/3.0);
    E_min  = 2.0*pow(R_minr,9.0)/45.0-pow(R_minr,3.0)/3.0;
    
    for(short i=0; i<hspecies; ++i)
    {
        R_star = R_min*D[i] - depthLR[0][i]*dr;
        R_cut  = D[i]/(R_star+alpha*D[i]);
        //coe    = 2.0*Pi*uWall[i];
        
        E_cut  = 2.0*pow(R_cut,9.0)/45.0-pow(R_cut,3.0)/3.0;
        E_min0 = E_min - E_cut;
        coe    = uWall[i]/E_min0;
        
        
        D3     = D[i]*D[i]*D[i];
        for(int k=LLIM; k<=ngrid_m; ++k)
        {
            Ext[i][k] = 0;
            if(k<depthLR[0][i] || k>depthLR[2][i]) Ext[i][k] += 2E10;
            
            if(k>=depthLR[0][i] && k<=depthLR[1][i])
            {
                RL  = dr*k + R_star;
                DR3L= D3/(RL*RL*RL);
                
                Ext[i][k] += coe*(2.0*DR3L*DR3L*DR3L/45.0-DR3L/3.0-E_cut);
            }
            
            if(k>=depthLR[3][i] && k<=depthLR[2][i])
            {
                RR  = dr*(ngrid-k) + R_star;
                DR3R= D3/(RR*RR*RR);
                
                Ext[i][k] += coe*(2.0*DR3R*DR3R*DR3R/45.0-DR3R/3.0-E_cut);
            }
            
            
        }
        
        
        
    }
    

    
    
    
}
