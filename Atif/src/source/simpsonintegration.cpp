//***********calculate the Simpson Integration************//
#include "clibrary.h"
#include "simpsonintegration.h"
extern double dr;
extern double errTol;
//aa means the start grid of integrandF
//bb means the end grid of integrandF
double SimpsonIntegration(double* integrandF,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr2_6;
    
    final_result = 0;
    if(ib==ia) return final_result;


    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        final_result += (integrandF[iaa] + integrandF[ibb]);
        final_result += (integrandF[ibb-1]*4.0);
        for(int k=(iaa+1); k<(ibb-2); k=(k+2))
        {
            final_result += (integrandF[k]*4.0);
            final_result += (integrandF[k+1]*2.0);
        }
        final_result = final_result*dr/3.0;
    }
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr2     = dr*dr;
        dr1_1   = 1/dr;
        dr2_6   = dr1_1*dr1_1/6;
        
        
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF[iaa0] + integrandF[ibb0])*0.5*dr;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            coe1 = (integrandF[ibb0+1] + integrandF[iaa0]*3.0 - integrandF[ibb0]*3.0 - integrandF[iaa0-1])*dr1_1*dr2_6/4.0;
            coe2 = (integrandF[iaa0-1] + integrandF[ibb0] - 2.0*integrandF[iaa0])*dr2_6;
            coe3 = (6.0*integrandF[ibb0] - 3.0*integrandF[iaa0] - 2.0*integrandF[iaa0-1] - integrandF[ibb0+1])*dr1_1/12;
            coe4 = integrandF[iaa0];
            
            final_result0 = dr2*dr2*coe1 + dr2*dr*coe2 + dr2*coe3 + dr*coe4;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
            
            if(iaa0 > aa)
            {
                coe1 = dr2_6*(integrandF[iaa0-1] - 2.0*integrandF[iaa0] + integrandF[ibb0]);
                coe2 = (integrandF[ibb0] - integrandF[iaa0-1])*dr1_1/4;
                coe3 = integrandF[iaa0];
            }
            else if(ibb0 < bb)
            {
                coe1 = dr2_6*(integrandF[iaa0] - 2.0*integrandF[ibb0] + integrandF[ibb0+1]);
                coe2 = (4.0*integrandF[ibb0] - integrandF[ibb0+1] - 3.0*integrandF[iaa0])*dr1_1/4;
                coe3 = integrandF[iaa0];
            }
            
            final_result0 = dr2*dr*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}





double SimpsonIntegration(double* integrandF,int aa,int bb,int ia,int ib,int i_R)
{
    int      iaa,ibb,iaa0,ibb0,i_Ra,i_Rb,i_Rk;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr1_6,dr2_6;
    
    final_result = 0;
    if(ib==ia) return final_result;
    
    dr2  = dr*dr;
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        i_Rk = ibb-i_R;
        final_result += (integrandF[iaa]*(iaa-i_R) + integrandF[ibb]*i_Rk);
        final_result += (integrandF[ibb-1]*(i_Rk-1)*4.0);
        for(int k=(iaa+1); k<(ibb-2); k=(k+2))
        {
            i_Rk = k - i_R;
            final_result += (integrandF[k]*i_Rk*4.0);
            final_result += (integrandF[k+1]*(i_Rk+1)*2.0);
        }
        final_result = final_result*dr2/3.0;
    }
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr1_1   = 1/dr;
        dr1_6   = dr1_1/6;
        dr2_6   = dr1_1*dr1_6;
        
        i_Ra = iaa0 - i_R;
        i_Rb = ibb0 - i_R;
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF[iaa0]*i_Ra + integrandF[ibb0]*i_Rb)*0.5*dr2;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            coe1 = (integrandF[ibb0+1]*(i_Rb+1) + integrandF[iaa0]*i_Ra*3.0 -
                    integrandF[ibb0]*i_Rb*3.0 - integrandF[iaa0-1]*(i_Ra-1))*dr2_6/4;
            coe2 = (integrandF[iaa0-1]*(i_Ra-1) + integrandF[ibb0]*i_Rb -
                    2.0*integrandF[iaa0]*i_Ra)*dr1_6;
            coe3 = (6.0*integrandF[ibb0]*i_Rb - 3.0*integrandF[iaa0]*i_Ra -
                    2.0*integrandF[iaa0-1]*(i_Ra-1) - integrandF[ibb0+1]*(i_Rb+1))/12;
            coe4 = integrandF[iaa0]*i_Ra*dr;
            //FF_IN[i] = dr3*(integrandF[ibb0+1]*(drr*i+dr)*(drr*i)*(drr*i-dr)/6.0 + 0.5*integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)*(drr*i-2*dr)
            //- 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i)*(drr*i-2*dr) - integrandF[iaa0-1]*(drr*i)*(drr*i-dr)*(drr*i-2*dr)/6.0);
            
            final_result0 = dr2*dr2*coe1 + dr2*dr*coe2 + dr2*coe3 + dr*coe4;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
            if(iaa0 > aa)
            {
                coe1 = dr1_6*(integrandF[iaa0-1]*(i_Ra-1) - 2.0*integrandF[iaa0]*i_Ra +
                              integrandF[ibb0]*i_Rb);
                coe2 = (integrandF[ibb0]*i_Rb - integrandF[iaa0-1]*(i_Ra-1))/4;
                coe3 = integrandF[iaa0]*i_Ra*dr;
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0-1]*(drr*i)*(drr*i-dr) - integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)
                //+ 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i));
            }
            else if(ibb0 < bb)
            {
                coe1 = dr1_6*(integrandF[iaa0]*i_Ra - 2.0*integrandF[ibb0]*i_Rb +
                              integrandF[ibb0+1]*(i_Rb+1));
                coe2 = (4.0*integrandF[ibb0]*i_Rb - integrandF[ibb0+1]*(i_Rb+1) -
                        3.0*integrandF[iaa0]*i_Ra)/4;
                coe3 = integrandF[iaa0]*i_Ra*dr;
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0]*(drr*i-dr)*(drr*i-2*dr) - integrandF[ibb0]*(drr*i)*(drr*i-2*dr)
                //+ 0.5*integrandF[ibb0+1]*(drr*i)*(drr*i-dr));
            }
            
            final_result0 = dr2*dr*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}



double SimpsonIntegration(double* integrandF,int aa,int bb,int ia,int ib,double D2_j,int i_R)
{
    int      iaa,ibb,iaa0,ibb0,i_Rk;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2,dr3,DD_2;
    double   coe1,coe2,coe3,coe4,dr1_1,i_D_Ra0,i_D_Rb0,i_D_Ra1,i_D_Rb1;
    
    dr2   = dr*dr;
    dr3   = dr2*dr;
    DD_2  = D2_j/dr2;
    
    final_result = 0;
    if(ib==ia) return final_result;

    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        i_Rk = ibb-i_R;
        final_result += (integrandF[iaa]*(DD_2-(iaa-i_R)*(iaa-i_R)) +
                         integrandF[ibb]*(DD_2-i_Rk*i_Rk));
        final_result += (integrandF[ibb-1]*(DD_2-(i_Rk-1)*(i_Rk-1))*4.0);
        for(int k=(iaa+1); k<(ibb-2); k=(k+2))
        {
            i_Rk = k - i_R;
            final_result += (integrandF[k]*(DD_2-i_Rk*i_Rk)*4.0);
            final_result += (integrandF[k+1]*(DD_2-(i_Rk+1)*(i_Rk+1))*2.0);
        }
        final_result = final_result*dr3/3.0;
    }
    
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        i_D_Ra0 = DD_2 - (iaa0 - i_R)*(iaa0 - i_R);
        i_D_Rb0 = DD_2 - (ibb0 - i_R)*(ibb0 - i_R);
        i_D_Ra1 = DD_2 - (iaa0-i_R-1)*(iaa0-i_R-1);
        i_D_Rb1 = DD_2 - (ibb0-i_R+1)*(ibb0-i_R+1);
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF[iaa0]*i_D_Ra0 + integrandF[ibb0]*i_D_Rb0)*0.5*dr3;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            dr1_1   = 1/dr;
            coe1 = (integrandF[ibb0+1]*i_D_Rb1 + integrandF[iaa0]*i_D_Ra0*3.0 -
                    integrandF[ibb0]*i_D_Rb0*3.0 - integrandF[iaa0-1]*i_D_Ra1)*dr1_1/24;
            coe2 = (integrandF[iaa0-1]*i_D_Ra1 + integrandF[ibb0]*i_D_Rb0 -
                    2.0*integrandF[iaa0]*i_D_Ra0)/6;
            coe3 = (6.0*integrandF[ibb0]*i_D_Rb0 - 3.0*integrandF[iaa0]*i_D_Ra0 -
                    2.0*integrandF[iaa0-1]*i_D_Ra1 - integrandF[ibb0+1]*i_D_Rb1)*dr/12;
            coe4 = integrandF[iaa0]*i_D_Ra0*dr2;
            //FF_IN[i] = dr3*(integrandF[ibb0+1]*(drr*i+dr)*(drr*i)*(drr*i-dr)/6.0 + 0.5*integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)*(drr*i-2*dr)
            //- 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i)*(drr*i-2*dr) - integrandF[iaa0-1]*(drr*i)*(drr*i-dr)*(drr*i-2*dr)/6.0);
            
            final_result0 = dr2*dr2*coe1 + dr3*coe2 + dr2*coe3 + dr*coe4;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
            if(iaa0 > aa)
            {
                coe1 = (integrandF[iaa0-1]*i_D_Ra1 - 2.0*integrandF[iaa0]*i_D_Ra0 +
                              integrandF[ibb0]*i_D_Rb0)/6;
                coe2 = (integrandF[ibb0]*i_D_Rb0 - integrandF[iaa0-1]*i_D_Ra1)*dr/4;
                coe3 = integrandF[iaa0]*i_D_Ra0*dr2;
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0-1]*(drr*i)*(drr*i-dr) - integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)
                //+ 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i));
            }
            else if(ibb0 < bb)
            {
                coe1 = (integrandF[iaa0]*i_D_Ra0 - 2.0*integrandF[ibb0]*i_D_Rb0 +
                              integrandF[ibb0+1]*i_D_Rb1)/6;
                coe2 = (4.0*integrandF[ibb0]*i_D_Rb0 - integrandF[ibb0+1]*i_D_Rb1 -
                        3.0*integrandF[iaa0]*i_D_Ra0)*dr/4;
                coe3 = integrandF[iaa0]*i_D_Ra0*dr2;
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0]*(drr*i-dr)*(drr*i-2*dr) - integrandF[ibb0]*(drr*i)*(drr*i-2*dr)
                //+ 0.5*integrandF[ibb0+1]*(drr*i)*(drr*i-dr));
            }
            
            final_result0 = dr3*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}


double SimpsonIntegration(int aa,int bb,int ia,int ib,double D2_j)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2,dr3,DD_2;
    double   coe1,coe2,coe3,i_D_Ra0,i_D_Rb0,i_D_Ra1;
    
    dr2   = dr*dr;
    dr3   = dr2*dr;
    DD_2  = D2_j/dr2;
    
    final_result = 0;
    if(ib==ia) return final_result;
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        final_result += ((DD_2-iaa*iaa) + (DD_2-ibb*ibb));
        final_result += ((DD_2-(ibb-1)*(ibb-1))*4.0);
        for(int k=(iaa+1); k<(ibb-2); k=(k+2))
        {
            final_result += (DD_2-k*k)*4.0;
            final_result += (DD_2-(k+1)*(k+1))*2.0;
        }
        final_result = final_result*dr3/3.0;
    }
    
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        i_D_Ra0 = DD_2 - iaa0*iaa0;
        i_D_Rb0 = DD_2 - ibb0*ibb0;
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (i_D_Ra0 + i_D_Rb0)*0.5*dr3;
        }
        else
        {
            i_D_Ra1 = DD_2 - (iaa0-1)*(iaa0-1);
            
            coe1 = (i_D_Ra1 - 2.0*i_D_Ra0 + i_D_Rb0)/6;
            coe2 = (i_D_Rb0 - i_D_Ra1)*dr/4;
            coe3 = i_D_Ra0*dr2;
            
            final_result0 = dr3*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}


double SimpsonIntegration(double* integrandF1,double* integrandF2,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr2_6;
    
    final_result = 0;
    if(ib==ia) return final_result;
    
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        final_result += (integrandF1[iaa]*integrandF2[iaa] + integrandF1[ibb]*integrandF2[ibb]);
        final_result += (integrandF1[ibb-1]*integrandF2[ibb-1]*4.0);
        for(int k=(iaa+1); k<(ibb-2); k=(k+2))
        {
            final_result += (integrandF1[k]*integrandF2[k]*4.0);
            final_result += (integrandF1[k+1]*integrandF2[k+1]*2.0);
        }
        final_result = final_result*dr/3.0;
    }
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr2     = dr*dr;
        dr1_1   = 1/dr;
        dr2_6   = dr1_1*dr1_1/6;
        
        
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF1[iaa0]*integrandF2[iaa0] + integrandF1[ibb0]*integrandF2[ibb0])*0.5*dr;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            coe1 = (integrandF1[ibb0+1]*integrandF2[ibb0+1] + integrandF1[iaa0]*integrandF2[iaa0]*3.0 -
                    integrandF1[ibb0]*integrandF2[ibb0]*3.0 - integrandF1[iaa0-1]*integrandF2[iaa0-1])*dr1_1*dr2_6/4.0;
            coe2 = (integrandF1[iaa0-1]*integrandF2[iaa0-1] + integrandF1[ibb0]*integrandF2[ibb0] -
                    2.0*integrandF1[iaa0]*integrandF2[iaa0])*dr2_6;
            coe3 = (6.0*integrandF1[ibb0]*integrandF2[ibb0] - 3.0*integrandF1[iaa0]*integrandF2[iaa0] -
                    2.0*integrandF1[iaa0-1]*integrandF2[iaa0-1]  - integrandF1[ibb0+1]*integrandF2[ibb0+1])*dr1_1/12;
            coe4 = integrandF1[iaa0]*integrandF2[iaa0];
            
            final_result0 = dr2*dr2*coe1 + dr2*dr*coe2 + dr2*coe3 + dr*coe4;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
            
            if(iaa0 > aa)
            {
                coe1 = dr2_6*(integrandF1[iaa0-1]*integrandF2[iaa0-1] - 2.0*integrandF1[iaa0]*integrandF2[iaa0] +
                              integrandF1[ibb0]*integrandF2[ibb0]);
                coe2 = (integrandF1[ibb0]*integrandF2[ibb0] - integrandF1[iaa0-1]*integrandF2[iaa0-1])*dr1_1/4;
                coe3 = integrandF1[iaa0]*integrandF2[iaa0];
            }
            else if(ibb0 < bb)
            {
                coe1 = dr2_6*(integrandF1[iaa0]*integrandF2[iaa0] - 2.0*integrandF1[ibb0]*integrandF2[ibb0] +
                              integrandF1[ibb0+1]*integrandF2[ibb0+1]);
                coe2 = (4.0*integrandF1[ibb0]*integrandF2[ibb0] - integrandF1[ibb0+1]*integrandF2[ibb0+1] -
                        3.0*integrandF1[iaa0]*integrandF2[iaa0])*dr1_1/4;
                coe3 = integrandF1[iaa0]*integrandF2[iaa0];
            }
            
            final_result0 = dr2*dr*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}



double SimpsonIntegration(double* integrandF1,double* integrandF2,int orient0,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,dr1_1,dr2_6;
    
    final_result = 0;
    if(ib==ia)
    {
        final_result = (dr*0.5*integrandF1[orient0]*integrandF2[orient0]*4.0/3.0);
        return final_result;
    }
    
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        final_result += (integrandF1[iaa]*integrandF2[iaa] + integrandF1[ibb]*integrandF2[ibb]);
        final_result += (integrandF1[ibb-1]*integrandF2[ibb-1]*4.0);
        for(int k=(iaa+1); k<(ibb-2); k=(k+2))
        {
            final_result += (integrandF1[k]*integrandF2[k]*4.0);
            final_result += (integrandF1[k+1]*integrandF2[k+1]*2.0);
        }
        final_result = final_result*dr/3.0;
    }
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        if((iaa0==ia) && (ibb0==ib))
        {
            final_result0 = (integrandF1[iaa0]*integrandF2[iaa0] + integrandF1[ibb0]*integrandF2[ibb0])*0.5*dr;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
        
            dr2     = dr*dr;
            dr1_1   = 1/dr;
            dr2_6   = dr1_1*dr1_1/6;
            
            if(iaa0 > ia)
            {
                coe1 = dr2_6*(integrandF1[iaa0-1]*integrandF2[iaa0-1] - 2.0*integrandF1[iaa0]*integrandF2[iaa0] +
                              integrandF1[ibb0]*integrandF2[ibb0]);
                coe2 = (integrandF1[ibb0]*integrandF2[ibb0] - integrandF1[iaa0-1]*integrandF2[iaa0-1])*dr1_1/4;
                coe3 = integrandF1[iaa0]*integrandF2[iaa0];
            }
            else if(ibb0 < ib)
            {
                coe1 = dr2_6*(integrandF1[iaa0]*integrandF2[iaa0] - 2.0*integrandF1[ibb0]*integrandF2[ibb0] +
                              integrandF1[ibb0+1]*integrandF2[ibb0+1]);
                coe2 = (4.0*integrandF1[ibb0]*integrandF2[ibb0] - integrandF1[ibb0+1]*integrandF2[ibb0+1] -
                        3.0*integrandF1[iaa0]*integrandF2[iaa0])*dr1_1/4;
                coe3 = integrandF1[iaa0]*integrandF2[iaa0];
            }
            
            final_result0 = dr2*dr*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}







                     

                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     

                     
                     
                     
                     
                     
                    

