//***********calculate the Gaussian Integration************//
#include "clibrary.h"
#include "gaussianintegration.h"
extern double dr;
extern double errTol;
//aa means the start grid of integrandF
//bb means the end grid of integrandF
double GaussianIntegrationZ0(double* integrandF,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr2_6;
    double   Gau2,gAA,gBB;
    double   gcoe1,gcoe2,temp1,temp2;
    
    final_result = 0;
    if(ib==ia) return final_result;
    gAA  = 0.5555555555556;
    gBB  = 0.8888888888889;
    Gau2 = 0.6;            //15/25
    
    gcoe1= dr*gAA*Gau2;
    gcoe2= dr*(2*gAA + gBB);

    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        temp1 = 0;
        temp2 = 0;
        for(int k=iaa; k<(ibb-1); k=(k+2))
        {
            temp1 += (integrandF[k] - 2*integrandF[k+1] + integrandF[k+2]);
            temp2 += integrandF[k+1];
        }
        final_result = (temp1*gcoe1 + temp2*gcoe2);
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



double GaussianIntegrationZ1(double* integrandF,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0,i_Ra,i_Rb;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr1_6,dr2_6;
    double   Gau2,gAA,gBB;
    double   gcoe1,gcoe2,temp1,temp2,temp3;
    
    final_result = 0;
    if(ib==ia) return final_result;
    
    gAA  = 0.5555555555556;
    gBB  = 0.8888888888889;
    Gau2 = 0.6;            //15/25
    dr2  = dr*dr;
    
    gcoe1= dr2*gAA*Gau2;
    gcoe2= dr2*(2*gAA + gBB);
    
    
    
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        for(int k=iaa; k<(ibb-1); k=(k+2))
        {
            temp1  = (integrandF[k] - 2*integrandF[k+1] + integrandF[k+2])*(k+1);
            temp2 += (temp1 + integrandF[k+2] - integrandF[k]);
            temp3 += (k+1)*integrandF[k+1];
        }
        final_result = (temp2*gcoe1 + temp3*gcoe2);
    }
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr1_1   = 1/dr;
        dr1_6   = dr1_1/6;
        dr2_6   = dr1_1*dr1_6;
        
        i_Ra = iaa0;
        i_Rb = ibb0;
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


double GaussianIntegrationZ2(double* integrandF,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2,dr3;
    double   coe1,coe2,coe3,coe4,dr1_1,i_D_Ra0,i_D_Rb0,i_D_Ra1,i_D_Rb1;
    double   Gau2,gAA,gBB;
    double   gcoe1,gcoe2,gcoe3,temp1,temp2,temp3,temp4,temp5,temp6;
    
    final_result = 0;
    if(ib==ia) return final_result;
    
    gAA  = 0.5555555555556;
    gBB  = 0.8888888888889;
    Gau2 = 0.6;            //15/25
    dr2  = dr*dr;
    dr3  = dr2*dr;
    
    gcoe1= dr3*gAA*Gau2;
    gcoe2= dr3*(2*gAA + gBB);
    gcoe3= 2*gcoe1;
    
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        temp4 = 0;
        temp5 = 0;
        temp6 = 0;
        for(int k=iaa; k<(ibb-1); k=(k+2))
        {
            temp1  = (k+1)*(k+1);
            temp2  = (integrandF[k] - 2*integrandF[k+1] + integrandF[k+2])*(temp1 + Gau2);
            temp3  = (k+1)*(integrandF[k+2] - integrandF[k]);
            temp4 += temp2;
            temp5 += temp1*integrandF[k+1];
            temp6 += (temp3 + integrandF[k+1]);
        }
        final_result = (temp4*gcoe1 + temp5*gcoe2 + temp6*gcoe3);
    }
    
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        i_D_Ra0 = iaa0*iaa0;
        i_D_Rb0 = ibb0*ibb0;
        i_D_Ra1 = (iaa0-1)*(iaa0-1);
        i_D_Rb1 = (ibb0+1)*(ibb0+1);
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


double GaussianIntegration(double* integrandF1,double* integrandF2,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr2_6;
    double   Gau2,gAA,gBB;
    double   gcoe1,gcoe2,temp1,temp2,temp3,temp4,temp5;
    
    final_result = 0;
    if(ib==ia) return final_result;
    gAA  = 0.5555555555556;
    gBB  = 0.8888888888889;
    Gau2 = 0.6;            //15/25
    
    gcoe1= dr*gAA*Gau2;
    gcoe2= dr*(2*gAA + gBB);
    
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        temp4 = 0;
        temp5 = 0;
        for(int k=iaa; k<(ibb-1); k=(k+2))
        {
            temp1  = integrandF1[k]*integrandF2[k];
            temp2  = integrandF1[k+1]*integrandF2[k+1];
            temp3  = integrandF1[k+2]*integrandF2[k+2];
            temp4 += (temp1 - 2*temp2 + temp3);
            temp5 += temp2;
        }
        final_result = (temp4*gcoe1 + temp5*gcoe2);
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



double GaussianIntegration(double* integrandF1,double* integrandF2,int orient0,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2;
    double   coe1,coe2,coe3,dr1_1,dr2_6;
    double   Gau2,gAA,gBB;
    double   gcoe1,gcoe2,temp1,temp2,temp3,temp4,temp5;
    
    final_result = 0;
    if(ib==ia)
    {
        final_result = (dr*0.5*integrandF1[orient0]*integrandF2[orient0]*4.0/3.0);
        return final_result;
    }
    
    gAA  = 0.5555555555556;
    gBB  = 0.8888888888889;
    Gau2 = 0.6;            //15/25
    
    gcoe1= dr*gAA*Gau2;
    gcoe2= dr*(2*gAA + gBB);
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        temp4 = 0;
        temp5 = 0;
        for(int k=iaa; k<(ibb-1); k=(k+2))
        {
            temp1  = integrandF1[k]*integrandF2[k];
            temp2  = integrandF1[k+1]*integrandF2[k+1];
            temp3  = integrandF1[k+2]*integrandF2[k+2];
            temp4 += (temp1 - 2*temp2 + temp3);
            temp5 += temp2;
        }
        final_result = (temp4*gcoe1 + temp5*gcoe2);
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



double GaussianIntegrationShell(double* integrandF,double* Psi_val,int aa,int bb,int ia,int ib)
{
    int      iaa,ibb,iaa0,ibb0,igg;//basic_n_mim,exp_v;
    double   final_result,final_result0,dr2,dr3;
    double   coe1,coe2,coe3,coe4,dr1_1,i_D_Ra0,i_D_Rb0,i_D_Ra1,i_D_Rb1;
    double   Gau1,Gau2,gAA,gBB,gg1,gg2,gg3;
    double   gcoe1,gcoe2,gcoe3,gcoe4,temp1,temp2,temp3,temp4;
        
    final_result = 0;
    if(ib==ia) return final_result;
    
    gAA  = 0.5555555555556;
    gBB  = 0.8888888888889;
    Gau1 = 0.7745966692415;
    Gau2 = 0.6;            //15/25
    dr2  = dr*dr;
    dr3  = dr2*dr;
    
    gcoe1= dr*gAA*Gau2*0.5;
    gcoe2= dr*gAA*Gau1*0.5;
    gcoe3= dr*gAA;
    gcoe4= dr*gBB;
    
    iaa  = ia;
    ibb  = ib;
    if((ib-ia+1)%2 == 0) ibb = ib -1;
    iaa0 = ibb;
    ibb0 = ib;
    final_result = 0;
    if((ibb-iaa) > 1)
    {
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        temp4 = 0;
        for(int k=iaa; k<(ibb-1); k=(k+2))
        {
            igg    = 3*k + 1;
            gg1    = Psi_val[igg];
            gg2    = Psi_val[igg+2];
            gg3    = Psi_val[igg+4];
            
            temp1 += (integrandF[k] - 2*integrandF[k+1] + integrandF[k+2])*(gg1 + gg3);
            temp2 += (integrandF[k+2] - integrandF[k])*(gg3 - gg1);
            temp3 += integrandF[k+1]*(gg1 + gg3);
            temp4 += integrandF[k+1]*gg2;
        }
        final_result = (temp1*gcoe1 + temp2*gcoe2 + temp3*gcoe3 + temp4*gcoe4);
    }
    
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        igg    = 3*iaa0;
        i_D_Ra0= Psi_val[igg];
        igg    = 3*ibb0;
        i_D_Rb0= Psi_val[igg];

        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF[iaa0]*i_D_Ra0 + integrandF[ibb0]*i_D_Rb0)*0.5*dr3;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            igg    = 3*(iaa0-1);
            i_D_Ra1= Psi_val[igg];
            igg    = 3*(ibb0+1);
            i_D_Rb1= Psi_val[igg];

            
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
                igg    = 3*(iaa0-1);
                i_D_Ra1= Psi_val[igg];

                coe1 = (integrandF[iaa0-1]*i_D_Ra1 - 2.0*integrandF[iaa0]*i_D_Ra0 +
                              integrandF[ibb0]*i_D_Rb0)/6;
                coe2 = (integrandF[ibb0]*i_D_Rb0 - integrandF[iaa0-1]*i_D_Ra1)*dr/4;
                coe3 = integrandF[iaa0]*i_D_Ra0*dr2;
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0-1]*(drr*i)*(drr*i-dr) - integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)
                //+ 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i));
            }
            else if(ibb0 < bb)
            {
                igg    = 3*(ibb0+1);
                i_D_Rb1= Psi_val[igg];

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





                     

                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     

                     
                     
                     
                     
                     
                    

