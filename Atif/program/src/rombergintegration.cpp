//***********calculate the Romberg Integration************//
#include "clibrary.h"
#include "rombergintegration.h"
extern double dr;
extern double errTol;
//aa means the start grid of integrandF
//bb means the end grid of integrandF
double RombergIntegration(double* integrandF,int aa,int bb,int ia,int ib)
{
    int      basic_n,iaa,ibb,iaa0,ibb0,pow2_deep;//basic_n_mim,exp_v;
    short    max_num_Rom,num_deep;
    double   basic_dr,final_result,final_result0;
    double   coe1,coe2,coe3,coe4,dr1_1,dr1_6,dr2_1,dr2_2,dr3_6;
    
    max_num_Rom = 8; //maximum extrapolation times of Romberg
    //basic_n_mim = 17; //minimum number of grids
    
    //exp_v    = (ib-ia)/(basic_n_mim-1);
    final_result = 0;
    if(ib==ia) return final_result;
    num_deep = (short) log2(ib-ia);
    num_deep = num_deep + 1;
    if(num_deep > max_num_Rom) num_deep = max_num_Rom;
    pow2_deep= pow(2,(num_deep-1));
    basic_n  = (ib-ia)/pow2_deep + 1;
    basic_dr = dr*pow2_deep;
    
    iaa  = ia;
    ibb  = ia + (basic_n-1)*pow2_deep;
    iaa0 = ibb;
    ibb0 = ib;
    
    if((ibb-iaa) > 1)
    {
        final_result = RombergIntegration(integrandF,basic_dr,num_deep,basic_n,iaa,ibb);
    }
    else
    {
        iaa0 = iaa;
        ibb0 = ibb;
    }
    
    while((ibb0-iaa0) > 1)
    {
        
        //exp_v    = ibb0-iaa0;
        num_deep = (short) log2(ibb0-iaa0);
        num_deep = num_deep + 1;
        if(num_deep > max_num_Rom) num_deep = max_num_Rom;
        pow2_deep= pow(2,(num_deep-1));
        basic_n  = (ibb0-iaa0)/pow2_deep + 1;
        basic_dr = dr*pow2_deep;
        
        iaa = iaa0;
        ibb = iaa0 + (basic_n-1)*pow2_deep;
        iaa0= ibb;
        
        final_result0 = RombergIntegration(integrandF,basic_dr,num_deep,basic_n,iaa,ibb);
        final_result += final_result0;
        
    }
    
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr1_1   = 1/dr;
        dr1_6   = dr1_1/6;
        dr2_1   = 0.5*dr1_1;
        dr2_2   = dr1_1*dr2_1;
        dr3_6   = dr1_1*dr2_2/3.0;
        
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF[iaa0] + integrandF[ibb0])*0.5*dr;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            coe1 = (integrandF[ibb0+1] + integrandF[iaa0]*3.0 - integrandF[ibb0]*3.0 - integrandF[iaa0-1])*dr3_6/4;
            coe2 = (integrandF[iaa0-1] + integrandF[ibb0] - 2.0*integrandF[iaa0])*dr2_2/3;
            coe3 = (6.0*integrandF[ibb0] - 3.0*integrandF[iaa0] - 2.0*integrandF[iaa0-1] - integrandF[ibb0+1])*dr1_6/2;
            coe4 = integrandF[iaa0];
            //FF_IN[i] = dr3*(integrandF[ibb0+1]*(drr*i+dr)*(drr*i)*(drr*i-dr)/6.0 + 0.5*integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)*(drr*i-2*dr)
            //- 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i)*(drr*i-2*dr) - integrandF[iaa0-1]*(drr*i)*(drr*i-dr)*(drr*i-2*dr)/6.0);
            
            final_result0 = pow(dr,4)*coe1 + pow(dr,3)*coe2 + dr*dr*coe3 + dr*coe4;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
            if(iaa0 > aa)
            {
                coe1 = dr2_2*(integrandF[iaa0-1] - 2.0*integrandF[iaa0] + integrandF[ibb0])/3;
                coe2 = (integrandF[ibb0] - integrandF[iaa0-1])*dr2_1/2;
                coe3 = integrandF[iaa0];
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0-1]*(drr*i)*(drr*i-dr) - integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)
                //+ 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i));
            }
            else if(ibb0 < bb)
            {
                coe1 = dr2_2*(integrandF[iaa0] - 2.0*integrandF[ibb0] + integrandF[ibb0+1])/3;
                coe2 = (4.0*integrandF[ibb0] - integrandF[ibb0+1] - 3.0*integrandF[iaa0])*dr2_1/2;
                coe3 = integrandF[iaa0];
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0]*(drr*i-dr)*(drr*i-2*dr) - integrandF[ibb0]*(drr*i)*(drr*i-2*dr)
                //+ 0.5*integrandF[ibb0+1]*(drr*i)*(drr*i-dr));
            }
            
            final_result0 = pow(dr,3)*coe1 + dr*dr*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}








double RombergIntegration(double* integrandF1,double* integrandF2,int aa,int bb,int ia,int ib)
{
    int      basic_n,iaa,ibb,iaa0,ibb0,pow2_deep;//basic_n_mim,exp_v;
    short    max_num_Rom,num_deep;
    double   basic_dr,final_result,final_result0;
    double   coe1,coe2,coe3,coe4,dr1_1,dr1_6,dr2_1,dr2_2,dr3_6;
    
    max_num_Rom = 8; //maximum extrapolation times of Romberg
    //basic_n_mim = 17; //minimum number of grids
    
    //exp_v    = (ib-ia)/(basic_n_mim-1);
    final_result = 0;
    if(ib==ia) return final_result;
    num_deep = (short) log2(ib-ia);
    num_deep = num_deep + 1;
    if(num_deep > max_num_Rom) num_deep = max_num_Rom;
    pow2_deep= pow(2,(num_deep-1));
    basic_n  = (ib-ia)/pow2_deep + 1;
    basic_dr = dr*pow2_deep;
    
    iaa  = ia;
    ibb  = ia + (basic_n-1)*pow2_deep;
    iaa0 = ibb;
    ibb0 = ib;
    
    if((ibb-iaa) > 1)
    {
        final_result = RombergIntegration(integrandF1,integrandF2,basic_dr,num_deep,basic_n,iaa,ibb);
    }
    else
    {
        iaa0 = iaa;
        ibb0 = ibb;
    }
    
    
    while((ibb0-iaa0) > 1)
    {
        
        //exp_v    = ibb0-iaa0;
        num_deep = (short) log2(ibb0-iaa0);
        num_deep = num_deep + 1;
        if(num_deep > max_num_Rom) num_deep = max_num_Rom;
        pow2_deep= pow(2,(num_deep-1));
        basic_n  = (ibb0-iaa0)/pow2_deep + 1;
        basic_dr = dr*pow2_deep;
        
        iaa = iaa0;
        ibb = iaa0 + (basic_n-1)*pow2_deep;
        iaa0= ibb;
        
        final_result0 = RombergIntegration(integrandF1,integrandF2,basic_dr,num_deep,basic_n,iaa,ibb);
        final_result += final_result0;
        
    }
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr1_1   = 1/dr;
        dr1_6   = dr1_1/6;
        dr2_1   = 0.5*dr1_1;
        dr2_2   = dr1_1*dr2_1;
        dr3_6   = dr1_1*dr2_2/3.0;
        
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF1[iaa0]*integrandF2[iaa0] + integrandF1[ibb0]*integrandF2[ibb0])*0.5*dr;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            coe1 = (integrandF1[ibb0+1]*integrandF2[ibb0+1] + integrandF1[iaa0]*integrandF2[iaa0]*3.0 -
                    integrandF1[ibb0]*integrandF2[ibb0]*3.0 - integrandF1[iaa0-1]*integrandF2[iaa0-1])*dr3_6/4;
            coe2 = (integrandF1[iaa0-1]*integrandF2[iaa0-1] + integrandF1[ibb0]*integrandF2[ibb0] -
                    2.0*integrandF1[iaa0]*integrandF2[iaa0])*dr2_2/3;
            coe3 = (6.0*integrandF1[ibb0]*integrandF2[ibb0] - 3.0*integrandF1[iaa0]*integrandF2[iaa0] -
                    2.0*integrandF1[iaa0-1]*integrandF2[iaa0-1] - integrandF1[ibb0+1]*integrandF2[ibb0+1])*dr1_6/2;
            coe4 = integrandF1[iaa0]*integrandF2[iaa0];
            //FF_IN[i] = dr3*(integrandF[ibb0+1]*(drr*i+dr)*(drr*i)*(drr*i-dr)/6.0 + 0.5*integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)*(drr*i-2*dr)
            //- 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i)*(drr*i-2*dr) - integrandF[iaa0-1]*(drr*i)*(drr*i-dr)*(drr*i-2*dr)/6.0);
            
            final_result0 = pow(dr,4)*coe1 + pow(dr,3)*coe2 + dr*dr*coe3 + dr*coe4;
        }
        else
        {
            coe1 = 0;
            coe2 = 0;
            coe3 = 0;
            if(iaa0 > aa)
            {
                coe1 = dr2_2*(integrandF1[iaa0-1]*integrandF2[iaa0-1] - 2.0*integrandF1[iaa0]*integrandF2[iaa0] +
                              integrandF1[ibb0]*integrandF2[ibb0])/3;
                coe2 = (integrandF1[ibb0]*integrandF2[ibb0] - integrandF1[iaa0-1]*integrandF2[iaa0-1])*dr2_1/2;
                coe3 = integrandF1[iaa0]*integrandF2[iaa0];
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0-1]*(drr*i)*(drr*i-dr) - integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)
                //+ 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i));
            }
            else if(ibb0 < bb)
            {
                coe1 = dr2_2*(integrandF1[iaa0]*integrandF2[iaa0] - 2.0*integrandF1[ibb0]*integrandF2[ibb0] +
                              integrandF1[ibb0+1]*integrandF2[ibb0+1])/3;
                coe2 = (4.0*integrandF1[ibb0]*integrandF2[ibb0] - integrandF1[ibb0+1]*integrandF2[ibb0+1] -
                        3.0*integrandF1[iaa0]*integrandF2[iaa0])*dr2_1/2;
                coe3 = integrandF1[iaa0]*integrandF2[iaa0];
                //FF_IN[i] = dr2*(0.5*integrandF[iaa0]*(drr*i-dr)*(drr*i-2*dr) - integrandF[ibb0]*(drr*i)*(drr*i-2*dr)
                //+ 0.5*integrandF[ibb0+1]*(drr*i)*(drr*i-dr));
            }
            
            final_result0 = pow(dr,3)*coe1 + dr*dr*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}

double RombergIntegration(double* integrandF,int aa,int bb,int ia,int ib,int i_R)
{
    int      basic_n,iaa,ibb,iaa0,ibb0,pow2_deep,i_Ra,i_Rb;//basic_n_mim,exp_v;
    short    max_num_Rom,num_deep;
    double   basic_dr,final_result,final_result0;
    double   coe1,coe2,coe3,coe4,dr1_1,dr1_6,dr2_6;
    
    max_num_Rom = 8; //maximum extrapolation times of Romberg
    //basic_n_mim = 17; //minimum number of grids
    
    //exp_v    = (ib-ia)/(basic_n_mim-1);
    final_result = 0;
    if(ib==ia) return final_result;
    num_deep = (short) log2(ib-ia);
    num_deep = num_deep + 1;
    if(num_deep > max_num_Rom) num_deep = max_num_Rom;
    pow2_deep= pow(2,(num_deep-1));
    basic_n  = (ib-ia)/pow2_deep + 1;
    basic_dr = dr*pow2_deep;
    
    iaa  = ia;
    ibb  = ia + (basic_n-1)*pow2_deep;
    iaa0 = ibb;
    ibb0 = ib;
    
    if((ibb-iaa) > 1)
    {
        final_result = RombergIntegration(integrandF,basic_dr,num_deep,basic_n,iaa,ibb,i_R);
    }
    else
    {
        iaa0 = iaa;
        ibb0 = ibb;
    }
    
    
    while((ibb0-iaa0) > 1)
    {
        
        //exp_v    = ibb0-iaa0;
        num_deep = (short) log2(ibb0-iaa0);
        num_deep = num_deep + 1;
        if(num_deep > max_num_Rom) num_deep = max_num_Rom;
        pow2_deep= pow(2,(num_deep-1));
        basic_n  = (ibb0-iaa0)/pow2_deep + 1;
        basic_dr = dr*pow2_deep;
        
        iaa = iaa0;
        ibb = iaa0 + (basic_n-1)*pow2_deep;
        iaa0= ibb;
        
        final_result0 = RombergIntegration(integrandF,basic_dr,num_deep,basic_n,iaa,ibb,i_R);
        final_result += final_result0;
        
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
            final_result0 = (integrandF[iaa0]*i_Ra + integrandF[ibb0]*i_Rb)*0.5*dr*dr;
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
            
            final_result0 = pow(dr,4)*coe1 + pow(dr,3)*coe2 + dr*dr*coe3 + dr*coe4;
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
            
            final_result0 = pow(dr,3)*coe1 + dr*dr*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}



double RombergIntegration(double* integrandF,int aa,int bb,int ia,int ib,double D2_j,int i_R)
{
    int      basic_n,iaa,ibb,iaa0,ibb0,pow2_deep;//basic_n_mim,exp_v;
    short    max_num_Rom,num_deep;
    double   basic_dr,final_result,final_result0,dr2,DD_2;
    double   coe1,coe2,coe3,coe4,dr1_1,dr1_6,i_D_Ra0,i_D_Rb0,i_D_Ra1,i_D_Rb1;
    
    dr2   = dr*dr;
    DD_2  = D2_j/dr2;
    
    max_num_Rom = 8; //maximum extrapolation times of Romberg
    //basic_n_mim = 17; //minimum number of grids
    
    //exp_v    = (ib-ia)/(basic_n_mim-1);
    final_result = 0;
    if(ib==ia) return final_result;
    num_deep = (short) log2(ib-ia);
    num_deep = num_deep + 1;
    if(num_deep > max_num_Rom) num_deep = max_num_Rom;
    pow2_deep= pow(2,(num_deep-1));
    basic_n  = (ib-ia)/pow2_deep + 1;
    basic_dr = dr*pow2_deep;
    
    iaa  = ia;
    ibb  = ia + (basic_n-1)*pow2_deep;
    iaa0 = ibb;
    ibb0 = ib;
    
    if((ibb-iaa) > 1)
    {
        final_result = RombergIntegration(integrandF,basic_dr,num_deep,basic_n,iaa,ibb,D2_j,i_R);
    }
    else
    {
        iaa0 = iaa;
        ibb0 = ibb;
    }
    
    
    while((ibb0-iaa0) > 1)
    {
        
        //exp_v    = ibb0-iaa0;
        num_deep = (short) log2(ibb0-iaa0);
        num_deep = num_deep + 1;
        if(num_deep > max_num_Rom) num_deep = max_num_Rom;
        pow2_deep= pow(2,(num_deep-1));
        basic_n  = (ibb0-iaa0)/pow2_deep + 1;
        basic_dr = dr*pow2_deep;
        
        iaa = iaa0;
        ibb = iaa0 + (basic_n-1)*pow2_deep;
        iaa0= ibb;
        
        final_result0 = RombergIntegration(integrandF,basic_dr,num_deep,basic_n,iaa,ibb,D2_j,i_R);
        final_result += final_result0;
        
    }
    
    
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        dr1_1   = 1/dr;
        dr1_6   = dr1_1/6;
        
        i_D_Ra0 = DD_2 - (iaa0 - i_R)*(iaa0 - i_R);
        i_D_Rb0 = DD_2 - (ibb0 - i_R)*(ibb0 - i_R);
        i_D_Ra1 = DD_2 - (iaa0-i_R-1)*(iaa0-i_R-1);
        i_D_Rb1 = DD_2 - (ibb0-i_R+1)*(ibb0-i_R+1);
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (integrandF[iaa0]*i_D_Ra0 + integrandF[ibb0]*i_D_Rb0)*0.5*dr*dr2;
        }
        else if((iaa0>aa) && (ibb0<bb))
        {
            coe1 = (integrandF[ibb0+1]*i_D_Rb1 + integrandF[iaa0]*i_D_Ra0*3.0 -
                    integrandF[ibb0]*i_D_Rb0*3.0 - integrandF[iaa0-1]*i_D_Ra1)*dr1_6/4;
            coe2 = (integrandF[iaa0-1]*i_D_Ra1 + integrandF[ibb0]*i_D_Rb0 -
                    2.0*integrandF[iaa0]*i_D_Ra0)/6;
            coe3 = (6.0*integrandF[ibb0]*i_D_Rb0 - 3.0*integrandF[iaa0]*i_D_Ra0 -
                    2.0*integrandF[iaa0-1]*i_D_Ra1 - integrandF[ibb0+1]*i_D_Rb1)*dr/12;
            coe4 = integrandF[iaa0]*i_D_Ra0*dr2;
            //FF_IN[i] = dr3*(integrandF[ibb0+1]*(drr*i+dr)*(drr*i)*(drr*i-dr)/6.0 + 0.5*integrandF[iaa0]*(drr*i+dr)*(drr*i-dr)*(drr*i-2*dr)
            //- 0.5*integrandF[ibb0]*(drr*i+dr)*(drr*i)*(drr*i-2*dr) - integrandF[iaa0-1]*(drr*i)*(drr*i-dr)*(drr*i-2*dr)/6.0);
            
            final_result0 = dr2*dr2*coe1 + dr*dr2*coe2 + dr2*coe3 + dr*coe4;
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
            
            final_result0 = dr2*dr*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}

double RombergIntegration(int aa,int bb,int ia,int ib,double D2_j)
{
    int      basic_n,iaa,ibb,iaa0,ibb0,pow2_deep;//basic_n_mim,exp_v;
    short    max_num_Rom,num_deep;
    double   basic_dr,final_result,final_result0,dr2,DD_2;
    double   coe1,coe2,coe3,i_D_Ra0,i_D_Rb0,i_D_Ra1;
    
    dr2   = dr*dr;
    DD_2  = D2_j/dr2;
    
    max_num_Rom = 8; //maximum extrapolation times of Romberg
    //basic_n_mim = 17; //minimum number of grids
    
    //exp_v    = (ib-ia)/(basic_n_mim-1);
    final_result = 0;
    if(ib==ia) return final_result;
    num_deep = (short) log2(ib-ia);
    num_deep = num_deep + 1;
    if(num_deep > max_num_Rom) num_deep = max_num_Rom;
    pow2_deep= pow(2,(num_deep-1));
    basic_n  = (ib-ia)/pow2_deep + 1;
    basic_dr = dr*pow2_deep;
    
    iaa  = ia;
    ibb  = ia + (basic_n-1)*pow2_deep;
    iaa0 = ibb;
    ibb0 = ib;
    
    if((ibb-iaa) > 1)
    {
        final_result = RombergIntegration(basic_dr,num_deep,basic_n,iaa,ibb,D2_j);
    }
    else
    {
        iaa0 = iaa;
        ibb0 = ibb;
    }
    
    while((ibb0-iaa0) > 1)
    {
        
        //exp_v    = ibb0-iaa0;
        num_deep = (short) log2(ibb0-iaa0);
        num_deep = num_deep + 1;
        if(num_deep > max_num_Rom) num_deep = max_num_Rom;
        pow2_deep= pow(2,(num_deep-1));
        basic_n  = (ibb0-iaa0)/pow2_deep + 1;
        basic_dr = dr*pow2_deep;
        
        iaa = iaa0;
        ibb = iaa0 + (basic_n-1)*pow2_deep;
        iaa0= ibb;
        
        final_result0 = RombergIntegration(basic_dr,num_deep,basic_n,iaa,ibb,D2_j);
        final_result += final_result0;
        
    }
    
    final_result0 = 0;
    if((ibb0-iaa0) == 1)
    {
        i_D_Ra0 = DD_2 - iaa0*iaa0;
        i_D_Rb0 = DD_2 - ibb0*ibb0;
        if((iaa0==aa) && (ibb0==bb))
        {
            final_result0 = (i_D_Ra0 + i_D_Rb0)*0.5*dr*dr2;
        }
        else
        {
            i_D_Ra1 = DD_2 - (iaa0-1)*(iaa0-1);
            
            coe1 = (i_D_Ra1 - 2.0*i_D_Ra0 + i_D_Rb0)/6;
            coe2 = (i_D_Rb0 - i_D_Ra1)*dr/4;
            coe3 = i_D_Ra0*dr2;
            
            final_result0 = dr2*dr*coe1 + dr2*coe2 + dr*coe3;
        }
    }
    
    final_result += final_result0;
    
    return final_result;
}




//num_deep means the maximum deepth before interpolation
//basic_dr is the step size in the first layer
double RombergIntegration(double* integrandF,double basic_dr,short num_deep,int basic_n,int iaa,int ibb)
{
    short    i_deep;
    double   dr_i,DD,final_result,error;
    int      i_grid_0,n_num,n_num_i,cc;
    short*   pow2;
    double*  T;
    
    error = errTol*1.0E-6;
    
    pow2     = new short[num_deep]();
    T        = new double[num_deep]();
    
    
    pow2[num_deep-1] = 1;
    for(short i=(num_deep-2); i>=0; --i)
    {
        pow2[i] = 2*pow2[i+1];
    }
    
    n_num    = basic_n - 1;
    i_grid_0 = iaa;
    T[0] = (integrandF[iaa] + integrandF[ibb])*0.5;
    for(int j=1; j<n_num; ++j)
    {
        i_grid_0 += pow2[0];
        T[0]     += integrandF[i_grid_0];
    }
    T[0] = T[0]*basic_dr;
    
    
    i_deep = 1;
    dr_i   = basic_dr;
    n_num_i= n_num;
    DD     = 100;
    while(i_deep < num_deep)
    {
        dr_i    = dr_i*0.5;
        i_grid_0= iaa + pow2[i_deep];
        pow2[i_deep] = pow2[i_deep]*2;
        
        T[i_deep] = integrandF[i_grid_0];
        for(int j=1; j<n_num_i; ++j)
        {
            i_grid_0  += pow2[i_deep];
            T[i_deep] += integrandF[i_grid_0];
        }
        T[i_deep] = T[i_deep-1]*0.5 + T[i_deep]*dr_i;
        
        cc = 4;
        for(short k=1; k<=i_deep; ++k)
        {
            DD = (T[i_deep+1-k] - T[i_deep-k])/(cc - 1);
            T[i_deep-k] = T[i_deep+1-k] + DD;
            cc = 4*cc;
            
            //std::cout<<"i_deep="<<i_deep<<" k="<<k<<" DD="<<DD<<std::endl;
            
            if(fabs(DD) < error)
            {
                
                final_result = T[i_deep-k];

                delete [] pow2;
                delete [] T;
                return final_result;
            }
        }
        ++i_deep;
        n_num_i *= 2;
    }
    final_result = T[0];
    
    delete [] pow2;
    delete [] T;
    
    return final_result;
}



double RombergIntegration(double* integrandF1,double* integrandF2,double basic_dr,short num_deep,int basic_n,int iaa,int ibb)
{
    short    i_deep;
    double   dr_i,DD,final_result,error;
    int      i_grid_0,n_num,n_num_i,cc;
    short*   pow2;
    double*  T;
    
    error = errTol*1.0E-6;
    
    pow2     = new short[num_deep]();
    T        = new double[num_deep]();
    
    
    pow2[num_deep-1] = 1;
    for(short i=(num_deep-2); i>=0; --i)
    {
        pow2[i] = 2*pow2[i+1];
    }
    
    n_num    = basic_n - 1;
    i_grid_0 = iaa;
    T[0] = (integrandF1[iaa]*integrandF2[iaa] + integrandF1[ibb]*integrandF2[ibb])*0.5;
    for(int j=1; j<n_num; ++j)
    {
        i_grid_0 += pow2[0];
        T[0]     += integrandF1[i_grid_0]*integrandF2[i_grid_0];
    }
    T[0] = T[0]*basic_dr;
    
    
    i_deep = 1;
    dr_i   = basic_dr;
    n_num_i= n_num;
    DD     = 100;
    while(i_deep < num_deep)
    {
        dr_i    = dr_i*0.5;
        i_grid_0= iaa + pow2[i_deep];
        pow2[i_deep] = pow2[i_deep]*2;
        
        T[i_deep] = integrandF1[i_grid_0]*integrandF2[i_grid_0];
        for(int j=1; j<n_num_i; ++j)
        {
            i_grid_0  += pow2[i_deep];
            T[i_deep] += integrandF1[i_grid_0]*integrandF2[i_grid_0];
        }
        T[i_deep] = T[i_deep-1]*0.5 + T[i_deep]*dr_i;
        
        cc = 4;
        for(short k=1; k<=i_deep; ++k)
        {
            DD = (T[i_deep+1-k] - T[i_deep-k])/(cc - 1);
            T[i_deep-k] = T[i_deep+1-k] + DD;
            cc = 4*cc;
            
            if(fabs(DD) < error)
            {
                
                final_result = T[i_deep-k];
                
                delete [] pow2;
                delete [] T;
                return final_result;
            }
        }
        ++i_deep;
        n_num_i *= 2;
    }
    final_result = T[0];
    
    delete [] pow2;
    delete [] T;
    
    return final_result;
}



double RombergIntegration(double* integrandF,double basic_dr,short num_deep,int basic_n,int iaa,int ibb,int i_R)
{
    short    i_deep;
    double   dr_i,DD,final_result,error;
    int      i_grid_0,n_num,n_num_i,cc;
    short*   pow2;
    double*  T;
    
    error = errTol*1.0E-6;
    //RR    = dr*i_R;
    
    pow2     = new short[num_deep]();
    T        = new double[num_deep]();
    
    
    pow2[num_deep-1] = 1;
    for(short i=(num_deep-2); i>=0; --i)
    {
        pow2[i] = 2*pow2[i+1];
    }
    
    n_num    = basic_n - 1;
    i_grid_0 = iaa;
    T[0] = (integrandF[iaa]*(iaa-i_R) + integrandF[ibb]*(ibb-i_R))*0.5;
    for(int j=1; j<n_num; ++j)
    {
        i_grid_0 += pow2[0];
        T[0]     += (integrandF[i_grid_0]*(i_grid_0-i_R));
    }
    T[0] = T[0]*basic_dr*dr;
    
    
    i_deep = 1;
    dr_i   = basic_dr;
    n_num_i= n_num;
    DD     = 100;
    while(i_deep < num_deep)
    {
        dr_i    = dr_i*0.5;
        i_grid_0= iaa + pow2[i_deep];
        pow2[i_deep] = pow2[i_deep]*2;
        
        T[i_deep] = integrandF[i_grid_0]*(i_grid_0-i_R);
        for(int j=1; j<n_num_i; ++j)
        {
            i_grid_0  += pow2[i_deep];
            T[i_deep] += (integrandF[i_grid_0]*(i_grid_0-i_R));
        }
        T[i_deep] = T[i_deep-1]*0.5 + T[i_deep]*dr_i*dr;
        
        cc = 4;
        for(short k=1; k<=i_deep; ++k)
        {
            DD = (T[i_deep+1-k] - T[i_deep-k])/(cc - 1);
            T[i_deep-k] = T[i_deep+1-k] + DD;
            cc = 4*cc;
            
            if(fabs(DD) < error)
            {
                final_result = T[i_deep-k];
                
                delete [] pow2;
                delete [] T;
                return final_result;
            }
        }
        ++i_deep;
        n_num_i *= 2;
    }
    final_result = T[0];
    
    delete [] pow2;
    delete [] T;
    
    return final_result;
}


double RombergIntegration(double* integrandF,double basic_dr,short num_deep,int basic_n,int iaa,int ibb,
                          double D2_j,int i_R)
{
    short    i_deep;
    double   dr_i,DD,final_result,error,dr2,DD_2;
    int      i_grid_0,n_num,n_num_i,cc,delta_R;
    short*   pow2;
    double*  T;
    
    error = errTol*1.0E-6;
    dr2   = dr*dr;
    DD_2  = D2_j/dr2;
    //RR    = dr*i_R;
    
    pow2     = new short[num_deep]();
    T        = new double[num_deep]();
    
    
    pow2[num_deep-1] = 1;
    for(short i=(num_deep-2); i>=0; --i)
    {
        pow2[i] = 2*pow2[i+1];
    }
    
    n_num    = basic_n - 1;
    i_grid_0 = iaa;
    T[0] = (integrandF[iaa]*(DD_2-(iaa-i_R)*(iaa-i_R)) +
            integrandF[ibb]*(DD_2-(ibb-i_R)*(ibb-i_R)))*0.5;
    for(int j=1; j<n_num; ++j)
    {
        i_grid_0 += pow2[0];
        delta_R   = i_grid_0-i_R;
        T[0]     += (integrandF[i_grid_0]*(DD_2-delta_R*delta_R));
    }
    T[0] = T[0]*basic_dr*dr2;
    
    
    i_deep = 1;
    dr_i   = basic_dr;
    n_num_i= n_num;
    DD     = 100;
    while(i_deep < num_deep)
    {
        dr_i    = dr_i*0.5;
        i_grid_0= iaa + pow2[i_deep];
        pow2[i_deep] = pow2[i_deep]*2;
        
        T[i_deep] = integrandF[i_grid_0]*(DD_2-(i_grid_0-i_R)*(i_grid_0-i_R));
        for(int j=1; j<n_num_i; ++j)
        {
            i_grid_0  += pow2[i_deep];
            delta_R    = i_grid_0-i_R;
            T[i_deep] += (integrandF[i_grid_0]*(DD_2-delta_R*delta_R));
        }
        T[i_deep] = T[i_deep-1]*0.5 + T[i_deep]*dr_i*dr2;
        
        cc = 4;
        for(short k=1; k<=i_deep; ++k)
        {
            DD = (T[i_deep+1-k] - T[i_deep-k])/(cc - 1);
            T[i_deep-k] = T[i_deep+1-k] + DD;
            cc = 4*cc;
            
            if(fabs(DD) < error)
            {
                final_result = T[i_deep-k];
                
                delete [] pow2;
                delete [] T;
                return final_result;
            }
        }
        ++i_deep;
        n_num_i *= 2;
    }
    final_result = T[0];
    
    delete [] pow2;
    delete [] T;
    
    return final_result;
}

double RombergIntegration(double basic_dr,short num_deep,int basic_n,int iaa,int ibb,double D2_j)
{
    short    i_deep;
    double   dr_i,DD,final_result,error,dr2,DD_2;
    int      i_grid_0,n_num,n_num_i,cc;
    short*   pow2;
    double*  T;
    
    error = errTol*1.0E-6;
    dr2   = dr*dr;
    DD_2  = D2_j/dr2;
    //RR    = dr*i_R;
    
    pow2     = new short[num_deep]();
    T        = new double[num_deep]();
    
    
    pow2[num_deep-1] = 1;
    for(short i=(num_deep-2); i>=0; --i)
    {
        pow2[i] = 2*pow2[i+1];
    }
    
    n_num    = basic_n - 1;
    i_grid_0 = iaa;
    T[0] = ((DD_2-iaa*iaa)+(DD_2-ibb*ibb))*0.5;
    for(int j=1; j<n_num; ++j)
    {
        i_grid_0 += pow2[0];
        T[0]     += (DD_2-i_grid_0*i_grid_0);
    }
    T[0] = T[0]*basic_dr*dr2;
    
    
    i_deep = 1;
    dr_i   = basic_dr;
    n_num_i= n_num;
    DD     = 100;
    while(i_deep < num_deep)
    {
        dr_i    = dr_i*0.5;
        i_grid_0= iaa + pow2[i_deep];
        pow2[i_deep] = pow2[i_deep]*2;
        
        T[i_deep] = DD_2-i_grid_0*i_grid_0;
        for(int j=1; j<n_num_i; ++j)
        {
            i_grid_0  += pow2[i_deep];
            T[i_deep] += (DD_2-i_grid_0*i_grid_0);
        }
        T[i_deep] = T[i_deep-1]*0.5 + T[i_deep]*dr_i*dr2;
        
        cc = 4;
        for(short k=1; k<=i_deep; ++k)
        {
            DD = (T[i_deep+1-k] - T[i_deep-k])/(cc - 1);
            T[i_deep-k] = T[i_deep+1-k] + DD;
            cc = 4*cc;
            
            if(fabs(DD) < error)
            {
                final_result = T[i_deep-k];
                
                delete [] pow2;
                delete [] T;
                return final_result;
            }
        }
        ++i_deep;
        n_num_i *= 2;
    }
    final_result = T[0];
    
    delete [] pow2;
    delete [] T;
    
    return final_result;
}



                     

                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     

                     
                     
                     
                     
                     
                    

