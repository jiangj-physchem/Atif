//***********potential from electrostatistic correlation****************//
#include "clibrary.h"
#include "psichargeshell.h"
#include "constantnum.h"

extern double dr;
extern double BJ;
extern int    DMAXg;
extern int    DMAX;
extern short  nspecies;
//using namespace std;

void PsiChargeShell(double* B,float* Z,double*** Psi_IJ)
{
    double  PiBJ,temp1,temp2,R_t2,R_t3,R_t,R_t0;
    double  temp32,temp33,temp42,temp43,temp3,temp4,temp5;
    double*   SubK;
    double**  coe11_ij;
    double**  coe12_ij;
    double**  coe2_ij;
    double**  coe3_ij;
    double**  coe4_ij;
    double**  coe5_ij;
    double**  temp11_ij;
    double**  temp12_ij;
    double**  temp2_ij;
    double**  DBIJ;
    int**     BIJ;
    short   hspecies;
    int     kin,kip,kk;
    
    hspecies = nspecies - 1;
    
    BIJ     = new int*[hspecies]();
    SubK    = new double[3]();
    DBIJ    = new double*[hspecies]();
    coe11_ij= new double*[hspecies]();
    coe12_ij= new double*[hspecies]();
    coe2_ij = new double*[hspecies]();
    coe3_ij = new double*[hspecies]();
    coe4_ij = new double*[hspecies]();
    coe5_ij = new double*[hspecies]();
    temp2_ij= new double*[hspecies]();
    temp11_ij= new double*[hspecies]();
    temp12_ij= new double*[hspecies]();
    
    SubK[0]  = 0.0;
    SubK[1]  = 0.2254033307585*dr;
    SubK[2]  = 0.7745966692415*dr;
    
    for(short i=0; i<hspecies; ++i)
    {
        //B[i] = D[i]*0.5;
        BIJ[i]     = new int[hspecies]();
        DBIJ[i]    = new double[hspecies]();
        coe11_ij[i]= new double[hspecies]();
        coe12_ij[i]= new double[hspecies]();
        coe2_ij[i] = new double[hspecies]();
        coe3_ij[i] = new double[hspecies]();
        coe4_ij[i] = new double[hspecies]();
        coe5_ij[i] = new double[hspecies]();
        temp2_ij[i]= new double[hspecies]();
        temp11_ij[i]= new double[hspecies]();
        temp12_ij[i]= new double[hspecies]();
    }
    PiBJ= Pi*BJ;
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = PiBJ*Z[i]*Z[j];
            temp2 = temp1/(6*B[i]*B[j]);
            temp3 = fabs(B[i]-B[j]);
            temp4 = B[i] + B[j];
            temp32= temp3*temp3;
            temp33= temp32*temp3;
            temp42= temp4*temp4;
            temp43= temp42*temp4;
            
            temp5 = temp2*(3*temp42*temp3-3*temp4*temp32-temp43+temp33) - 2*temp1*temp3;
            
            temp11_ij[i][j]= (temp1*temp32)/B[i] + temp5;
            temp12_ij[i][j]= (temp1*temp32)/B[j] + temp5;
            
            temp2_ij[i][j] = -(temp2*temp43);
            coe11_ij[i][j] = -(temp1/B[i]);
            coe12_ij[i][j] = -(temp1/B[j]);
            coe2_ij[i][j]  = 2*temp1;
            coe3_ij[i][j]  = 3*temp2*temp42;
            coe4_ij[i][j]  = -(3*temp2*temp4);
            coe5_ij[i][j]  = temp2;
            
            BIJ[i][j] = round(temp4/dr);
            DBIJ[i][j]= temp3;
            //DBIJ[i][j]= round(temp3/dr);
        }
    }
    
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            kin = -BIJ[i][j];
            kip = BIJ[i][j];
            
            if(B[i] > B[j])
            {
                for(int k=kin; k<kip; ++k)
                {
                    R_t0= k*dr;
                    kk  = 3*(k + DMAXg);
                    for(int k1=kk; k1<(kk+3); ++k1)
                    {
                        R_t = R_t0 + SubK[k1-kk];
                        if(R_t < 0) R_t = -R_t;
                        R_t2 = R_t*R_t;
                        R_t3 = R_t*R_t2;
                        
                        if(R_t < DBIJ[i][j])
                        {
                            Psi_IJ[i][j][k1] = temp11_ij[i][j] + coe11_ij[i][j]*R_t2 + coe2_ij[i][j]*R_t;
                        }
                        else
                        {
                            Psi_IJ[i][j][k1] = temp2_ij[i][j] + coe3_ij[i][j]*R_t + coe4_ij[i][j]*R_t2 + coe5_ij[i][j]*R_t3;
                        }
                        
                    }
                }
                
                kk  = 3*(kip + DMAXg);
                R_t = kip*dr;
                R_t2= R_t*R_t;
                R_t3= R_t*R_t2;
                
                Psi_IJ[i][j][kk] = temp2_ij[i][j] + coe3_ij[i][j]*R_t + coe4_ij[i][j]*R_t2 + coe5_ij[i][j]*R_t3;
            }
            else
            {
                for(int k=kin; k<kip; ++k)
                {//loop
                    R_t0= k*dr;
                    kk  = 3*(k + DMAXg);
                    for(int k1=kk; k1<(kk+3); ++k1)
                    {
                        R_t = R_t0 + SubK[k1-kk];
                        if(R_t < 0) R_t = -R_t;
                        R_t2 = R_t*R_t;
                        R_t3 = R_t*R_t2;
                        
                        if(R_t < DBIJ[i][j])
                        {
                            Psi_IJ[i][j][k1] = temp12_ij[i][j] + coe12_ij[i][j]*R_t2 + coe2_ij[i][j]*R_t;
                        }
                        else
                        {
                            Psi_IJ[i][j][k1] = temp2_ij[i][j] + coe3_ij[i][j]*R_t + coe4_ij[i][j]*R_t2 + coe5_ij[i][j]*R_t3;
                        }
                    }
                }//loop
                
                kk  = 3*(kip + DMAXg);
                R_t = kip*dr;
                R_t2= R_t*R_t;
                R_t3= R_t*R_t2;
                
                Psi_IJ[i][j][kk] = temp2_ij[i][j] + coe3_ij[i][j]*R_t + coe4_ij[i][j]*R_t2 + coe5_ij[i][j]*R_t3;

            }
            
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] coe11_ij[i];
        delete [] coe12_ij[i];
        delete [] coe2_ij[i];
        delete [] coe3_ij[i];
        delete [] coe4_ij[i];
        delete [] coe5_ij[i];
        delete [] temp11_ij[i];
        delete [] temp12_ij[i];
        delete [] temp2_ij[i];
        delete [] BIJ[i];
        delete [] DBIJ[i];
    }
    delete [] coe11_ij;
    delete [] coe12_ij;
    delete [] coe2_ij;
    delete [] coe3_ij;
    delete [] coe4_ij;
    delete [] coe5_ij;
    delete [] temp11_ij;
    delete [] temp12_ij;
    delete [] temp2_ij;
    delete [] BIJ;
    delete [] DBIJ;
    delete [] SubK;
}

void PsiChargeShellJJ(double* B,double* TB,float* Z,double*** Psi_IJ)
{
    double    PiBJ,R_t,R_t0,temp1,temp2,temp3,temp4;
    double*   SubK;
    double**  coe_ij;
    double**  temp_ij;
    double**  DBIJ;
    double**  BIJ;
    int**     TBIJ;
    short   hspecies;
    int     kin,kip,kk;
    
    hspecies = nspecies - 1;
    
    TBIJ    = new int*[hspecies]();
    SubK    = new double[3]();
    DBIJ    = new double*[hspecies]();
    BIJ     = new double*[hspecies]();
    coe_ij  = new double*[hspecies]();
    temp_ij = new double*[hspecies]();
    

    
    SubK[0]  = 0.0;
    SubK[1]  = 0.2254033307585*dr;
    SubK[2]  = 0.7745966692415*dr;
    
    for(short i=0; i<hspecies; ++i)
    {
        TBIJ[i]    = new int[hspecies]();
        DBIJ[i]    = new double[hspecies]();
        BIJ[i]     = new double[hspecies]();
        coe_ij[i]  = new double[hspecies]();
        temp_ij[i] = new double[hspecies]();
    }
    PiBJ= Pi*BJ;
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = (PiBJ*Z[i]*Z[j])/(6*B[i]*B[j]);
            temp2 = fabs(B[i]-B[j]);
            temp3 = B[i] + B[j];
            temp4 = TB[i] + TB[j];
            
            TBIJ[i][j] = round(temp4/dr);
            temp4      = TBIJ[i][j]*dr;
            
            temp_ij[i][j]= temp1*(temp3-temp4)*(temp3-temp4)*(temp3-temp4);
            coe_ij[i][j] = temp1;
            
            DBIJ[i][j] = temp2;
            BIJ[i][j]  = temp3;
        }
    }
    
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            kin = -TBIJ[i][j];
            kip = TBIJ[i][j];
            
            for(int k=kin; k<kip; ++k)
            {
                R_t0= k*dr;
                kk  = 3*(k + DMAXg);
                for(int k1=kk; k1<(kk+3); ++k1)
                {
                    R_t = R_t0 + SubK[k1-kk];
                    if(R_t < 0) R_t = -R_t;
                    
                    temp1 = R_t - BIJ[i][j];
                    temp2 = R_t - DBIJ[i][j];
                    
                    temp3 = temp1*temp1*temp1;
                    temp4 = temp2*temp2*temp2;
                    
                    if(R_t < DBIJ[i][j])
                    {
                        Psi_IJ[i][j][k1] = temp_ij[i][j] + coe_ij[i][j]*(temp3-temp4);
                    }
                    else
                    {
                        Psi_IJ[i][j][k1] = temp_ij[i][j] + coe_ij[i][j]*temp3;
                    }
                    
                }
            }
            
            kk  = 3*(kip + DMAXg);
            R_t = kip*dr;
            temp1 = R_t - BIJ[i][j];
            temp3 = temp1*temp1*temp1;
            
            Psi_IJ[i][j][kk] = temp_ij[i][j] + coe_ij[i][j]*temp3;
            
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] coe_ij[i];
        delete [] temp_ij[i];
        delete [] TBIJ[i];
        delete [] DBIJ[i];
        delete [] BIJ[i];
    }
    delete [] coe_ij;
    delete [] temp_ij;
    delete [] TBIJ;
    delete [] DBIJ;
    delete [] BIJ;
    delete [] SubK;
}





void PsiChargeShellDirk(double gama,float* D,float* Z,double*** Psi_IJ)
{
    double  PiBJ,temp1,temp2,R_t2,R_t3,R_t,R_t0,gammar;
    double  temp32,temp33,temp42,temp3,temp4;
    double*   B;
    double*   R;
    double*   SubK;
    double**  coe1_ij;
    double**  coe2_ij;
    double**  coe3_ij;
    double**  temp1_ij;
    int**     RIJ;
    short   hspecies;
    int     kin,kip,kk,igammar;
    
    hspecies = nspecies - 1;
    
    B       = new double[hspecies]();
    R       = new double[hspecies]();
    RIJ     = new int*[hspecies]();
    SubK    = new double[3]();
    coe1_ij = new double*[hspecies]();
    coe2_ij = new double*[hspecies]();
    coe3_ij = new double*[hspecies]();
    temp1_ij= new double*[hspecies]();
    
    SubK[0]  = 0.0;
    SubK[1]  = 0.2254033307585*dr;
    SubK[2]  = 0.7745966692415*dr;
    
    gammar   = 0.5/gama;
    igammar  = round(gammar/dr);
    gammar   = igammar*dr;
    for(short i=0; i<hspecies; ++i)
    {
        R[i] = D[i]*0.5;
        B[i] = R[i] + gammar;
        RIJ[i]     = new int[hspecies]();
        coe1_ij[i] = new double[hspecies]();
        coe2_ij[i] = new double[hspecies]();
        coe3_ij[i] = new double[hspecies]();
        temp1_ij[i]= new double[hspecies]();
    }
    PiBJ= Pi*BJ;
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = PiBJ*Z[i]*Z[j];
            temp2 = temp1/(6*B[i]*B[j]);
            temp3 = R[i] + R[j];
            temp4 = B[i] + B[j];
            temp32= temp3*temp3;
            temp33= temp32*temp3;
            temp42= temp4*temp4;
            
            temp1_ij[i][j] = temp2*(3*temp4*temp32-3*temp42*temp3-temp33);
            
            coe1_ij[i][j]  = temp2;
            coe2_ij[i][j]  = -(3*temp2*temp4);
            coe3_ij[i][j]  = 3*temp2*temp42;

            RIJ[i][j] = round(temp3/dr);
        }
    }
    
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            kin = -RIJ[i][j];
            kip = RIJ[i][j];
            
            for(int k=kin; k<kip; ++k)
            {
                R_t0= k*dr;
                kk  = 3*(k + DMAX);
                for(int k1=kk; k1<(kk+3); ++k1)
                {
                    R_t = R_t0 + SubK[k1-kk];
                    if(R_t < 0) R_t = -R_t;
                    R_t2 = R_t*R_t;
                    R_t3 = R_t*R_t2;
                    
                    Psi_IJ[i][j][k1] = temp1_ij[i][j] + coe3_ij[i][j]*R_t + coe2_ij[i][j]*R_t2 + coe1_ij[i][j]*R_t3;
                }
            }
            
            kk  = 3*(kip + DMAX);
            R_t = kip*dr;
            R_t2= R_t*R_t;
            R_t3= R_t*R_t2;
            
            Psi_IJ[i][j][kk] = temp1_ij[i][j] + coe3_ij[i][j]*R_t + coe2_ij[i][j]*R_t2 + coe1_ij[i][j]*R_t3;
            
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] coe1_ij[i];
        delete [] coe2_ij[i];
        delete [] coe3_ij[i];
        delete [] temp1_ij[i];
        delete [] RIJ[i];
    }
    delete [] coe1_ij;
    delete [] coe2_ij;
    delete [] coe3_ij;
    delete [] temp1_ij;
    delete [] RIJ;
    delete [] R;
    delete [] B;
    delete [] SubK;
}

/*
void PsiChargeShell(float* D,float* Z,double*** Psi_IJ)
{
    double  PiBJ,temp1,temp2,R_t2,R_t3,R_t;
    double  temp32,temp33,temp42,temp43,temp3,temp4,temp5;
    float*   B;
    double**  coe11_ij;
    double**  coe12_ij;
    double**  coe2_ij;
    double**  coe3_ij;
    double**  coe4_ij;
    double**  coe5_ij;
    double**  temp11_ij;
    double**  temp12_ij;
    double**  temp2_ij;
    int**     BIJ;
    int**     DBIJ;
    short   hspecies;
    int     kin,kip,iR_t,k_abs;
    
    hspecies = nspecies - 1;
    
    B       = new float[hspecies]();
    BIJ     = new int*[hspecies]();
    DBIJ    = new int*[hspecies]();
    coe11_ij= new double*[hspecies]();
    coe12_ij= new double*[hspecies]();
    coe2_ij = new double*[hspecies]();
    coe3_ij = new double*[hspecies]();
    coe4_ij = new double*[hspecies]();
    coe5_ij = new double*[hspecies]();
    temp2_ij= new double*[hspecies]();
    temp11_ij= new double*[hspecies]();
    temp12_ij= new double*[hspecies]();
    
    for(short i=0; i<hspecies; ++i)
    {
        B[i] = D[i]*0.5;
        BIJ[i]     = new int[hspecies]();
        DBIJ[i]    = new int[hspecies]();
        coe11_ij[i]= new double[hspecies]();
        coe12_ij[i]= new double[hspecies]();
        coe2_ij[i] = new double[hspecies]();
        coe3_ij[i] = new double[hspecies]();
        coe4_ij[i] = new double[hspecies]();
        coe5_ij[i] = new double[hspecies]();
        temp2_ij[i]= new double[hspecies]();
        temp11_ij[i]= new double[hspecies]();
        temp12_ij[i]= new double[hspecies]();
    }
    PiBJ= Pi*BJ;
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            temp1 = PiBJ*Z[i]*Z[j];
            temp2 = temp1/(6*B[i]*B[j]);
            temp3 = fabs(B[i]-B[j]);
            temp4 = B[i] + B[j];
            temp32= temp3*temp3;
            temp33= temp32*temp3;
            temp42= temp4*temp4;
            temp43= temp42*temp4;
            
            temp5 = temp2*(3*temp42*temp3-3*temp4*temp32-temp43+temp33) - 2*temp1*temp3;
            
            temp11_ij[i][j]= (temp1*temp32)/B[i] + temp5;
            temp12_ij[i][j]= (temp1*temp32)/B[j] + temp5;
            
            temp2_ij[i][j] = -(temp2*temp43);
            coe11_ij[i][j] = -(temp1/B[i]);
            coe12_ij[i][j] = -(temp1/B[j]);
            coe2_ij[i][j]  = 2*temp1;
            coe3_ij[i][j]  = 3*temp2*temp42;
            coe4_ij[i][j]  = -(3*temp2*temp4);
            coe5_ij[i][j]  = temp2;
            
            BIJ[i][j] = round(temp4/dr);
            DBIJ[i][j]= round(temp3/dr);
        }
    }
    
    for(short i=0; i<hspecies; ++i)
    {
        for(short j=0; j<hspecies; ++j)
        {
            kin = -BIJ[i][j];
            kip = BIJ[i][j];
            
            if(B[i] > B[j])
            {
                for(int k=kin; k<=kip; ++k)
                {
                    k_abs= k;
                    if(k < 0) k_abs = -k_abs;
                    R_t  = k_abs*dr;
                    iR_t = k + DMAX;
                    
                    R_t2 = R_t*R_t;
                    R_t3 = R_t*R_t2;
                    if(k_abs < DBIJ[i][j])
                    {
                        Psi_IJ[i][j][iR_t] = temp11_ij[i][j] + coe11_ij[i][j]*R_t2 + coe2_ij[i][j]*R_t;
                    }
                    else
                    {
                        Psi_IJ[i][j][iR_t] = temp2_ij[i][j] + coe3_ij[i][j]*R_t + coe4_ij[i][j]*R_t2 + coe5_ij[i][j]*R_t3;
                    }
                }
            }
            else
            {
                for(int k=kin; k<=kip; ++k)
                {
                    k_abs= k;
                    if(k < 0) k_abs = -k_abs;
                    R_t  = k_abs*dr;
                    iR_t = k + DMAX;
                    
                    R_t2 = R_t*R_t;
                    R_t3 = R_t*R_t2;
                    if(k_abs < DBIJ[i][j])
                    {
                        Psi_IJ[i][j][iR_t] = temp12_ij[i][j] + coe12_ij[i][j]*R_t2 + coe2_ij[i][j]*R_t;
                    }
                    else
                    {
                        Psi_IJ[i][j][iR_t] = temp2_ij[i][j] + coe3_ij[i][j]*R_t + coe4_ij[i][j]*R_t2 + coe5_ij[i][j]*R_t3;
                    }
                }
            }
            
        }
    }
    
    
    
    for(short i=0; i<hspecies; ++i)
    {
        delete [] coe11_ij[i];
        delete [] coe12_ij[i];
        delete [] coe2_ij[i];
        delete [] coe3_ij[i];
        delete [] coe4_ij[i];
        delete [] coe5_ij[i];
        delete [] temp11_ij[i];
        delete [] temp12_ij[i];
        delete [] temp2_ij[i];
        delete [] BIJ[i];
        delete [] DBIJ[i];
    }
    delete [] coe11_ij;
    delete [] coe12_ij;
    delete [] coe2_ij;
    delete [] coe3_ij;
    delete [] coe4_ij;
    delete [] coe5_ij;
    delete [] temp11_ij;
    delete [] temp12_ij;
    delete [] temp2_ij;
    delete [] BIJ;
    delete [] DBIJ;
    delete [] B;
}
*/
