//********************************************************************//
//  This code is for the polyelectrolytes density functional theory:  //
// This code is workable for one or two kind of charged or neutral    //
//  diblock copolymer with or without salts; asymmetric or symmetric; //
//Authors: Jian Jiang  (copyright owner: Jian Jiang)                  //
//Date: June 3, 2017                                                  //
//********************************************************************//

#include "main.h"

using namespace std;
using namespace systemset;



//*** Global variable ***********************************************//
double errTol; //the error tolerence
double BJ;     //the Bjerrum length
double kapa;   //the inverse of Deby length
double boundary;
double dr;     //step length of the discretization
double size;   //the size of the simulation box
double depthW;
double C_psi;
short  nspecies;
short* mp;//monomer # on copolymer i
short* nb;//# of blocks in copolymer i;
short** mb;//# of monomers in each block
int ngrid; //the number of grids
int ngrid_b; //the number of grids: the left bounday of bulk region
int ngrid_m; //the number of grids: the middle between two surfaces
int ULIM; //the maxium uper limit of intergral
int LLIM; //the minimum lower limit of intergral
int DMAX;
int DMAXg;
short neutsys;
float ZMax;

double* rhoB_S;
double* temp_S;
double* DL1;
double* DR1;
double* DL2;
double* DR2;
double* GRC1;
double* GRC2;




//******************************************************************//

int main()
{
    //**************************declare variables*******************//
    ios::sync_with_stdio(false);
    cin.tie(0);
    //srand((unsigned int)time(0));//set seed for rand()
    //**************************************************************//
    double    mixCoe,maxItera;
    double    nm,mix_f,mix_f0,err,temp,temp1,temp2;
    double    phi,phi0,sigma,sigma0,epsilonS,epsilonW,lunit,rhot,size0,errTol0,errTol1,errTol2,errTol3;
    double    mixCoeMax0,mixCoeMax1,mixCoeMin0,mixCoeMin1,err_temp1,err_std,std_temp;
    double    alpha,time_backup,eta,R_e1,R_e2,Ener_tot; //,variable
    //double    P,P1,P2;
    double    f,dcharge,deltaPhi,phL_L,sigma_t,sigma_t0;
    double    coe0,coe1,coe2,coe3,coe4,dr0,coe11;
    double    fsigma,rhoBM1,rhoBM2,eta_t,stiffness_b,gama,P_bulk,QQ;
    float     SortTem;
    //double*   Ener_Save;
    double*   etar;
    double*   Psi;
    double*   mu;
    double*   rhoB;
    double*   rhoB1;
    double*   epsilon_b;
    float*    D;
    float*    Z;
    float*    uWall;
    double*   B;
    double*   TB;
    double**  ATT;
    double**  Ext;
    double**  pairEner;
    double**  pairEner1;
    double**  rho;
    double**  rho1;
    double**  rho_temp1;
    double**  rho_temp2;
    double*** BesselZero;
    double*** Psi_IJ;
    short     iCode,iSort0,iSort,iDMin,nbb,nspecies1,ii2,jj2,ii3,nb1,nb2,mp1,mp2,nblocks,hspecies;//,copy_yes;
    short     AccelKeep,nblocks1;
    int       thresh_copy,thresh_in,thresh_out,thresh_nor;
    int       iloop,ngridm1,DMAX1;//,num_max;
    int       iter,iback,cstep,nBessel,NRe;
    //int        depth1,depth2;
    short*    DSort;
    float**   zb;
    int** depthLR;
    int*  LLI;
    int*  ULI;

    char fileName[1000] = {};

    SysInfo    sysInfo;
    IteraInfo  iteraInfo;
    MethodGeom MethG0;
    Species*  species;
    PolyInfo* pInfo;
    
    string filePath;
    string fileILoop= "temporary_iloop_";
    string fileTemp = "temporary_density_";
    string fileConverge = "converge_record_";
    string fileIterative = "iterative_";
    string* MODEL;
    string* MethG;
    //ofstream density_temp;
    //ofstream iloop_temp;
    ofstream iterative_out;
    ofstream contact_theorem;
    ofstream vari_press;
    ofstream parameters_s;
    ofstream converge_rec;
    ifstream inFile;
    clock_t  start_back, finish_back;
    //*************************declare variables*******************//
    pInfo= new PolyInfo[2];
    MODEL= new string[2];
    MethG= new string[4];
    mp   = new short[2]();
    epsilon_b = new double[2]();
    readSystem(pInfo,fileName,MethG0);
    
    MODEL[0] = pInfo[0].MODEL;
    MODEL[1] = pInfo[1].MODEL;
    MethG[0] = MethG0.Method;
    MethG[1] = MethG0.Geomet;
    MethG[2] = MethG0.nSurfs;
    MethG[3] = MethG0.exPoten;
    
    //bending potential
    epsilon_b[0]=  pInfo[0].epsilon_b;
    epsilon_b[1]=  pInfo[1].epsilon_b;
    nspecies1 = 5 + pInfo[0].nb + pInfo[1].nb;
    //nspecies: 0 --  (pInfo[0].nb-1)            polymer 1
    //nspecies: pInfo[0].nb -- (pInfo[1].nb-1)   polymer 2
    //nspecies: pInfo[1].nb -- (pInfo[1].nb+1)   salts
    //nspecies: pInfo[1].nb+2 -- (pInfo[1].nb+3) counterions
    //nspecies: pInfo[1].nb+4                    hard sphere
    
    nb    = new short[2]();
    nb[0] = pInfo[0].nb;
    nb[1] = pInfo[1].nb;
    nb1   = pInfo[0].nb;
    nb2   = pInfo[1].nb;
    nblocks1= nb1 + nb2;
    
    mb   = new short*[2]();
    mb[0]= new short[nb1]();
    mb[1]= new short[nb2]();
    
    zb   = new float*[2]();
    zb[0]= new float[nb1]();
    zb[1]= new float[nb2]();
    
    
    
    rhoB1    = new double[nspecies1]();
    pairEner1= new double*[nspecies1]();
    species  = new Species[nspecies1];
    for(short i=0; i<nspecies1; ++i)
    {
        pairEner1[i] = new double[nspecies1]();
    }
    //input
    
    readSystem(nspecies1,pInfo,species,sysInfo,iteraInfo,rhoB1,pairEner1,eta_t,mb,zb,fileName,filePath);
    
    ii2 = 0;
    for(short i=0; i<(nspecies1-1); ++i)
    {
        if(rhoB1[i] > 1.0E-10) ++ii2;
    }
    nspecies = ii2 + 1;
    
    ii2 = 0;
    ii3 = 0;
    for(short i=0; i<nb1; ++i)
    {
        if(rhoB1[i] > 1.0E-10) ++ii2;
        if(rhoB1[i] < 1.0E-10) mb[0][i] = 0;
    }
    for(short i=nb1; i<nblocks1; ++i)
    {
        if(rhoB1[i] > 1.0E-10) ++ii3;
        if(rhoB1[i] < 1.0E-10) mb[1][i-nb1] = 0;
    }
    nb[0] = ii2;
    nb[1] = ii3;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    hspecies= nspecies - 1;
    
    LLI  = new int[nspecies]();
    ULI  = new int[nspecies]();
    mu   = new double[hspecies]();
    rhoB = new double[nspecies]();
    rhoB_S=new double[hspecies]();
    
    etar = new double[hspecies]();
    uWall= new float[hspecies]();
    D    = new float[nspecies]();
    DSort= new short[nspecies]();
    Z    = new float[hspecies]();
    B    = new double[hspecies]();
    TB   = new double[hspecies]();
    
    
    pairEner = new double*[hspecies]();
    for(short i=0; i<hspecies; ++i)
    {
        pairEner[i] = new double[hspecies]();
    }
    
    sigma   = 0.0;
    dcharge = 0.0;
    sigma0  = 0.0;
    phi     = 0.0;
    phi0    = 0.0;
    size0   = 0.0;
    iCode   = 0;
    
    
    depthW  = 1.25;
    nm      = 1E-9;
    dr0     = sysInfo.stepLen;
    size0   = sysInfo.size;
    epsilonS= sysInfo.permiS;
    epsilonW= sysInfo.permiW;
    lunit   = sysInfo.lenth;
    errTol  = iteraInfo.errTol;
    mixCoe  = iteraInfo.mixCoe;
    maxItera= iteraInfo.maxItera;
    cstep   = iteraInfo.cstep;
    alpha   = 1.0;
    mixCoeMax0= 0.01;
    mixCoeMax1= 0.2;
    mixCoeMin0= 1E-4;
    mixCoeMin1= 5E-4;
    errTol1   = errTol*100;
    errTol2   = errTol*1000;
    errTol3   = errTol*10000;
    
    
    coe0 = kB*sysInfo.T/e0;
    coe11 = nm*nm/(lunit*lunit);// surface charged density from e0/lunit^2 to e0/nm^2
    coe1 = e0/(lunit*lunit);// surface charged density from e0/lunit^2 to C/m^2
    coe2 = lunit/nm;// length from lunit to nm
    coe3 = 1.0E-3/(lunit*lunit*lunit)/Na; //density from N/lunit^3 to Molar
    coe4 = lunit*lunit*lunit/(nm*nm*nm);
    dr   = dr0;
    
    dr   =round(dr*100000)/100000.0;
    
    phL_L   = iteraInfo.phL_L/coe0;
    //phi0    = sysInfo.phi/coe0;
    sigma0  = sysInfo.phi/coe1;

    
    f = (epsilonS - epsilonW)/(epsilonS + epsilonW);
    
    for(short i=0; i<2; ++i)
    {
        mp[i]  = pInfo[i].mp;
    }
    
    ii2 = 0;
    jj2 = 0;
    for(short i=0; i<(nspecies1-1); ++i)
    {
        if(rhoB1[i] > 1.0E-10)
        {
            jj2 = 0;
            for(short j=0; j<(nspecies1-1); ++j)
            {
                if(rhoB1[j] > 1.0E-10)
                {
                    pairEner[ii2][jj2] = pairEner1[i][j];
                    ++jj2;
                }
            }
            Z[ii2]    = species[i].Z;
            uWall[ii2]= species[i].uWall;
            
            
            ++ii2;
        }
    }
    
    
    if((ii2 != hspecies) || (jj2 != hspecies))
    {
        cerr<<"the actual # of hspecies is wrong"<<endl;
        exit(0);
    }
    
    ii2 = 0;
    for(short i=0; i<nspecies1; ++i)
    {
        if(rhoB1[i] > 1.0E-10)
        {
            D[ii2]    = species[i].D;
            rhoB[ii2] = rhoB1[i];
            ++ii2;
        }
    }
    
    if(ii2 != nspecies)
    {
        cerr<<"the actual # of nspecies is wrong"<<endl;
        exit(0);
    }
    
    for(short i=0; i<nspecies; ++i)
    {
        LLI[i]  = round(D[i]*0.5/dr);
        DSort[i]= i;
    }
    
    LLIM = LLI[0];
    for(short i=1; i<nspecies; ++i)
    {
        if(LLIM > LLI[i]) LLIM = LLI[i];
    }
    
    ZMax = 0;
    nbb  = nb[0] + nb[1];//numbers of species of polymers
    for(short i=nbb; i<hspecies; ++i)
    {
        if((fabs(Z[i])>ZMax) && (rhoB[i]>1.0E-10)) ZMax = fabs(Z[i]);
    }
    
    neutsys = 0;
    for(short i=0; i<hspecies; ++i)
    {
        rhoB_S[i]= rhoB[i];
        if((Z[i]!=0.0) && (rhoB[i]>1.0E-10)) neutsys = 1;
    }

    if((ZMax==0) && (neutsys==1))
    {
        cerr<<"the renormalization valence is wrong, should revise renormalization code; ZMax="<<ZMax<<endl;
        exit(0);
    }
    
    
    for(short i=0; i<(nspecies-1); ++i)
    {
        iSort  = DSort[i];
        SortTem= D[iSort];
        iDMin  = i;
        iSort0 = iSort;
        for(short j=(i+1); j<nspecies; ++j)
        {
            iSort = DSort[j];
            if(SortTem>D[iSort])
            {
                iDMin  = j;
                iSort0 = iSort;
                SortTem= D[iSort];
            }
        }
        DSort[iDMin]= DSort[i];
        DSort[i]    = iSort0;
        
    }
    iSort= DSort[nspecies-1];
    DMAX = round(D[iSort]/dr);
    

    //////////////////////////////////////////////////////////////////////////////////
    
    BJ  = 0.25*e0*e0/(Pi*lunit*epsilonS*epsilon0*kB*sysInfo.T);
    
    
    rhot = 0;
    eta  = 0;
    for(short i=0; i<hspecies; ++i)
    {
        rhot    = rhot + Z[i]*Z[i]*rhoB[i];
        eta     = eta + rhoB[i]*D[i]*D[i]*D[i]*Pi/6.0;
    }
    if(MethG[0] == "scft") rhoB[hspecies] = 6*(eta_t - eta)/(Pi*D[hspecies]*D[hspecies]*D[hspecies]);
    if(MethG[0] == "dft") rhoB[hspecies] = 0;
    
    
    for(short i=0; i<hspecies; ++i)
    {
        etar[i] = eta_t*(rhoB[i]*D[i]*D[i]*D[i]*Pi/6.0)/eta; //we assume the close packing is eta_t=0.74
    }
    kapa= sqrt(4*Pi*BJ*rhot);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    mp1    = mp[0];
    mp2    = mp[1];
    nblocks= nb[0] + nb[1];
    rhoBM1 = 0.0;
    temp_S = new double[nblocks]();
    for(short i=0; i<nb[0]; ++i)
    {
        rhoBM1   = rhoBM1 + rhoB[i];
        temp_S[i]= rhoB[i];
    }
    rhoBM2 = 0.0;
    for(short i=nb[0]; i<nblocks; ++i)
    {
        rhoBM2    = rhoBM2 + rhoB[i];
        temp_S[i]= rhoB[i];
    }
    
    ATT   = new double*[nblocks]();
    for(short i=0; i<nblocks; ++i)
    {
        ATT[i]  = new double[nblocks]();
    }
    /////////////////////////////////////////////////////////////////////
    
    gama  = 0;
    P_bulk= 0;
    DMAXg = DMAX;
    if(neutsys==1)
    {
        BulkGama(rhoB,D,Z,gama);
        for(short i=0; i<hspecies; ++i)
        {
            B[i]    = 0.5*(D[i] + MethG0.a/gama);
            TB[i]   = 0.5*(D[i] + MethG0.b/gama);
        }
        DMAXg = round(D[iSort]/dr) + round(MethG0.b/(dr*gama));
        DMAX1 = DMAXg*2;
        
        Psi_IJ = new double**[hspecies]();
        for(short i=0; i<hspecies; ++i)
        {
            Psi_IJ[i] = new double*[hspecies]();
            for(short j=0; j<hspecies; ++j)
            {
                Psi_IJ[i][j] = new double[3*DMAX1+1]();
            }
        }
        
        PsiChargeShellJJ(B,TB,Z,Psi_IJ);
        //PsiChargeShell(B,Z,Psi_IJ);
        //PsiChargeShellDirk(gama,D,Z,Psi_IJ);
    }
        
    //calculate bulk chemical potential
    if(MethG[0] == "scft") BulkChemPotential(rhoB,D,Z,pairEner,mu);
    if(MethG[0] == "dft") BulkChemPotentialDFT(gama,rhoB,D,B,TB,Z,pairEner,ATT,Psi_IJ,mu,P_bulk);

    ///////////////////////////////////////////////////////////////////////////////////////////
    iloop  = 0;
    //num_max= 0;
    FileOpen(iloop,cstep,vari_press,converge_rec,parameters_s,contact_theorem,fileIterative,fileILoop,
             fileTemp,fileConverge,filePath,MethG[2]);
    //Initialize the lattice
    inFile.open(fileILoop.c_str(),ios_base::in);
    if(inFile) inFile>>iloop;
    inFile.close();
    
    
    ///////////////////set the simulation box/////////////////////////
    if(MethG[2] == "single")
    {
        stiffness_b = DMAX*dr;
        if(neutsys==1) stiffness_b = stiffness_b + 1.0/gama;
        
        if(rhoBM1 > 1E-10)
        {
            R_e1        = 0;
            if(epsilon_b[0] <= 0)
            {
                R_e1 = sqrt(mp1);
            }
            else
            {
                temp = 0.0;
                temp2= exp(-2*epsilon_b[0]);
                temp1=(temp2*(1+epsilon_b[0])+epsilon_b[0]-1)/((1-temp2)*epsilon_b[0]);
                for(short i=1; i<(mp1-1); ++i)
                {
                    temp = temp + (mp1-1-i)*pow(temp1,i);
                }
                R_e1 = mp1 - 1 + 2.0*temp;
                
                R_e1 = sqrt(R_e1);
            }
            stiffness_b = R_e1;
            if(stiffness_b < 6) stiffness_b = 6;
        }
        
        
        if(rhoBM2 > 1E-10)
        {
            R_e2 = 0;
            if(epsilon_b[1] <= 0)
            {
                R_e2 = sqrt(mp2);
            }
            else
            {
                temp = 0.0;
                temp2= exp(-2*epsilon_b[1]);
                temp1=(temp2*(1+epsilon_b[1])+epsilon_b[1]-1)/((1-temp2)*epsilon_b[1]);
                for(short i=1; i<(mp2-1); ++i)
                {
                    temp = temp + (mp2-1-i)*pow(temp1,i);
                }
                R_e2 = mp2 - 1 + 2.0*temp;
                
                R_e2 = sqrt(R_e2);
            }
            if(stiffness_b < R_e2) stiffness_b = R_e2;
        }
        stiffness_b = stiffness_b*0.5;
        stiffness_b = ceil(stiffness_b);
        
        NRe     = round(stiffness_b/dr);
        ngrid   = 10*NRe;
        ngrid_m = 5*NRe;
        ngrid_b = 4*NRe;
        ngridm1 = ngrid_m + 1;
        size    = ngrid*dr;
        ///////////////////set the simulation box/////////////////////////
    }
    else if(MethG[2] == "two")
    {
        size    = size0;
        ngrid   = round(size/dr);
        ngrid_m = ngrid/2;
        ngrid_b = ngrid_m;
        ngridm1 = ngrid_m + 1;
    }
    else
    {
        cerr<<"fail to set the simulation box"<<endl;
        exit(0);
    }
    
    
    depthLR= new int*[4]();
    for(short i=0; i<4; ++i)
    {
        depthLR[i]= new int[hspecies]();
    }
    
    
    nBessel    = 2*round(1.0/dr) + 1;
    BesselZero = new double**[2]();
    for(short i=0; i<2; ++i)
    {
        BesselZero[i] = new double*[nBessel]();
        for(short j=0; j<nBessel; ++j)
        {
            BesselZero[i][j] = new double[nBessel]();
        }
    }
    
        
    Ext = new double*[hspecies]();
    rho = new double*[nspecies]();
    rho1= new double*[nspecies]();
    rho_temp1 = new double*[nspecies]();
    rho_temp2 = new double*[nspecies]();
    Psi  = new double[ngrid+1]();
    for(short i=0; i<nspecies; ++i)
    {
        if(i != hspecies) Ext[i] = new double[ngrid_m+1]();
        rho[i] = new double[ngrid+1]();
        rho1[i]= new double[ngrid+1]();
        rho_temp1[i] = new double[ngrid_m+1]();
        rho_temp2[i] = new double[ngrid_m+1]();
        ULI[i] = round((size-D[i]*0.5)/dr);
    }
    ULIM = ULI[0];
    for(short i=1; i<nspecies; ++i)
    {
        if(ULIM < ULI[i])  ULIM = ULI[i];
    }
    
    iter = 0;
    for(int k=0; k<=ngrid; ++k)
    {
        Psi[k]  = 0;
    }
        
        
    alpha   = 1.2; //cut off
    for(short i=0; i<hspecies; ++i)
    {
        depthLR[0][i] = round(D[i]*0.5/dr);
        depthLR[1][i] = round(D[i]*alpha/dr);
        depthLR[2][i] = round((size-D[i]*0.5)/dr);
        depthLR[3][i] = round((size-D[i]*alpha)/dr);
    }
    //External potential
    if(MethG[3] == "sw")
    {
        ExternalPotential(uWall,depthLR,Ext);
    }
    else if(MethG[3] == "lj")
    {
       ExternalPotential(alpha,D,uWall,depthLR,Ext);
    }
    else
    {
        cerr<<"fail to set the external potential"<<endl;
        exit(0);
    }
    
   
    
    BesselFunction(nBessel,epsilon_b,BesselZero);
    
    
    
    if(rhoBM1 > 1E-10)
    {
        DL1  = new double[mp1]();
        DR1  = new double[mp1]();
        GRC1 = new double[nBessel]();
    }
    
    if(rhoBM2 > 1E-10)
    {
        DL2  = new double[mp2]();
        DR2  = new double[mp2]();
        GRC2 = new double[nBessel]();
    }

    RenormGLGR(D,rhoB,BesselZero,MODEL);
    
    //Initialize the lattice
    sigma_t = 0;
    dcharge = 0;
    sigma_t0= 0;
    inFile.open(fileTemp.c_str(),ios_base::in);
    
    
    if(MethG[0] == "scft")
    {
        EulerLagrange(LLI,ULI,D,Z,pairEner,mu,rhoB,BesselZero,MODEL);
    }
    else if(MethG[0] == "dft")
    {
        EulerLagrangeDFT(gama,LLI,ULI,D,TB,Z,pairEner,mu,rhoB,ATT,BesselZero,Psi_IJ,MODEL);
    }
    else
    {
        cerr<<"fail to set the theoretical method"<<endl;
        exit(0);
    }

    
    
    if(!inFile)
    {
        //Initialize the system
        
        Initialization(LLI,ULI,rhoB,rho);
        //Initialization(gama,f,LLI,ULI,D,TB,etar,Z,rhoB,mu,Ext,ATT,BesselZero,Psi_IJ,rho,
        //               pairEner,MODEL,MethG);
        
        dcharge =  (sigma0-sigma_t0)/cstep;
        if(dcharge == 0) cstep = 1;
        sigma_t0 = sigma0 - dcharge*(cstep - 1);
        
    }
    else
    {
        Initialization(sigma_t0,dcharge,cstep,iter,LLI,ULI,rhoB,rho,inFile);
    }
    inFile.close();
    
    parameters_s<<"------------------------------------  system parameters  --------------------------------------------"<<endl;
    parameters_s<<"  "<<endl;
    parameters_s<<"------------------------------  parameters for simulation box  --------------------------------------"<<endl;
    parameters_s<<"number of surfaces: "<<MethG[2]<<";  size of box: "<<size<<" [unit];  unit: "<<lunit<<" [m]"<<endl;
    if(rhoBM1 > 1E-10)
    {
        parameters_s<<"   "<<endl;
        parameters_s<<"---------------------------------  parameters for polymer 1  ----------------------------------------"<<endl;
        parameters_s<<"monomer concentrations: "<<rhoBM1<<" [1/unit^3] or "<<(rhoBM1*coe3)<<" [M];  number of blocks: "
                    <<nb[0]<<";  polymerization: "<<mp[0]<<endl;
    }
    if(rhoBM2 > 1E-10)
    {
        parameters_s<<"   "<<endl;
        parameters_s<<"---------------------------------  parameters for polymer 2  ----------------------------------------"<<endl;
        parameters_s<<"monomer concentration: "<<rhoBM2<<" [1/unit^3] or "<<(rhoBM2*coe3)<<" [M];  number of blocks: "
                    <<nb[1]<<";  polymerization: "<<mp[1]<<endl;
    }
    if((species[nblocks1].Z !=0)&&((rhoB1[nblocks1]+rhoB1[nblocks1+1])>1E-10))
    {
        parameters_s<<"   "<<endl;
        parameters_s<<"-----------------------------------  parameters for salts  ------------------------------------------"<<endl;
        parameters_s<<"salt concentrations: "<<(rhoB1[nblocks1]+rhoB1[nblocks1+1])<<" [1/unit^3] or "
                    <<(coe3*rhoB1[nblocks1]/species[nblocks1].Z)<<" [M];  Valencies: "<<species[nblocks1].Z
                    <<"  :  "<<species[nblocks1+1].Z<<endl;
    }
    if(neutsys==1)
    {
        parameters_s<<"   "<<endl;
        parameters_s<<"----------------------------------  parameters for others  ------------------------------------------"<<endl;
        parameters_s<<"surface charge density: "<<(sigma0*coe1)<<" [C/m^2] or "<<(coe11*sigma0)<<" [e/nm^2]"<<endl;
        parameters_s<<"Bjerrum length: "<<BJ<<" [unit]"<<";  Deybe Screening length: "<<(1.0/kapa)<<" [unit]"<<endl;;
    }
    parameters_s<<"  "<<endl;
    
    //loop for changing the surface charge density
    for(int i_num=iloop; i_num<cstep; ++i_num)//loop for changing the surface charge density
    {
        //errTol0 = errTol*10;
        errTol0 = errTol;
        //if(i_num == cstep) errTol0 = errTol;
        //phi = dcharge*i_num;
        sigma = sigma_t0 + dcharge*i_num;
        boundary= -4*Pi*BJ*sigma*dr;
        
        iterative_out.open(fileIterative.c_str());
        iterative_out.setf(ios::fixed);
        iterative_out.precision(8);
        if(!iterative_out.is_open())
        {
            cerr<<"fail to open the file for iterative"<<endl;
            exit(0);
        }
        
        
        //err       = 1.0;
        mix_f      = mixCoe;
        iback      = 0;
        start_back = clock();
        err    = 0;
        
        thresh_in  = 0;
        thresh_out = 0;
        thresh_nor = 0;
        thresh_copy= 0;
        err_temp1  = 0;
        std_temp   = 1;
        
        AccelKeep  = 0;
        deltaPhi   = 0;
        C_psi      = 1;
        
        
        do
        {
            //EulerLagrange(phL_L,f,eta_t,LLI,ULI,D,Z,etar,Psi,pairEner,mu,rhoB,Ext,BesselZero,
            //rho,rho1,MODEL);
            if(MethG[0] == "scft")
            {
                EulerLagrange(sigma,f,eta_t,deltaPhi,err,LLI,ULI,D,etar,Z,Psi,pairEner,mu,rhoB,Ext,
                              BesselZero,rho,rho1,MODEL);
            }
            
            if(MethG[0] == "dft")
            {
                EulerLagrangeDFT(gama,sigma,f,eta_t,deltaPhi,err,LLI,ULI,D,TB,etar,Z,Psi,pairEner,mu,rhoB,Ext,
                                 ATT,BesselZero,Psi_IJ,rho,rho1,MODEL);
            }
            /////////////////////////////////////////////////////////////////////////////////////
            fsigma = 0;
            for(short i=0; i<hspecies; ++i)
            {
                if(Z[i] != 0)
                {
                    temp    = SimpsonIntegration(rho1[i],0,ngridm1,0,ngrid_b);
                    fsigma += temp*Z[i];
                }
            }
            
            sigma_t = -fsigma;
            //////////////////////////////////////////////////////////////////////////////////////
            
            
            err = 0;
            for(short i=0; i<hspecies; ++i)
            {
                if(rhoB[i] > 1E-10)
                {
                    for(int k=LLI[i]; k<=ngrid_b; ++k)
                    {
                        temp = fabs(rho1[i][k]-rho[i][k])/rhoB[i];
                        if(err < temp) err = temp;
                    }
                }
            }
            
            //improve the convergence: how to choose the mixture parameter-----start
            if(err >= errTol)
            {
                ++iter;
                
                
                if(thresh_copy > 1)
                {
                    //test = (err_temp1-err)/err_temp1;
                    if(err <= err_temp1)
                    {
                        //err_std = (err_temp1-err)/err_temp1;
                        err_std = (err_temp1-err)/err;
                        //err_std = err_temp1-err;
                        
                        if((err<errTol1) && (err_std>std_temp))
                        {
                            if(AccelKeep == 0) mix_f += 5*mixCoeMin0;
                            
                            ++AccelKeep;
                            if(AccelKeep > 5) AccelKeep = 0;
                            
                        }
                        else if((err<errTol2) && (err_std>std_temp))
                        {
                            if(AccelKeep == 0) mix_f += 4*mixCoeMin0;
                            
                            ++AccelKeep;
                            if(AccelKeep > 5) AccelKeep = 0;
                            
                        }
                        else if((err<errTol3) && (err_std>std_temp))
                        {
                            if(AccelKeep == 0) mix_f += 3*mixCoeMin0;
                            
                            ++AccelKeep;
                            if(AccelKeep > 5) AccelKeep = 0;
                            
                        }
                        else if(err_std > std_temp)
                        {
                            if(AccelKeep == 0) mix_f += 2*mixCoeMin0;
                            
                            ++AccelKeep;
                            if(AccelKeep > 5) AccelKeep = 0;
                            
                        }
                        else
                        {
                            if(AccelKeep == 0) mix_f += mixCoeMin0;
                            
                            ++AccelKeep;
                            if(AccelKeep > 5) AccelKeep = 0;
                        }
                        
                        ++thresh_in;
                        thresh_nor = 0;
                        thresh_out = 0;
                        std_temp   = err_std;
                        
                        if((mix_f > mixCoeMax1) && (err<errTol1)) mix_f = mixCoeMax1;
                        if((mix_f > mixCoeMax0) && (err>errTol1)) mix_f = mixCoeMax0;
                        //if(mix_f < mixCoeMin1) mix_f = mixCoeMin1;
                    }
                    else if((thresh_nor==0) && (thresh_in>2 || thresh_out>0))
                    {
                        if(thresh_out == 0)
                        {
                            for(short i=0; i<nspecies; ++i)
                            {
                                for(int k=LLI[i]; k<=ngrid_m; ++k)
                                {
                                    rho[i][k] = rho_temp2[i][k];
                                    rho[i][ngrid-k] = rho[i][k];
                                }
                            }
                            mix_f = mixCoe;
                            err_temp1 = 0;
                        }
                        
                        thresh_nor = 0;
                        if(thresh_out > 50) thresh_nor = 1;
                        if(thresh_out > 5) err_temp1 = err;
                        ++thresh_out;
                        thresh_in  = 0;
                        std_temp   = 1;
                        AccelKeep  = 0;
                        
                    }
                    else
                    {
                        mix_f  = mix_f - 2*mixCoeMin0;
                        mix_f0 = mixCoe/err;
                        if(mix_f0 < mix_f) mix_f = mix_f0;
                        if(mix_f  > mixCoeMax0) mix_f = mixCoeMax0;
                        if(mix_f  < mixCoeMin1) mix_f = mixCoeMin1;
                        std_temp   = 1;
                        thresh_out = 0;
                        thresh_in  = 0;
                        AccelKeep  = 0;
                    }
                    
                }
                else
                {
                    mix_f = mixCoe/err;
                    if(mix_f > mixCoeMax0) mix_f = mixCoeMax0;
                    if(mix_f < mixCoeMin1) mix_f = mixCoeMin1;
                    std_temp  = 1;
                    thresh_out= 0;
                    AccelKeep = 0;
                }
                
                //back up
                if(thresh_out ==0)
                {
                    for(short i=0; i<nspecies; ++i)
                    {
                        for(int k=LLI[i]; k<=ngrid_m; ++k)
                        {
                            rho_temp2[i][k] = rho_temp1[i][k];
                            rho_temp1[i][k] = rho[i][k];
                        }
                    }
                    
                    //err_temp2 = err_temp1;
                    err_temp1 = err;
                    
                    ++thresh_copy;
                }
                //improve the convergence: how to choose the mixture parameter-----end
                
                for(short i=0; i<nspecies; ++i)
                {
                    for(int k=LLI[i]; k<=ngrid_b; ++k)
                    {
                        rho[i][k] = (1-mix_f)*rho[i][k] + mix_f*rho1[i][k];
                        rho[i][ngrid-k] = rho[i][k];
                    }
                }
            }
        
            iterative_out<<"iter=  "<<iter<<" err= "<<err<<" mix_f="<<mix_f<<" deltaPhi= "<<(deltaPhi*coe0)<<" sigma_t="<<(sigma_t*coe1)<<endl;
            

            //if(iter==2) exit(0);
            //back_up
            finish_back = clock();
            time_backup = (double)(finish_back-start_back)/CLOCKS_PER_SEC;
            time_backup = time_backup/3600.0;
            //if(((time_backup-iback*0.25) > 0.25) || (time_backup>3.98) || (iter == maxItera))
            
            if(((time_backup-iback*1) > 1) || (iter == maxItera))
            {
                FileOutPut(sigma_t0,dcharge,cstep,iCode,iter,i_num,Z,rho,fileILoop,fileTemp);
                ++iback;
            }
            
        }
        while((err >= errTol0) && (iter <= maxItera));
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        if(MethG[2] == "single")
        {
            double C_rho;
            C_rho = 0;
            for(short i=0; i<hspecies; ++i)
            {
                C_rho += rho[i][LLI[i]];
            }
            
            QQ = sigma0*sigma0*2*Pi*BJ;
            contact_theorem<<"P_bulk= "<<P_bulk<<" Contact="<<(C_rho - QQ)<<endl;
        }
        
        
        iterative_out.unsetf(ios::fixed);
        iterative_out.close();
        
        //if((err > errTol0) && (iter > maxItera))
        if(err > errTol0)
        {
            iter = 0;
            converge_rec<<"fail to get convergence result:  size= "<<size<<" err= "<<err<<endl;

            FileOutPut(i_num,lunit,coe3,rhoB,rho,Psi,filePath); //i_num
            FileOutPut(sigma_t0,dcharge,cstep,iCode,0,i_num+1,Z,rho,fileILoop,fileTemp);
        }
        else
        {
            iter = 0;
            
            ///////////////////////////////////////////////////////////////////////////////////////
            //output
            FileOutPut(i_num,lunit,coe3,rhoB,rho,Psi,filePath); //i_num
            
            FileOutPut(sigma_t0,dcharge,cstep,iCode,0,i_num+1,Z,rho,fileILoop,fileTemp);
            
            if(MethG[2] == "two")
            {
                if(MethG[0] == "scft") EnergyCalculation(sigma,f,eta_t,LLI,ULI,D,Z,rhoB,pairEner,rho,Psi,Ener_tot);
                    
                if(MethG[0] == "dft") EnergyCalculation(sigma,gama,f,LLI,ULI,D,Z,rhoB,TB,pairEner,ATT,rho,Psi,Psi_IJ,Ener_tot);
                
                vari_press<<size<<" "<<sigma<<" "<<Ener_tot<<endl;
            }
            
        }
        
        //break;
    }//loop for changing the surface charge density
        
        
        
        for(short i=0; i<nspecies; ++i)
        {
            if(i != hspecies) delete [] Ext[i];
            delete [] rho[i];
            delete [] rho1[i];
            delete [] rho_temp1[i];
            delete [] rho_temp2[i];
        }
        delete [] rho;
        delete [] rho1;
        delete [] rho_temp1;
        delete [] rho_temp2;
        delete [] Psi;
        delete [] Ext;
        delete [] temp_S;
    
        
        for(short i=0; i<4; ++i)
        {
            delete [] depthLR[i];
        }
        delete [] depthLR;

    //}//loop for changing the surface condition
    FileClose(vari_press,converge_rec,parameters_s,contact_theorem,MethG[2]);
   

    ifstream file_r;
    int ch_file,ch_EOF;
    file_r.open(fileConverge.c_str(),ios::in);
    ch_file=file_r.get();
    ch_EOF =EOF;
    file_r.close();
    if((ch_file == ch_EOF) && (remove(fileConverge.c_str()) != 0))
    {
        cerr<<"Remove operation failed"<<endl;
        exit(0);
    }
    
    
    
    
    delete [] LLI;
    delete [] ULI;
    delete [] rhoB;
    delete [] rhoB_S;
    delete [] rhoB1;
    delete [] etar;
    delete [] uWall;
    delete [] D;
    delete [] DSort;
    delete [] B;
    delete [] TB;
    delete [] Z;
    delete [] mu;
    delete [] species;
    delete [] pInfo;
    delete [] MODEL;
    delete [] MethG;
    delete [] mp;
    delete [] nb;
    delete [] epsilon_b;
    //delete [] Ener_Save;
    
    delete [] mb[0];
    delete [] mb[1];
    delete [] mb;
    
    delete [] zb[0];
    delete [] zb[1];
    delete [] zb;
    
    for(short i=0; i<nblocks; ++i)
    {
        delete [] ATT[i];
    }
    delete [] ATT;

    for(short i=0; i<hspecies; ++i)
    {
        delete [] pairEner[i];
    }
    delete [] pairEner;
    
    for(short i=0; i<nspecies1; ++i)
    {
        delete [] pairEner1[i];
    }
    delete [] pairEner1;
    
    
    for(short i=0; i<2; ++i)
    {
        for(short j=0; j<nBessel; ++j)
        {
            delete [] BesselZero[i][j];
        }
        delete [] BesselZero[i];
    }
    delete [] BesselZero;
    
    if(neutsys==1)
    {
        for(short i=0; i<hspecies; ++i)
        {
            for(short j=0; j<hspecies; ++j)
            {
                delete [] Psi_IJ[i][j];
            }
            delete [] Psi_IJ[i];
        }
        delete [] Psi_IJ;
    }
    
    if(rhoBM1 > 1E-10)
    {
        delete [] DL1;
        delete [] DR1;
        delete [] GRC1;
    }
    
    if(rhoBM2 > 1E-10)
    {
        delete [] DL2;
        delete [] DR2;
        delete [] GRC2;
    }
    
    
    return 0;
}



