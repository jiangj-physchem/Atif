//***********solve the possion equation****************//
#include "fileoutput.h"
//#include "simpsonintegration.h"
//#include "constantnum.h"

extern double dr;
extern int ngrid;
extern int ngrid_m; //the number of grids: the middle between two surfaces
extern short  nspecies;

void FileOpen(int num_min,int num_max,ofstream& vari_press,ofstream& converge_rec,ofstream& inum_sigma,
              ofstream& contact_theorem,string& fileIterative,string& fileILoop,string& fileTempDen,
              string& filePath)
{
    short  iDirect;
    string fileTemp;
    string fileDash  = "_";
    string fileVar_Press= "variable_pressure_";
    string fileNum_Sigma= "inum_sigma_";
    string fileContT= "contact_theorem_";
    string fileConverge = "converge_record_";
    string suffix=".dat";
    ostringstream os1;
    ostringstream os2;
    
    
    os1<<num_min;
    os2<<num_max;
    
    fileIterative = fileIterative + os1.str() + fileDash + os2.str() + suffix;
    fileILoop     = fileILoop + os1.str() + fileDash + os2.str() + suffix;
    fileTempDen   = fileTempDen + os1.str() + fileDash + os2.str() + suffix;
    fileVar_Press = fileVar_Press + os1.str() + fileDash + os2.str() + suffix;
    fileConverge  = fileConverge + os1.str() + fileDash + os2.str() + suffix;
    fileNum_Sigma = fileNum_Sigma + os1.str() + fileDash + os2.str() + suffix;
    fileContT     = fileContT + os1.str() + fileDash + os2.str() + suffix;
    
    iDirect    = 0;
    fileTemp   = "";
    filePath   = "../" + filePath;
    fileTemp   = filePath + fileVar_Press;
    vari_press.open(fileTemp.c_str(),ios::app);
    //vari_press.open(fileTemp.c_str());
    while(!vari_press.is_open() && iDirect < 10)
    {
        filePath = "../" + filePath;
        fileTemp = filePath + fileVar_Press;
        vari_press.open(fileTemp.c_str(),ios::app);
        ++iDirect;
    }
    if(iDirect < 10)
    {
        vari_press.setf(ios::fixed);
        vari_press.precision(8);
    }
    else
    {
        cerr<<"fail to open the file for variable vs pressure file"<<endl;
        exit(0);
    }
    
    fileTemp   = filePath + fileConverge;
    converge_rec.open(fileTemp.c_str(),ios::app);
    converge_rec.setf(ios::fixed);
    converge_rec.precision(8);
    if(!converge_rec.is_open())
    {
        cerr<<"fail to open the file for converge record"<<endl;
        exit(0);
    }
    
    
    fileTemp   = filePath + fileContT;
    contact_theorem.open(fileTemp.c_str(),ios::app);
    contact_theorem.setf(ios::fixed);
    contact_theorem.precision(8);
    if(!contact_theorem.is_open())
    {
        cerr<<"fail to open the file for contact theorem record"<<endl;
        exit(0);
    }
    
    fileTemp   = filePath + fileNum_Sigma;
    inum_sigma.open(fileTemp.c_str(),ios::app);
    inum_sigma.setf(ios::fixed);
    inum_sigma.precision(8);
    if(!inum_sigma.is_open())
    {
        cerr<<"fail to open the file for inum_sigma record"<<endl;
        exit(0);
    }
    
    
    os1.str("");
    os2.str("");
    os1.clear();
    os2.clear();
    
}

void FileClose(ofstream& vari_press,ofstream& converge_rec,ofstream& inum_sigma,ofstream& contact_theorem)
{
    vari_press.setf(ios::fixed);
    vari_press.precision(8);
    
    converge_rec.setf(ios::fixed);
    converge_rec.precision(8);
    
    inum_sigma.setf(ios::fixed);
    inum_sigma.precision(8);
    
    contact_theorem.setf(ios::fixed);
    contact_theorem.precision(8);
}




void FileOutPut(int num_file,double lunit,double coe3,double* rhoB,double** rho,double* Psi,string filePath)
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    double  R,nm;
    double* rho_put;
    
    ostringstream os;
    string fileDensity  = "density_potential_";
    string suffix=".dat";
    ofstream density_potential;
    
    rho_put = new double[nspecies]();
    nm  = 1E-9;
    
    os<<num_file;
    
    fileDensity= fileDensity + os.str() + suffix;
    fileDensity= filePath + fileDensity;
    density_potential.open(fileDensity.c_str());
    density_potential.setf(ios::fixed);
    density_potential.precision(8);
    if(!density_potential.is_open())
    {
        cerr<<"fail to open the file for density profile and electric potential"<<endl;
        exit(0);
    }
    os.str("");
    os.clear();
    
    
    
    for(int k=0; k<=ngrid_m; ++k)
    {
        for(short i=0; i<nspecies; ++i)
        {
            if(rhoB[i] > 1E-10)
            {
                //rho_put[i] = rho[i][k]/rhoB[i];
                rho_put[i] = rho[i][k]*coe3;
            }
            else
            {
                rho_put[i] = 0;
            }
        }
        
        R = k*dr;
        
        density_potential<<(R*lunit/nm)<<" ";
        for(short i=0; i<nspecies; ++i)
        {
            density_potential<<rho_put[i]<<" ";
        }
        density_potential<<Psi[k]<<endl;
    }
    
    density_potential.unsetf(ios::fixed);
    density_potential.close();
    
    
    
    delete [] rho_put;

}



void FileOutPut(double& sigma_t,double& dcharge,int& cstep,short& iCode,int iter,
                float* Z,double** rho,ofstream& density_temp,string& fileTemp)
{
    //double ngridm1;
    //double fsigma,temp;
    //fsigma  = 0;
    //ngridm1 = ngrid_m + 1;
    
    //for(short i=0; i<nspecies; ++i)
    //{
        //if(Z[i] != 0)
        //{
            //temp    = SimpsonIntegration(rho[i],0,ngridm1,0,ngrid_m);
            //fsigma += temp*Z[i];
        //}
    //}
    //sigma_t = -fsigma;
    
    density_temp.open(fileTemp.c_str());
    density_temp.setf(ios::fixed);
    density_temp.precision(12);
    if(density_temp.is_open())
    {
        iCode = 1;
        density_temp<<sigma_t<<" "<<dcharge<<" "<<cstep<<" "<<iter<<" "<<ngrid<<endl;
        for(int k=0; k<=ngrid_m; ++k)
        {
            for(short i=0; i<(nspecies-1); ++i)
            {
                density_temp<<rho[i][k]<<" ";
            }
            density_temp<<rho[nspecies-1][k]<<endl;
        }
        
    }
    else
    {
        cerr<<"fail to open the file for temporary density profile"<<endl;
        exit(0);
    }
    density_temp.unsetf(ios::fixed);
    density_temp.close();
    
}



void FileOutPut(double* Psi,double** rho,ofstream& density_temp,string& fileTemp)
{

    
    density_temp.open(fileTemp.c_str());
    density_temp.setf(ios::fixed);
    density_temp.precision(12);
    if(density_temp.is_open())
    {
        for(int k=0; k<=ngrid_m; ++k)
        {
            density_temp<<(k*dr)<<" ";
            for(short i=0; i<nspecies; ++i)
            {
                density_temp<<rho[i][k]<<" ";
            }
            density_temp<<Psi[k]<<endl;
        }
        
    }
    else
    {
        cerr<<"fail to open the file for temporary density profile"<<endl;
        exit(0);
    }
    density_temp.unsetf(ios::fixed);
    density_temp.close();
    
}


void FileOutPut(double& sigma_t,double& dcharge,int& cstep,short& iCode,int iter,int i_num,float* Z,
                double** rho,ofstream& iloop_temp,string& fileILoop,ofstream& density_temp,string& fileTemp)
{
    iloop_temp.open(fileILoop.c_str());
    if(iloop_temp.is_open())
    {
        iloop_temp<<i_num<<endl;
    }
    else
    {
        cerr<<"fail to open the file for temporary iloop"<<endl;
        exit(0);
    }
    iloop_temp.close();
    
    density_temp.open(fileTemp.c_str());
    density_temp.setf(ios::fixed);
    density_temp.precision(12);
    if(density_temp.is_open())
    {
        iCode = 1;
        density_temp<<sigma_t<<" "<<dcharge<<" "<<cstep<<" "<<iter<<" "<<ngrid<<endl;
        for(int k=0; k<=ngrid_m; ++k)
        {
            for(short i=0; i<(nspecies-1); ++i)
            {
                density_temp<<rho[i][k]<<" ";
            }
            density_temp<<rho[nspecies-1][k]<<endl;
        }
        
    }
    else
    {
        cerr<<"fail to open the file for temporary density profile"<<endl;
        exit(0);
    }
    density_temp.unsetf(ios::fixed);
    density_temp.close();
    
}
