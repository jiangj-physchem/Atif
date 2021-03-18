//define systemset.h
#ifndef SYSTEMSET_H_
#define SYSTEMSET_H_
#include "clibrary.h"
#include "constantnum.h"
#include <ctype.h>

extern short  nspecies;
extern short* nb;//# of blocks in copolymer i;
using namespace std;


namespace systemset
{
    struct Species        //coordinate for short
    {
        float Z;    //valency
        float D;    //diameter
        float uWall;//interaction with wall; unit: kBT
        
        Species()
        {
            Z = 0.0;
            D = 0.0;
            uWall = 0.0;
        }
    };
    struct SysInfo //dielectric constant assigned on the field link
    {
        float size;  //separation of two plates
        double phi; //the surface electric potential
        float T;     //Temperature
        float permiS; //dilectric constant for solution
        float permiW; //dilectric constant for surfaces
        float stepLen; //step length of the numerical method
        double lenth; //length unit
        //short iCode;//initialize
        
        SysInfo()
        {
            size = 0.0;
            phi  = 0.0;
            T = 0.0;
            permiS = 0.0;
            permiW = 0.0;
            stepLen = 0.0;
            lenth = 0.0;
        }
    };
    struct IteraInfo
    {
        int    cstep;
        double mixCoe;  //mixture coefficient of picard iteration
        double phL_L;
        double errTol; //the error tolerence
        double maxItera; //maximum iteration steps
        
        IteraInfo()
        {
            cstep  = 0;
            mixCoe = 0.0;
            phL_L  = 0.0;
            errTol = 0.0;
            maxItera = 0.0;
        }
    };
    struct PolyInfo
    {
        short mp;   //polymerization
        short nb;   //number of blocks
        double rhopm;//monomer density
        string MODEL;
        double epsilon_b;
        
        PolyInfo()
        {
            mp   = 0;
            nb   = 0;
            rhopm= 0;
            MODEL= "";
            epsilon_b=0;
        }
    };
    
    struct MethodGeom
    {
        string Method;
        string Geomet;
        string nSurfs;
        string exPoten;
        double a;
        double b;
        
        MethodGeom()
        {
            Method= "";
            Geomet= "";
            nSurfs= "";
            exPoten= "";
            a = 1;
            b = 0;
        }
    };
    
    void readSystem(short& nspecies1,PolyInfo* pInfo,Species* species,SysInfo& sysInfo,IteraInfo& iteraInfo,
                    double* rhoB,double** pairEner,double& eta0,short** mb,float** zb,char fileName[],string& filePath);
    
    void readSystem(PolyInfo* pInfo,char fileName[],MethodGeom& MethG);
}
#endif // SYSTEMSET_H_
