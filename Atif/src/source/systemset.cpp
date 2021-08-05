//**********************Read the system parameters********************//
#include "systemset.h"

namespace systemset
{
    void readSystem(short& nspecies1,PolyInfo* pInfo,Species* species,SysInfo& sysInfo,IteraInfo& iteraInfo,
                    double* rhoB,double** pairEner,double& eta0,short** mb,float** zb,char fileName[],string& filePath)
    {
        ios::sync_with_stdio(false);
        cin.tie(0);
        string s;
        string mark[12];
        //string markPath = "../program/src/mark.dat";
        char* temp2;
        short* i22;
        short i,j,k,n,i2,i4;
        short len;
        short ilen;
        float ZsaltPo,ZsaltNe;
        double rhoCounP,rhoCounN,rhoTemp,rhoSalt,eta,D3,lunit,coe3;
        int    nChange0,nChange1;
        short  nbb,iDent4,iDent5; //iDent4,iDent5,
        ifstream inFile;
        ostringstream os;

        i4  = 0;
        nbb = nb[0] + nb[1];
        
        
        i22 = new short[2]();

        mark[0] = "METHOD:";
        mark[1] = "POLYMER:";
        mark[2] = "SEQUENCE:";
        mark[3] = "SIZE:";
        mark[4] = "SALT_HS:";
        mark[5] = "WALL:";
        mark[6] = "ENERGY:";
        mark[7] = "VALENCY:";
        mark[8] = "DIAMETER:";
        mark[9] = "PERMITEMLEN:";
        mark[10] = "ITERATIVE:";
        mark[11] = "FILEPATH:";

        inFile.open(fileName,ios_base::in);
        if(inFile)
        {
            while(getline(inFile,s))
            {
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[2].c_str())==0)
                {
                    for(short i0=0; i0<2; ++i0)
                    {
                        i22[i0] = nb[i0];
                        if(getline(inFile,s))
                        {
                            ilen = 0;
                            //strcpy(temp1,s.c_str());
                            len   = s.length();
                            temp2 = new char[len+1] ();
                            i=0;
                            j=0;
                            k=0;
                            n=1;
                            do{
                                if(k!=0)  ++ilen;
                                if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                                {
                                    temp2[i] = s[ilen];
                                    ++i;
                                    n=0;
                                }
                                
                                if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                    || k==(len-1))&&(n==0))
                                {
                                    if(j < nb[i0])
                                    {
                                       mb[i0][j] = atoi(temp2);
                                    }
                                    else
                                    {
                                        cout<<"Warning! mb[i0][j]; the # of blocks for polymer "<<(i0+1)<<": the Column number is larger than nb[i0], the input is cut from nb[i0]. Please check it"<<endl;
                                        break;
                                    }
                                    
                                    memset(temp2,0,(i+1)*sizeof(char));
                                    ++j;
                                    i=0;
                                    n=1;
                                }
                                ++k;
                            }while(k < len);
                            i22[i0] = j;
                            j=0;
                            k=0;
                            
                            delete [] temp2;
                        }
                        else
                        {
                            cerr<<"Error in systemset.cpp: fail to read PolymerInfo:"<<endl;
                            exit(0);
                        }
                        
                        if(i22[i0] < nb[i0])
                        {
                            cout<<"Warning! mb[i0][j]; the # of blocks for polymer "<<(i0+1)<<": i22[i0]="<<i22[i0]<<" < nb[i0]. The rest Parameters are automatic completion. Please check it"<<endl;
                            i2 = i22[i0];
                            for(short i3=i2; i3<nb[i0]; ++i3)
                            {
                                mb[i0][i3] = mb[i0][i2-1];
                            }
                        }

                    }
                    
                    
                    for(short i0=0; i0<2; ++i0)
                    {
                        i22[i0] = nb[i0];
                        if(getline(inFile,s))
                        {
                            ilen = 0;
                            //strcpy(temp1,s.c_str());
                            len   = s.length();
                            temp2 = new char[len+1] ();
                            i=0;
                            j=0;
                            k=0;
                            n=1;
                            do{
                                if(k!=0)  ++ilen;
                                if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                                {
                                    temp2[i] = s[ilen];
                                    ++i;
                                    n=0;
                                }
                                
                                if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                    || k==(len-1))&&(n==0))
                                {
                                    if(j < nb[i0])
                                    {
                                        zb[i0][j] = atof(temp2);
                                        nChange0 = pow(10,i-2);
                                        if(nChange0 < 1) nChange0 = 1;
                                        nChange1 = round(zb[i0][j]*nChange0);
                                        zb[i0][j] = ((double) nChange1)/nChange0;
                                    }
                                    else
                                    {
                                        cout<<"Warning! zb[i0][j]; the # of blocks for polymer "<<(i0+1)<<": the Column number is larger than nb[i0], the input is cut from nb[i0]. Please check it"<<endl;
                                        break;
                                    }

                                    
                                    memset(temp2,0,(i+1)*sizeof(char));
                                    ++j;
                                    i=0;
                                    n=1;
                                }
                                ++k;
                            }while(k < len);
                            i22[i0] = j;
                            j=0;
                            k=0;
                            
                            delete [] temp2;
                        }
                        else
                        {
                            cerr<<"Error in systemset.cpp: fail to read PolymerInfo:"<<endl;
                            exit(0);
                        }
                        
                        if(i22[i0] < nb[i0])
                        {
                            cout<<"Warning! zb[i0][j]; the # of blocks for polymer "<<(i0+1)<<": i22[i0]="<<i22[i0]<<" < nb[i0]. The rest Parameters are automatic completion. Please check it"<<endl;
                            i2 = i22[i0];
                            for(short i3=i2; i3<nb[i0]; ++i3)
                            {
                                zb[i0][i3] = zb[i0][i2-1];
                            }
                        }
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[3].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }
                            
                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j==0) //sysInfo.size = atof(temp2);
                                {
                                    sysInfo.size = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(sysInfo.size*nChange0);
                                    sysInfo.size = ((double) nChange1)/nChange0;
                                }
                                if(j==1)
                                {
                                    sysInfo.stepLen = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(sysInfo.stepLen*nChange0);
                                    sysInfo.stepLen = ((double) nChange1)/nChange0;
                                }
                                
                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        delete [] temp2;
                    }
                    else
                    {
                        cerr<<"Error in systemset.cpp: fail to read Size:"<<endl;
                        exit(0);
                    }
                }
                
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[4].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';'
                               && s[ilen] != ':')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }

                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || s[ilen] == ':' || k==(len-1))&&(n==0))
                            {
                                if(j==0)// rhoSalt  = atof(temp2);
                                {
                                    rhoSalt = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(rhoSalt*nChange0);
                                    rhoSalt  = ((double) nChange1)/nChange0;
                                }
                                if(j==1)// rhoB[nspecies1-1]  = atof(temp2);
                                {
                                    eta0 = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(eta0*nChange0);
                                    eta0 = ((double) nChange1)/nChange0;
                                }

                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        delete [] temp2;
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read Salt:"<<endl;
                       exit(0);
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[5].c_str())==0)
                {
                    i4 = 0;
                    i2 = nspecies1 + 1;
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }

                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j==0)
                                {
                                    //sysInfo.phi = atof(temp2);
                                    sysInfo.phi = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(sysInfo.phi*nChange0);
                                    sysInfo.phi = ((double) nChange1)/nChange0;
                                }
                                else if(j <= nspecies1)
                                {
                                    //species[j-1].uWall = atof(temp2);
                                    species[j-1].uWall = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(species[j-1].uWall*nChange0);
                                    species[j-1].uWall = ((double) nChange1)/nChange0;
                                }
                                else
                                {
                                    if(i4 == 0)
                                    {
                                        cout<<"Warning! the external potential: the Column number is larger than nspecies1, the input is cut from nspecies1. Please check it"<<endl;

                                    }
                                    i4 = 1;
                                    break;
                                }

                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        i2 = j;
                        j=0;
                        k=0;
                        delete [] temp2;
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read Wall energy"<<endl;
                       exit(0);
                    }
                    
                    if(i2 <= nspecies1)
                    {
                        cout<<"Warning! the external potential: i2="<<i2<<" < nspecies1. The rest external potential Parameters are automatic completion. Please check it"<<endl;
                        for(short i3=(i2-1); i3<nspecies1; ++i3)
                        {
                            species[i3].uWall = species[i2-2].uWall;
                        }
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[6].c_str())==0)
                {
                    i4   = 0;
                    i2   = nspecies1;
                    for(short i0=0; i0<nspecies1; ++i0)
                    {
                        if(getline(inFile,s))
                        {
                            ilen = 0;
                            //strcpy(temp1,s.c_str());
                            len   = s.length();
                            temp2 = new char[len+1] ();
                            i=0;
                            j=0;
                            k=0;
                            n=1;
                            do{
                                if(k!=0)  ++ilen;
                                if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                                {
                                    temp2[i] = s[ilen];
                                    ++i;
                                    n=0;
                                }
                                
                                if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                    || k==(len-1))&&(n==0))
                                {
                                    if(j < nspecies1)
                                    {
                                        //pairEner[i0][j] = atof(temp2);
                                        pairEner[i0][j] = atof(temp2);
                                        nChange0 = pow(10,i-2);
                                        if(nChange0 < 1) nChange0 = 1;
                                        nChange1 = round(pairEner[i0][j]*nChange0);
                                        pairEner[i0][j] = ((double) nChange1)/nChange0;
                                    }
                                    else
                                    {
                                        if(i4 == 0)
                                        {
                                            cout<<"Warning! the pairEner input: the Column number is larger than nspecies1, the input is cut from nspecies1. Please check it"<<endl;
                                        }
                                        i4 = 1;
                                        break;
                                    }
                                    memset(temp2,0,(i+1)*sizeof(char));
                                    ++j;
                                    i=0;
                                    n=1;
                                }
                                ++k;
                            }while(k < len);
                            if(i0 == 0) i2 = j;
                            j=0;
                            k=0;
                            delete [] temp2;
                        }
                        else
                        {
                            cerr<<"Error in systemset.cpp: fail to read Pairwise Energy:"<<endl;
                            exit(0);
                        }
                        if((i2 < nspecies1)&&(i0 >= i2)) break;
                    }
                    
                    if(i2 < nspecies1)
                    {
                        cout<<"Warning! the pairEner input: the actual Row number i2="<<i2<<" < nspecies1. The rest Energy Parameters are automatic completion. Please check it"<<endl;
                        i=0;
                        for(short i3=i2; i3<nspecies1; ++i3)
                        {
                            ++i;
                            j=0;
                            for(short j3=i2; j3<nspecies1; ++j3)
                            {
                                ++j;
                                pairEner[i3][j3] = pairEner[i3-i][j3-j];
                            }
                        }
                    }
                    
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[7].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';'
                               && s[ilen] != ':')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }

                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || s[ilen] == ':' || k==(len-1))&&(n==0))
                            {
                                //species[nbb+j].Z = atof(temp2);
                                species[nbb+j].Z = atof(temp2);
                                nChange0 = pow(10,i-2);
                                if(nChange0 < 1) nChange0 = 1;
                                nChange1 = round(species[nbb+j].Z*nChange0);
                                species[nbb+j].Z = ((double) nChange1)/nChange0;
                                
                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        
                        delete [] temp2;
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read Valency:"<<endl;
                       exit(0);
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[8].c_str())==0)
                {
                    i2 = nspecies1;
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }

                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j < nspecies1)
                                {
                                    //species[j].D = atof(temp2);
                                    species[j].D = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(species[j].D*nChange0);
                                    species[j].D = ((double) nChange1)/nChange0;
                                }
                                else
                                {
                                    cout<<"Warning! the diamter parameter: the Column number is larger than nspecies1, the input is cut from nspecies1. Please check it"<<endl;
                                    break;
                                }
                                
                                
                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        i2 = j;
                        j=0;
                        k=0;
                        delete [] temp2;
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read Diameter:"<<endl;
                    }
                    
                    
                    if(i2 < nspecies1)
                    {
                        cout<<"Warning! the external potential: i2="<<i2<<" < nspecies1. The rest external potential Parameters are automatic completion. Please check it"<<endl;
                        for(short i3=i2; i3<nspecies1; ++i3)
                        {
                            species[i3].D = species[i2-1].D;
                        }
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[9].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }

                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j==0) //sysInfo.permiS = atof(temp2);
                                {
                                    sysInfo.permiS = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(sysInfo.permiS*nChange0);
                                    sysInfo.permiS = ((double) nChange1)/nChange0;
                                }
                                if(j==1) //sysInfo.permiW = atof(temp2);
                                {
                                    sysInfo.permiW = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(sysInfo.permiW*nChange0);
                                    sysInfo.permiW = ((double) nChange1)/nChange0;
                                }
                                if(j==2) //sysInfo.T = atof(temp2);
                                {
                                    sysInfo.T = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(sysInfo.T*nChange0);
                                    sysInfo.T = ((double) nChange1)/nChange0;
                                }
                                if(j==3) sysInfo.lenth = atof(temp2);

                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        delete [] temp2;
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read PermiTempLen:"<<endl;
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[10].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }

                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j==0) iteraInfo.maxItera = atof(temp2);
                                if(j==1) iteraInfo.mixCoe = atof(temp2);
                                if(j==2) iteraInfo.cstep  = atoi(temp2);
                                if(j==3) iteraInfo.phL_L  = atof(temp2);
                                if(j==4) iteraInfo.errTol = atof(temp2);

                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        delete [] temp2;
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read Iterative:"<<endl;
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[11].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                       //strcpy(temp1,s.c_str());
                       len   = s.length();
                       filePath = "";
                       for(ilen=0;ilen<len;++ilen)
                       {
                            if(s[ilen] !=' ' && s[ilen] != '\r' &&
                               s[ilen] != '\n' && s[ilen] != '\t')
                            {
                                os<<s[ilen];
                                filePath += os.str();
                                os.str("");
                                os.clear();
                            }
                       }
                    }
                    else
                    {
                       cerr<<"Error in systemset.cpp: fail to read filePath:"<<endl;
                    }
                }
 /////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
        else
        {
            cerr<<"Error in systemset.cpp: fail to open the input file"<<endl;
            exit(0);
        }
        inFile.close();
        
        
        rhoCounP = 0.0;
        rhoCounN = 0.0;
        if(pInfo[0].rhopm > 0)
        {
            for(short i0=0;i0<nb[0];++i0)
            {
                species[i0].Z = zb[0][i0];
                rhoB[i0] = pInfo[0].rhopm*(((float) mb[0][i0])/((float) pInfo[0].mp));
                if(zb[0][i0] < 0)
                {
                    if(species[nbb+2].Z <= 0)
                    {
                        cerr<<"systemset.cpp: Counterion Valency Wrong"<<endl;
                        exit(0);
                    }
                    
                    rhoTemp = rhoB[i0]*zb[0][i0]/species[nbb+2].Z;
                    rhoCounP += fabs(rhoTemp);
                }
                
                if(zb[0][i0] > 0)
                {
                    if(species[nbb+3].Z >= 0)
                    {
                        cerr<<"systemset.cpp: Counterion Valency Wrong"<<endl;
                        exit(0);
                    }
                    
                    rhoTemp = rhoB[i0]*zb[0][i0]/species[nbb+3].Z;
                    rhoCounN += fabs(rhoTemp);
                }
            }
        }
        
        if(pInfo[1].rhopm > 0)
        {
            for(short i0=0;i0<nb[1];++i0)
            {
                species[i0+nb[0]].Z = zb[1][i0];
                rhoB[i0+nb[0]] = pInfo[1].rhopm*(((float) mb[1][i0])/((float) pInfo[1].mp));
                if(zb[1][i0] < 0)
                {
                    if(species[nbb+2].Z <= 0)
                    {
                        cerr<<"systemset.cpp: Counterion Valency Wrong"<<endl;
                        exit(0);
                    }
                    
                    rhoTemp = rhoB[i0+nb[0]]*zb[1][i0]/species[nbb+2].Z;
                    rhoCounP += fabs(rhoTemp);
                }
                
                if(zb[1][i0] > 0)
                {
                    if(species[nbb+3].Z >= 0)
                    {
                        cerr<<"systemset.cpp: Counterion Valency Wrong"<<endl;
                        exit(0);
                    }
                    
                    rhoTemp = rhoB[i0+nb[0]]*zb[1][i0]/species[nbb+3].Z;
                    rhoCounN += fabs(rhoTemp);
                }
            }
        }
        
        rhoB[nbb+2] = rhoCounP;
        rhoB[nbb+3] = rhoCounN;
        
        ZsaltPo = fabs(species[nbb].Z);
        ZsaltNe = fabs(species[nbb+1].Z);
        
        rhoB[nbb]  = rhoSalt*ZsaltNe;
        rhoB[nbb+1]= rhoSalt*ZsaltPo;
        
        species[nspecies1-1].Z = 0;
        
        
        eta  = 0;
        lunit= sysInfo.lenth;
        coe3 = 1.0E-3/(lunit*lunit*lunit)/Na; //density from N/lunit^3 to Molar
        for(short i1=0; i1<(nspecies1-1); ++i1)
        {
            
            D3     = species[i1].D*species[i1].D*species[i1].D;
            rhoB[i1]= rhoB[i1]/coe3;
            eta  = eta + rhoB[i1]*D3*Pi/6.0;
            if(rhoB[i1] < 1E-10) species[i1].Z = 0;
        }
        
        D3 = species[nspecies1-1].D*species[nspecies1-1].D*species[nspecies1-1].D;
        rhoB[nspecies1-1] = 6*(eta0 - eta)/(Pi*D3);
        
        
        iDent4 = 0;
        iDent5 = 0;
        if(species[nbb].uWall != species[nbb+2].uWall) iDent4 = 1;
        if(species[nbb+1].uWall != species[nbb+3].uWall) iDent5 = 1;
        
        for(short i1=0; i1<(nspecies1-1); ++i1)
        {
            if(pairEner[nbb][i1] != pairEner[nbb+2][i1]) iDent4 = 1;
            
            if(pairEner[nbb+1][i1] != pairEner[nbb+3][i1]) iDent5 = 1;
        }
        
        if((species[nbb].Z==species[nbb+2].Z) && (species[nbb].D==species[nbb+2].D) && iDent4 == 0)
        {
            rhoB[nbb] = rhoB[nbb] + rhoB[nbb+2];
            species[nbb+2].Z = 0;
            rhoB[nbb+2] = 0;
        }
        
        if((species[nbb+1].Z==species[nbb+3].Z) && (species[nbb+1].D==species[nbb+3].D) && iDent5 == 0)
        {
            rhoB[nbb+1] = rhoB[nbb+1] + rhoB[nbb+3];
            species[nbb+3].Z = 0;
            rhoB[nbb+3] = 0;
        }
        
        delete [] i22;
        
        //delete temp1;
    }
    
    
    
    void readSystem(PolyInfo* pInfo,char fileName[],MethodGeom& MethG)
    {
        ios::sync_with_stdio(false);
        cin.tie(0);
        string s;
        string mark[12];
        //string markPath = "../program/src/mark.dat";
        char* temp2;
        short i,j,k,n;
        short len;
        short ilen;
        int   nChange0,nChange1;
        
        ifstream inFile;
        ostringstream os;
        
        mark[0] = "METHOD:";
        mark[1] = "POLYMER:";
        mark[2] = "SEQUENCE:";
        mark[3] = "SIZE:";
        mark[4] = "SALT_HS:";
        mark[5] = "WALL:";
        mark[6] = "ENERGY:";
        mark[7] = "VALENCY:";
        mark[8] = "DIAMETER:";
        mark[9] = "PERMITEMLEN:";
        mark[10] = "ITERATIVE:";
        mark[11] = "FILEPATH:";
        
        cout<<"Please type in the path and name of the input file, e.g., ~/input/input.dat"<<endl;
        scanf("%s",fileName);
        inFile.open(fileName,ios_base::in);
        //inFile.open("input.dat",ios_base::in);
        if(inFile)
        {
            while(getline(inFile,s))
            {
                /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[0].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }
                            
                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j==0)
                                {
                                    MethG.Method = "";
                                    for(short i1=0; i1<i; ++i1)
                                    {
                                        if(isupper(temp2[i1])) temp2[i1] = tolower(temp2[i1]);
                                        os<<temp2[i1];
                                        MethG.Method += os.str();
                                        os.str("");
                                        os.clear();
                                    }
                                }
                                if(j==1)
                                {
                                    MethG.Geomet = "";
                                    for(short i1=0; i1<i; ++i1)
                                    {
                                        if(isupper(temp2[i1])) temp2[i1] = tolower(temp2[i1]);
                                        os<<temp2[i1];
                                        MethG.Geomet += os.str();
                                        os.str("");
                                        os.clear();
                                    }
                                }
                                if(j==2)
                                {
                                    MethG.nSurfs = "";
                                    for(short i1=0; i1<i; ++i1)
                                    {
                                        if(isupper(temp2[i1])) temp2[i1] = tolower(temp2[i1]);
                                        os<<temp2[i1];
                                        MethG.nSurfs += os.str();
                                        os.str("");
                                        os.clear();
                                    }
                                }
                                
                                if(j==3)
                                {
                                    MethG.exPoten = "";
                                    for(short i1=0; i1<i; ++i1)
                                    {
                                        if(isupper(temp2[i1])) temp2[i1] = tolower(temp2[i1]);
                                        os<<temp2[i1];
                                        MethG.exPoten += os.str();
                                        os.str("");
                                        os.clear();
                                    }
                                }

                                
                                if(j==4) MethG.a = atof(temp2);
                                if(j==5) MethG.b = atof(temp2);

                                
                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        
                        delete [] temp2;
                    }
                    else
                    {
                        cerr<<"Error in systemset.cpp: fail to read Method:"<<endl;
                        exit(0);
                    }
                }
                
                
                /////////////////////////////////////////////////////////////////////////////////////////////
                if(strcmp(s.c_str(),mark[1].c_str())==0)
                {
                    if(getline(inFile,s))
                    {
                        ilen = 0;
                        //strcpy(temp1,s.c_str());
                        len   = s.length();
                        temp2 = new char[len+1] ();
                        i=0;
                        j=0;
                        k=0;
                        n=1;
                        do{
                            if(k!=0)  ++ilen;
                            if(s[ilen]!=' ' && s[ilen] != ',' && s[ilen] != ';')
                            {
                                temp2[i] = s[ilen];
                                ++i;
                                n=0;
                            }
                            
                            if((s[ilen]==' ' || s[ilen] == ','|| s[ilen] == ';'
                                || k==(len-1))&&(n==0))
                            {
                                if(j==0) //pInfo[0].rhopm= atof(temp2);
                                {
                                    pInfo[0].rhopm = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(pInfo[0].rhopm*nChange0);
                                    pInfo[0].rhopm = ((double) nChange1)/nChange0;
                                }
                                if(j==1) pInfo[0].mp   = atoi(temp2);
                                if(j==2) pInfo[0].nb   = atoi(temp2);
                                if(j==3)
                                {
                                    pInfo[0].MODEL = "";
                                    for(short i1=0; i1<i; ++i1)
                                    {
                                        if(isupper(temp2[i1])) temp2[i1] = tolower(temp2[i1]);
                                        os<<temp2[i1];
                                        pInfo[0].MODEL += os.str();
                                        os.str("");
                                        os.clear();
                                    }
                                }
                                if(j==4) //pInfo[0].epsilon_b= atof(temp2);
                                {
                                    pInfo[0].epsilon_b = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(pInfo[0].epsilon_b*nChange0);
                                    pInfo[0].epsilon_b = ((double) nChange1)/nChange0;
                                }
                                if(j==5) //pInfo[1].rhopm= atof(temp2);
                                {
                                    pInfo[1].rhopm = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(pInfo[1].rhopm*nChange0);
                                    pInfo[1].rhopm = ((double) nChange1)/nChange0;
                                }
                                if(j==6) pInfo[1].mp   = atoi(temp2);
                                if(j==7) pInfo[1].nb   = atoi(temp2);
                                if(j==8)
                                {
                                    pInfo[1].MODEL = "";
                                    for(short i1=0; i1<i; ++i1)
                                    {
                                        if(isupper(temp2[i1])) temp2[i1] = tolower(temp2[i1]);
                                        os<<temp2[i1];
                                        pInfo[1].MODEL += os.str();
                                        os.str("");
                                        os.clear();
                                    }
                                }
                                if(j==9) //pInfo[1].epsilon_b= atof(temp2);
                                {
                                    pInfo[1].epsilon_b = atof(temp2);
                                    nChange0 = pow(10,i-2);
                                    if(nChange0 < 1) nChange0 = 1;
                                    nChange1 = round(pInfo[1].epsilon_b*nChange0);
                                    pInfo[1].epsilon_b = ((double) nChange1)/nChange0;
                                }
                                
                                memset(temp2,0,(i+1)*sizeof(char));
                                ++j;
                                i=0;
                                n=1;
                            }
                            ++k;
                        }while(k < len);
                        j=0;
                        k=0;
                        
                        delete [] temp2;
                    }
                    else
                    {
                        cerr<<"Error in systemset.cpp: fail to read PolymerInfo:"<<endl;
                        exit(0);
                    }
                }
                
                
                
                
            }
        }
        else
        {
            cerr<<"Error in systemset.cpp: fail to open the input file"<<endl;
            exit(0);
        }
        inFile.close();
        
    }
    
}
    
