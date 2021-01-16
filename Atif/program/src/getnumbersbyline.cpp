//*****************for get numbers line by line************//
#include "getnumbersbyline.h"

void GetNumbersByLine(ifstream& inFile,double* getNum)
{

    ios::sync_with_stdio(false);
    cin.tie(0);
    string s;
    char* temp2;
    short i,j,k,n;
    short len;
    short ilen;
 
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
                getNum[j] = atof(temp2);
                
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
        cerr<<"Error in getnumbersbyline.cpp: fail to read the numbers:"<<endl;
        exit(0);
    }
}

